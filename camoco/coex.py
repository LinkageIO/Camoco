"""Coexpression Networks."""

import logging
import psutil

import numpy as np
import pandas as pd
import minus80 as m80
import locuspocus as lp
import camoco.PCCUP as PCCUP

from collections import Counter
from itertools import chain
from scipy.special import comb
from tinydb import where
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram

from .exceptions import (
    CoexExistsError,
)

log = logging.getLogger(__name__)


class Coex(m80.Freezable):

    def __init__(self, name: str, rootdir: str = None) -> None:
        # Init the minus80 freezer
        super().__init__(name, rootdir=rootdir)

        # cold cache
        self._loci = None
        self._expr = None
        self._coex = None
        self._degree = None

        self.metadata = self.m80.doc.table("metadata")

    # -----------------------------------------
    #       Properties
    # -----------------------------------------

    @property
    def loci(self) -> lp.Loci:
        # lazy evaluation
        if self._loci is None:
            self._loci = lp.Loci(self.m80.name, rootdir=self.m80.thawed_dir)
        return self._loci

    @property
    def expr(self) -> pd.DataFrame:
        if self._expr is None:
            self._expr = self.m80.col["expr"]
        return self._expr

    @property
    def coex(self) -> pd.DataFrame:
        if self._coex is None:
            self._coex = self.m80.col["coex"]
        return self._coex

    @property
    def degree(self) -> pd.DataFrame:
        if self._degree is None:
            self._degree = self.m80.col["degree"]
        return self._degree

    @property
    def num_genes(self):
        return len(self.expr.index)

    @property
    def num_accessions(self):
        return len(self.expr.columns)


    # -----------------------------------------
    #       Methods
    # -----------------------------------------


    # -----------------------------------------
    #       Internal Methods
    # -----------------------------------------

    def _calculate_coexpression(self,significance_threshold=3.0):
        """
            Generates pairwise PCCs for gene expression profiles in self.expr.
            Also calculates pairwise gene distance.
        """
        # Calculate the PCCs -------------------------------------------------------------------------------------------
        log.info("Calculating Coexpression")

        # Estimate if we have enough RAM
        num_bytes_needed = comb(len(self.expr), 2) * 8
        if num_bytes_needed > psutil.virtual_memory().available:
            raise MemoryError("Not enough RAM to calculate co-expression network")
        # pass in a contigious array to the cython function to calculate PCCs
        pccs = PCCUP.pair_correlation(
            np.ascontiguousarray(
                # PCCUP expects floats
                self.expr.to_numpy()
            )
        )

        log.info("Applying Fisher Transform")
        pccs[pccs >= 1.0] = 0.9999999
        pccs[pccs <= -1.0] = -0.9999999
        pccs = np.arctanh(pccs)

        # Do a PCC check to make sure they are not all NaNs
        if not any(np.logical_not(np.isnan(pccs))):
            raise ValueError(
                "Not enough data is available to reliably calculate co-expression, "
                "please ensure you have more than 10 accessions to calculate correlation coefficient"
            )

        log.info("Calculating Mean and STD")
        # Sometimes, with certain datasets, the NaN mask overlap
        # completely for the two genes expression data making its PCC a nan.
        # This affects the mean and std fro the gene.
        pcc_mean = np.ma.masked_array(pccs, np.isnan(pccs)).mean()
        self.metadata.insert({
            'name': 'pcc_mean',
            'val':  pcc_mean,
            'type':'coexpression_parameter'
        })
        pcc_std = np.ma.masked_array(pccs, np.isnan(pccs)).std()
        self.metadata.insert({
            'name': 'pcc_std',
            'val':  pcc_std,
            'type':'coexpression_parameter'
        })

        # Calculate Z Scores -------------------------------------------------------------------------------------------
        log.info("Finding adjusted scores")
        pccs = (pccs - pcc_mean) / pcc_std
        # 3. Build the dataframe
        log.info("Build the dataframe and set the significance threshold")
        self.metadata.insert({
            'name': "significance_threshold", 
            'val':  significance_threshold,
            'type':'coexpression_parameter'
        })
        significant = pccs >= significance_threshold 

        # Save data ----------------------------------------------------------------------------------------------------
        self.m80.col["coex"] = pd.DataFrame({"score":pccs, "significant": significant})
        log.info("Done")

    def _calculate_degree(self):
        log.info("Building Degree")
        # Generate a df that starts all genes at 0
        names = self.expr.index.values
        degree = pd.DataFrame(0, index=names, columns=["Degree"])
        # Get the index and find the counts
        log.info("Calculating Gene degree")
        # Use the Cython function to convert the coex array index to the expr matrix index
        sigs = PCCUP.coex_expr_index(self.coex[self.coex.significant].index.values, self.num_genes)
        sigs = list(Counter(chain(*sigs)).items())
        if len(sigs) > 0:
            # Translate the expr indexes to the gene names
            for i, node_degree in sigs:
                degree.loc[names[i]] = node_degree
        self.m80.col["degree"] = degree


    def _calculate_gene_hierarchy(self, method='single'):
        """
            Calculate the hierarchical gene distance for the Expr matrix
            using the coex data.

            Notes
            -----
            This is kind of expenive.
        """
        import fastcluster

        # We need to recreate the original PCCs
        log.info("Calculating hierarchical clustering using {}".format(method))
        if len(self.coex) == 0:
            raise ValueError("Cannot calculate leaves without coex")
        pcc_mean = self.metadata.search(where('name')=='pcc_mean').pop()['val']
        pcc_std = self.metadata.search(where('name')=='pcc_std').pop()['val']
        # Get score column and dump coex from memory for time being
        dists = self.coex.score
        # Subtract pccs from 1 so we do not get negative distances
        dists = (dists * pcc_std) + pcc_mean
        dists = np.tanh(dists)
        dists = 1 - dists
        # convert nan to 0's, linkage can only use finite values
        dists[np.isnan(dists)] = 0
        # Find the leaves from hierarchical clustering
        gene_link = fastcluster.linkage(dists, method=method)
        return gene_link

    def _calculate_leaves(self):
        """
            This calculates the leaves of the dendrogram from the coex
        """
        gene_link = self._calculate_gene_hierarchy(method=method)
        self.log("Finding the leaves")
        leaves = leaves_list(gene_link)

        # Put them in a dataframe and stow them
        self.leaves = pd.DataFrame(leaves, index=self._expr.index, columns=["index"])
        self._gene_link = gene_link
        self._bcolz("leaves", df=self.leaves)

        # Cleanup and reinstate the coex table
        return self

        pass

    def _calculate_clusters(self):
        pass

    # -----------------------------------------
    #       Factory Methods
    # -----------------------------------------

    @classmethod
    def from_DataFrame(
        cls,
        name: str,
        description: str,
        df: pd.DataFrame,
        loci: lp.Loci,
        rootdir: str = None,
        /,
        min_expr: float = 0.001,
        max_locus_missing_data: float = 0.2,
        max_accession_missing_data: float = 0.3,
        min_single_accession_expr: float = 0.08,
        normalize_expr_values: bool = True,
        significance: float = 3.0,
        force: bool = False,
    ) -> "Coex":
        """
        Create a Coex object from a Pandas DataFrame and
        a LocusPocus Loci object.

        Required Parameters
        -------------------
        name : str
            A name for the Coex object
        description : str
            A short description for the Coex object
        df : pandas.DataFrame
            A DataFrame containing expression values. Rows are Loci and
            columns are Accessions.
        loci : locuspocus.Loci
            A set of reference loci for each row in the df. The mapping
            here is the index values in the df correspond to top-level
            names (IDs) in the Loci object. NOTE: df index values without a
            corresponding Locus in the Loci object will be filtered.

        Optional Parameters
        -------------------
        min_expr : int (default: 0.001)
            FPKM (or equivalent) values under this threshold will be set to
            NaN and not used during correlation calculations.
        max_locus_missing_data : float (default: 0.2)
            Maximum percentage missing data a locus can have. Loci  with more
            than the amount of missing data are removed.
        max_accession_missing_data : float (default: 0.5)
            Maximum percentage missing data an accession (experiment) can
            have before it is removed.
        min_single_accession_expr : float (default: 0.08)
            Loci that do not have a single accession having an expression
            value above this threshold are removed from analysis. These are
            likely presence/absence and will not have a strong coexpression
            pattern.
        normalize_expr_values : bool (default: True)
            Normalize the expr values using np.arcsinh. This compresses
            larger values more than smaller values similar to a log transformation.
            Inverse Hyperbolic Sine is recommended for "count" type data from RNA-Seq.
        force : bool
            If true and a Coex already exists with the same name, it will
            be overwritten.
        """
        if force:
            m80.delete("Coex", name)
        # Sanity check
        if m80.exists("Coex", name):
            raise CoexExistsError(f"{name} already exists.")
        # Create an empty class
        self = cls(name)

        metadata = {
            'key': 'qc_parameters',
            'normalize_expr_values': normalize_expr_values,
            'min_expr': min_expr,
            'max_locus_missing_data': max_locus_missing_data,
            'max_accession_missing_data': max_accession_missing_data,
            'min_single_accession_expr': min_single_accession_expr,
        }
       
        for k,v in metadata.items():
            self.metadata.insert({
                'name':k,
                'val':v,
                'type':'qc_parameter',
            })

        # Make a copy of the input df
        df = df.copy()

        # Set the raw data to be the original df
        self.m80.col["raw_expr"] = df

        # Perform Quality Control on loci -----------------------------------------------------------------------------

        # Create a bool df for each: row, col to store if each passes QC
        qc_loci = pd.DataFrame({"has_id": True}, index=df.index)
        qc_accession = pd.DataFrame({"has_id": True}, index=df.columns)

        # include only loci that are in the refgen
        cur = loci.m80.db.cursor()
        qc_loci["pass_in_loci"] = np.array(
            [
                bool(x)
                for x, in cur.executemany(
                    "SELECT EXISTS ( SELECT 1 FROM loci WHERE name = ?)",
                    ((name,) for name in df.index),
                ).fetchall()
            ]
        )
        log.info(
            f"Found {sum(qc_loci['pass_in_loci'] == False)} df row IDs not in {loci.m80.name}"
        )

        # Set minimum expr threshold
        df[df < min_expr] = np.nan

        # Filter out loci with too much missing data
        qc_loci["pass_missing_data"] = df.apply(
            lambda x: ((sum(np.isnan(x))) < len(x) * max_locus_missing_data), axis=1
        )
        log.info(
            "Found {} loci with > {} missing data".format(
            sum(qc_loci["pass_missing_data"] == False),
            max_locus_missing_data,
        ))

        # Filter out loci with do not meet a minimum expr threshold in at least 1 accession
        qc_loci["pass_min_expression"] = df.apply(
            lambda x: any(x >= min_single_accession_expr), axis=1  # 1 is column
        )
        log.info(
            "Found {} loci which " "do not have one accession above {}".format(
            sum(qc_loci["pass_min_expression"] == False),
            min_single_accession_expr
        ))

        qc_loci["PASS_ALL"] = qc_loci.apply(lambda row: np.all(row), axis=1)
        df = df.loc[qc_loci["PASS_ALL"], :]

        # Perform quality control on accessions ------------------------------------------------------------------------
        # Filter out ACCESSIONS with too much missing data
        qc_accession["pass_missing_data"] = df.apply(
            lambda col: (
                ((sum(np.isnan(col)) / len(col)) <= max_accession_missing_data)
            ),
            axis=0,  # 0 is columns
        )
        log.info(
            "Found {} accessions with > {} missing data".format(
                sum(qc_accession["pass_missing_data"] == False),
                max_accession_missing_data,
        ))
        # Update the total QC passing column
        qc_accession["PASS_ALL"] = qc_accession.apply(lambda row: np.all(row), axis=1)
        df = df.loc[:, qc_accession["PASS_ALL"]]
       
        # Update the database
        log.info("Filtering loci to only those passing QC")
        lp.Loci.from_loci(self.m80.name,[loci[x] for x in df.index if x in loci], rootdir=self.m80.thawed_dir)
        self.m80.col["qc_accession"] = qc_accession
        self.m80.col["qc_loci"] = qc_loci
        self.m80.col["expr"] = df


        # Do some reporting --------------------------------------------------------------------------------------------
        log.info(f"Loci passing QC:\n{str(qc_loci.apply(sum, axis=0))}")
        log.info(f"Accessions passing QC:\n{str(qc_accession.apply(sum, axis=0))}")
        # Also report a breakdown by chromosome
        qc_loci = qc_loci[qc_loci["pass_in_loci"]]
        qc_loci["chromosome"] = [loci[x].chromosome for x in qc_loci.index]
        log.info(
            "Loci passing QC by chromosome:\n{}".format(
                 str(qc_loci.groupby("chromosome").aggregate(sum)),
        ))
        # update the df to reflect only locis/accession passing QC
        log.info("Kept: {} locis {} accessions".format(len(df.index), len(df.columns)))

        if normalize_expr_values:
            # Apply inverse hyperbolic sine
            # per https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3634062/
            df = df.apply(np.arcsinh, axis=0)
        # Remove DF entries that are not in Loci

        self._calculate_coexpression()
        self._calculate_degree()
        self._calculate_leaves()
        self._calculate_clusters()

        # Return an instance
        return self

    # -----------------------------------------
    #       Static Methods
    # -----------------------------------------
