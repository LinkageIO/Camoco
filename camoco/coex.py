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
from scipy.cluster.hierarchy import leaves_list

from .exceptions import (
    CoexExistsError,
)

log = logging.getLogger(__name__)


class Coex(m80.Freezable):
    """
    A Co-expression Network.
    """

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
    def num_loci(self):
        return len(self.expr.index)

    @property
    def num_accessions(self):
        return len(self.expr.columns)

    # -----------------------------------------
    #       Methods
    # -----------------------------------------

    def subnetwork(
        self,
        loci=None,
        sig_only=True,
        min_distance=None,
        filter_missing_gene_ids=True,
        trans_locus_only=False,
        names_as_index=True,
        names_as_cols=False,
    ):
        """
        Extract a subnetwork of edges exclusively between loci
        within the gene_list. Also includes various options for
        what information to report, see Parameters.

        Parameters
        ----------
        gene_list : lp.Loci
            The loci from which to extract a subnetwork.
            If gene_list is None, the function will assume
            gene_list is all loci in COB object (self).
        sig_only : bool
            A flag to include only significant interactions.
        min_distance : bool (default: None)
            If not None, only include interactions that are
            between loci that are a `min_distance` away from
            one another.
        filter_missing_gene_ids : bool (default: True)
            Filter out gene ids that are not in the current
            COB object (self).
        trans_locus_only : bool (default: True)
            Filter out gene interactions that are not in Trans,
            this argument requires that locus attr object has
            the 'parent_locus' key:val set to distinguish between
            cis and trans elements.
        names_as_index : bool (default: True)
            Include gene names as the index.
        names_as_cols : bool (default: False)
            Include gene names as two columns named 'gene_a' and 'gene_b'.

        Returns
        -------
        A pandas.DataFrame containing the edges. Columns
        include score, significant (bool), and inter-genic distance.
        """
        if loci is None:
            if sig_only:
                df = self.coex[self.coex.significant]
            else:
                df = self.coex
        else:
            # Extract the ids for each Gene
            gene_list = set(sorted(gene_list))
            ids = np.array([self._expr_index[x.id] for x in gene_list])
            if filter_missing_gene_ids:
                # filter out the Nones
                ids = np.array([x for x in ids if x is not None])
            if len(ids) == 0:
                df = pd.DataFrame(columns=["score", "significant", "distance"])
            else:
                # Grab the coexpression indices for the genes
                ids = PCCUP.coex_index(ids, num_genes)
                df = self._coex_DataFrame(ids=ids, sig_only=sig_only)
        if min_distance is not None:
            df = df[df.distance >= min_distance]
        if names_as_index or names_as_cols or trans_locus_only:
            names = self._expr.index.values
            ids = df.index.values
            if len(ids) > 0:
                ids = PCCUP.coex_expr_index(ids, num_genes)
                df.insert(0, "gene_a", names[ids[:, 0]])
                df.insert(1, "gene_b", names[ids[:, 1]])
                del ids
                del names
            else:
                df.insert(0, "gene_a", [])
                df.insert(0, "gene_b", [])
        if names_as_index and not names_as_cols:
            df = df.set_index(["gene_a", "gene_b"])
        if trans_locus_only:
            try:
                parents = {x.id: x.attr["parent_locus"] for x in gene_list}
            except KeyError as e:
                raise KeyError(
                    "Each locus must have 'parent_locus'"
                    " attr set to calculate trans only"
                )
            df["trans"] = [
                parents[gene_a] != parents[gene_b]
                for gene_a, gene_b in zip(
                    df.index.get_level_values(0), df.index.get_level_values(1)
                )
            ]
        return df


    def to_sparse_matrix(
        self, gene_list=None, min_distance=None, max_edges=None, remove_orphans=False
    ):
        """
        Convert the co-expression interactions to a
        scipy sparse matrix.

        Parameters
        -----
        gene_list: iter of Loci (default: None)
            If specified, return only the interactions among
            loci in the list. If None, use all genes.
        min_distance : int (default: None)
            The minimum distance between genes for which to consider
            co-expression interactions. This filters out cis edges.

        Returns
        -------
        A tuple (a,b) where 'a' is a scipy sparse matrix and
        'b' is a mapping from gene_id to index.
        """
        from scipy import sparse

        log.info("Getting genes")
        # first get the subnetwork in pair form
        log.info("Pulling edges")
        edges = self.subnetwork(
            gene_list=gene_list,
            min_distance=min_distance,
            sig_only=True,
            names_as_cols=True,
            names_as_index=False,
        )
        # Option to limit the number of edges
        if max_edges is not None:
            log.info("Filtering edges")
            edges = edges.sort_values(by="score", ascending=False)[
                0 : min(max_edges, len(edges))
            ]
        # Option to restrict gene list to only genes with edges
        if remove_orphans:
            log.info("Removing orphans")
            not_orphans = set(edges.gene_a).union(edges.gene_b)
            gene_list = [g for g in gene_list if g.id in not_orphans]
        # Create a gene index
        log.info("Creating Index")
        if gene_list == None:
            gene_list = list(self.refgen.iter_genes())
        else:
            gene_list = set(gene_list)
        gene_index = {g.id: i for i, g in enumerate(gene_list)}
        nlen = len(gene_list)
        # get the expression matrix indices for all the genes
        row = [gene_index[x] for x in edges.gene_a.values]
        col = [gene_index[x] for x in edges.gene_b.values]
        data = list(edges.score.values)
        # Make the values symmetric by doubling everything
        # Note: by nature we dont have cycles so we dont have to
        #   worry about the diagonal
        log.info("Making matrix symmetric")
        d = data + data
        r = row + col
        c = col + row
        log.info("Creating matrix")
        matrix = sparse.coo_matrix((d, (r, c)), shape=(nlen, nlen), dtype=None)
        return (matrix, gene_index)

    def mcl(
        self,
        gene_list=None,
        I=2.0,
        min_distance=None,
        min_cluster_size=0,
        max_cluster_size=10e10,
    ):
        """
        Returns clusters (as list) as designated by MCL (Markov Clustering).

        Parameters
        ----------
        gene_list : a gene iterable
            These are the loci which will be clustered
        I : float (default: 2.0)
            This is the inflation parameter passed into mcl.
        min_distance : int (default: None)
            The minimum distance between loci for which to consider
            co-expression interactions. This filters out cis edges.
        min_cluster_size : int (default: 0)
            The minimum cluster size to return. Filter out clusters smaller
            than this.
        max_cluster_size : float (default: 10e10)
            The maximum cluster size to return. Filter out clusters larger
            than this.

        Returns
        -------
        A list clusters containing a lists of loci within each cluster
        """
        import markov_clustering as mc

        matrix, gene_index = self.to_sparse_matrix(gene_list=gene_list)
        # Run MCL
        result = mc.run_mcl(matrix, inflation=I, verbose=True)
        clusters = mc.get_clusters(result)
        # MCL traditionally returns clusters by size with 0 being the largest
        clusters = sorted(clusters, key=lambda x: len(x), reverse=True)
        # Create a dictionary to map ids to gene names
        gene_id_index = {v: k for k, v in gene_index.items()}
        result = []
        for c in clusters:
            if len(c) < min_cluster_size or len(c) > max_cluster_size:
                continue
            # convert to loci
            loci = self.refgen.from_ids([gene_id_index[i] for i in c])
            result.append(loci)
        return result

    # -----------------------------------------
    #       Internal Methods
    # -----------------------------------------

    def _calculate_coexpression(self, significance_threshold: float = 3.0):
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
        # completely for the two loci expression data making its PCC a nan.
        # This affects the mean and std fro the gene.
        pcc_mean = np.ma.masked_array(pccs, np.isnan(pccs)).mean()
        self.metadata.insert(
            {"name": "pcc_mean", "val": pcc_mean, "type": "coexpression_parameter"}
        )
        pcc_std = np.ma.masked_array(pccs, np.isnan(pccs)).std()
        self.metadata.insert(
            {"name": "pcc_std", "val": pcc_std, "type": "coexpression_parameter"}
        )

        # Calculate Z Scores -------------------------------------------------------------------------------------------
        log.info("Finding adjusted scores")
        pccs = (pccs - pcc_mean) / pcc_std
        # 3. Build the dataframe
        log.info("Build the dataframe and set the significance threshold")
        self.metadata.insert(
            {
                "name": "significance_threshold",
                "val": significance_threshold,
                "type": "coexpression_parameter",
            }
        )
        significant = pccs >= significance_threshold

        log.info("Calculating pairwise gene distances")
        # Calculate the pairwise gene distances
        loci_positions = self.loci.m80.db.query("SELECT chromosome,start,end FROM loci")
        # Make sure types are as expected
        distances = np.absolute(PCCUP.pairwise_gene_distances(
            np.ascontiguousarray(loci_positions.chromosome.values),
            np.ascontiguousarray(loci_positions.start.values.astype(np.long)),
            np.ascontiguousarray(loci_positions.end.values.astype(np.long))
        ))

        # Save data ----------------------------------------------------------------------------------------------------
        self.m80.col["coex"] = pd.DataFrame({
            "score": pccs,
            "significant": significant, 
            "distance": distances
        })

    def _calculate_degree(self):
        log.info("Building Degree")
        # Generate a df that starts all loci at 0
        names = self.expr.index.values
        degree = pd.DataFrame(0, index=names, columns=["Degree"])
        # Get the index and find the counts
        log.info("Calculating Gene degree")
        # Use the Cython function to convert the coex array index to the expr matrix index
        sigs = PCCUP.coex_expr_index(
            self.coex[self.coex.significant].index.values, self.num_loci
        )
        sigs = list(Counter(chain(*sigs)).items())
        if len(sigs) > 0:
            # Translate the expr indexes to the gene names
            for i, node_degree in sigs:
                degree.loc[names[i]] = node_degree
        self.m80.col["degree"] = degree

    def _calculate_gene_hierarchy(
        self, method: str = "single"
    ) -> np.ndarray:
        """
        Calculate the hierarchical gene distance for the Expr matrix
        using the coex data.

        Parameters
        ----------
        method: a string that is passed to fastcluster used to calculate linkage metric

        """
        import fastcluster

        # We need to recreate the original PCCs
        log.info("Calculating hierarchical clustering using {}".format(method))
        if len(self.coex) == 0:
            raise ValueError("Cannot calculate leaves without coex")
        pcc_mean = self.metadata.search(where("name") == "pcc_mean").pop()["val"]
        pcc_std = self.metadata.search(where("name") == "pcc_std").pop()["val"]
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

    def _calculate_leaves(self, method: str = "single") -> None:
        """
        This calculates the leaves of the dendrogram from the coex
        """
        gene_link = self._calculate_gene_hierarchy(method=method)
        log.info("Finding the leaves")
        leaves = leaves_list(gene_link)

        # Put them in a dataframe and stow them
        leaves = pd.DataFrame(leaves, index=self.expr.index, columns=["index"])
        self.m80.col["leaves"] = leaves

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
            "key": "qc_parameters",
            "normalize_expr_values": normalize_expr_values,
            "min_expr": min_expr,
            "max_locus_missing_data": max_locus_missing_data,
            "max_accession_missing_data": max_accession_missing_data,
            "min_single_accession_expr": min_single_accession_expr,
        }

        for k, v in metadata.items():
            self.metadata.insert(
                {
                    "name": k,
                    "val": v,
                    "type": "qc_parameter",
                }
            )

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
            )
        )

        # Filter out loci with do not meet a minimum expr threshold in at least 1 accession
        qc_loci["pass_min_expression"] = df.apply(
            lambda x: any(x >= min_single_accession_expr), axis=1  # 1 is column
        )
        log.info(
            "Found {} loci which "
            "do not have one accession above {}".format(
                sum(qc_loci["pass_min_expression"] == False), min_single_accession_expr
            )
        )

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
            )
        )
        # Update the total QC passing column
        qc_accession["PASS_ALL"] = qc_accession.apply(lambda row: np.all(row), axis=1)
        df = df.loc[:, qc_accession["PASS_ALL"]]

        # Update the database
        log.info("Filtering loci to only those passing QC")
        lp.Loci.from_loci(
            self.m80.name,
            [loci[x] for x in df.index if x in loci],
            rootdir=self.m80.thawed_dir,
        )
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
            )
        )
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
