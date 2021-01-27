"""Coexpression Networks."""

import logging
import psutil

import numpy as np
import pandas as pd
import minus80 as m80
import locuspocus as lp
import camoco.PCCUP as PCCUP
import markov_clustering as mc
import statsmodels.api as sm

from collections import Counter
from itertools import chain
from tinydb import where

from scipy import sparse
from scipy.special import comb
from scipy.cluster.hierarchy import leaves_list
from scipy.stats import pearsonr


from typing import Dict, Optional

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
        self._clusters = None
        self._expr_name_index = None

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
    def num_loci(self) -> int:
        """Numer of locu in the network."""
        return len(self.expr.index)

    @property
    def num_accessions(self) -> int:
        """Number of accessions in network."""
        return len(self.expr.columns)

    @property
    def expr_name_index(self) -> Dict[str,int]:
        """Retreive the expr index based on locus name."""
        if self._expr_name_index is None:
            self._expr_name_index = {
                name:index for index,name in enumerate(self.expr.index)
            }
        return self._expr_name_index


    @property
    def zscore_cutoff(self) -> float:
        return self.metadata.search(
            where('name') == "significance_threshold"
        ).pop()["val"]

    @zscore_cutoff.setter
    def zscore_cutoff(self, val: float) -> None:
        self.coex.significant = self.coex.score >= float(val)
        self.m80.col["coex"] = self.coex
        self.metadata.upsert(
            {'val': float(val)},
            ((where('name')=='significance_threshold') &
            (where('type')=='coexpression_parameter'))
        )

    @property
    def clusters(self) -> lp.Ontology:
        if self._clusters is None:
            self._clusters = lp.Ontology("MCL", rootdir=self.m80.thawed_dir)
        return self._clusters

    # -----------------------------------------
    #       Methods
    # -----------------------------------------


    def subnetwork(
        self,
        loci=None,
        sig_only=True,
        min_distance=None,
        filter_missing_loci=False,
        trans_loci_only=False,
        names_as_index=True,
        names_as_cols=False,
    ):
        """
        Extract a subnetwork of edges exclusively between loci
        within the locus_list. Also includes various options for
        what information to report, see Parameters.

        Parameters
        ----------
        loci : lp.Loci
            The loci from which to extract a subnetwork.
            If loci is None, the function will assume
            locus_list is all loci in the network (self).
            NOTE: input loci will be associated with network loci
            based entirely on *NAME*. 
        sig_only : bool
            A flag to include only significant network interactions.
        min_distance : bool (default: None)
            If not None, only include interactions that are
            between loci that are a `min_distance` away from
            one another.
        filter_missing_loci : bool (default: False)
            Filter out loci that are not in the current
            Coex object (self). 
        trans_loci_only : bool (default: True)
            Filter out locus interactions that are not in Trans,
            this argument requires that locus attr object has
            the 'parent_locus' key:val set to distinguish between
            cis and trans elements.
        names_as_index : bool (default: True)
            Include locus names as the index.
        names_as_cols : bool (default: False)
            Include locus names as two columns named 'locus_a' and 'locus_b'.

        Returns
        -------
        A pandas.DataFrame containing the edges. Columns
        include score, significant (bool), and inter-genic distance.
        """
        if loci is None:
            df = self.coex
        else:
            # Extract the indices for each locus
            expr_indices = []
            for l in loci:
                try:
                    idx = self.expr_name_index[l.name]
                except KeyError as e:
                    if filter_missing_loci:
                        continue
                    else:
                        raise ValueError("Found a locus name not in the expr matrix") from e
                else:
                    expr_indices.append(idx)
            if len(expr_indices) == 0:
                df = pd.DataFrame(columns=["score", "significant", "distance"])
            else:
                # Grab the coexpression indices for the loci
                expr_indices = np.array(sorted(expr_indices))
                indices = PCCUP.coex_index(expr_indices, self.num_loci)
                df = self.coex.iloc[indices]
        if sig_only:
            df = df[df.significant] 
        if min_distance is not None:
            df = df[df.distance >= min_distance]
        if names_as_index or names_as_cols or trans_loci_only:
            names = self.expr.index.values
            coex_indices = df.index.values
            if len(coex_indices) > 0:
                expr_indices = PCCUP.coex_expr_index(coex_indices, self.num_loci)
                df.insert(0, "locus_a", names[expr_indices[:, 0]])
                df.insert(1, "locus_b", names[expr_indices[:, 1]])
            else:
                df.insert(0, "locus_a", [])
                df.insert(0, "locus_b", [])
        if names_as_index and not names_as_cols:
            df = df.set_index(["locus_a", "locus_b"])
        if trans_loci_only:
            try:
                parents = {x.id: x.attr["parent_locus"] for x in loci}
            except KeyError as e:
                raise KeyError(
                    "Each locus must have 'parent_locus'"
                    " attr set to calculate trans only"
                ) from e
            df["trans"] = [
                parents[locus_a] != parents[locus_b]
                for locus_a, locus_b in zip(
                    df.index.get_level_values(0), df.index.get_level_values(1)
                )
            ]
        return df

    def subnetwork_degree(
        self,
        loci: lp.Loci,
        trans_loci_only: bool = False
    ):
        """
        Calculate the degree among loci in a subnetwork

        Parameters
        ----------
        loci : iterable (co.Locus object)
            a list of loci for which to retrieve local degree for
        trans_loci_only: bool = False
            If set, only include degree among trans loci. 
            this argument requires that locus attr object has
            the 'parent_locus' key:val set to distinguish between
            cis and trans elements.

        """
        subnetwork = self.subnetwork(
            loci, sig_only=True, trans_loci_only=trans_loci_only
        )
        if trans_loci_only:
            # Filter out loci from the same parent locus
            subnetwork = subnetwork.iloc[subnetwork.trans]
        local_degree = pd.DataFrame(
            list(Counter(chain(*subnetwork.index.values)).items()),
            columns=["Gene", "Degree"],
        ).set_index("Gene")
        # We need to find loci not in the subnetwork and add them as degree 0
        # The code below is ~optimized~
        # DO NOT alter unless you know what you're doing :)
        degree_zero_loci = pd.DataFrame(
            [(locus.name, 0) for locus in loci if locus.name not in local_degree.index],
            columns=["Gene", "Degree"],
        ).set_index("Gene")
        return pd.concat([local_degree, degree_zero_loci])

    def to_sparse_matrix(
        self,
        loci: Optional[lp.Loci] = None,
        min_distance: Optional[float] = None,
        max_edges: Optional[int] = None,
        remove_orphans: bool = False
    ) -> sparse.coo_matrix:
        """
        Convert the co-expression interactions to a
        scipy sparse matrix.

        Parameters
        -----
        loci: iter of Loci (default: None)
            If specified, return only the interactions among
            loci in the list. If None, use all loci in the network.
        min_distance : int (default: None)
            The minimum distance between loci for which to consider
            co-expression interactions. This is useful for  filtering
            out potential cis edges.
        max_edges: int (default: None)
            Threshold the matrix to only contain the strongest `max_edges`
            based on score.
        remove_orphans: bool = False
            Remove loci with zero significant edges

        Returns
        -------
        A tuple (a,b) where 'a' is a scipy sparse matrix and
        'b' is a mapping from locus_name to index.

        """
        log.info("Getting locuss")
        # first get the subnetwork in pair form
        log.info("Pulling edges")
        edges = self.subnetwork(
            loci=loci,
            min_distance=min_distance,
            sig_only=True,
            names_as_cols=True,
            names_as_index=False,
            filter_missing_loci=False
        )
        # NOTE: perform this after the subnetwork call as it is more efficient
        if loci is None:
            loci = list(self.loci)
        # Option to limit the number of edges
        if max_edges is not None:
            log.info("Filtering edges")
            edges = edges.sort_values(by="score", ascending=False)[
                0 : min(max_edges, len(edges))
            ]
        # Option to restrict locus list to only locuss with edges
        if remove_orphans:
            log.info("Removing orphans")
            not_orphans = set(edges.locus_a).union(edges.locus_b)
            loci = [g for g in loci if g.name in not_orphans]
        # Create a locus index
        log.info("Creating Loci name Index")
        locus_index = {l.name:i for i,l in enumerate(loci)}
        # get the expression matrix indices for all the locuss
        row = [locus_index[x] for x in edges.locus_a.values]
        col = [locus_index[x] for x in edges.locus_b.values]
        data = list(edges.score.values)
        # Make the values symmetric by doubling everything
        # Note: by nature we dont have cycles so we dont have to
        #   worry about the diagonal
        log.info("Making matrix symmetric")
        d = data + data
        r = row + col
        c = col + row
        log.info("Creating matrix")
        nloci = len(locus_index)
        matrix = sparse.coo_matrix((d, (r, c)), shape=(nloci, nloci), dtype=None)
        return (matrix, locus_index)

    def mcl(
        self,
        loci: Optional[lp.Loci] = None,
        I: float = 2.0,
        min_distance: Optional[int] = None,
        min_cluster_size: Optional[int] = 0,
        max_cluster_size: Optional[int] = 10e10,
    ):
        """
        Returns clusters (as list) as designated by MCL (Markov Clustering).

        Parameters
        ----------
        loci : lp.Loci
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
        matrix, locus_index = self.to_sparse_matrix(loci=loci)
        # Run MCL
        result = mc.run_mcl(matrix, inflation=I, verbose=False)
        clusters = mc.get_clusters(result)
        # MCL traditionally returns clusters by size with 0 being the largest
        clusters = sorted(clusters, key=lambda x: len(x), reverse=True)
        # Create a dictionary to map ids to locus names
        locus_id_index = {v: k for k, v in locus_index.items()}
        result = []
        for c in clusters:
            if len(c) < min_cluster_size or len(c) > max_cluster_size:
                continue
            # convert to loci
            loci = [self.loci[locus_id_index[i]] for i in c]
            result.append(loci)
        return result

    def coexpression(self,
        locus_a: lp.Locus,
        locus_b: lp.Locus,
        recalculate: bool = False
    ) -> pd.Series:
        """
        Returns a coexpression z-score between two loci. This
        is the pearson correlation coefficient of the two locus'
        expression profiles across the accessions (experiments).

        Parameters
        ----------
        locus_a : lp.Locus
            The first locus
        locus_b : lp.Locus
            The second locus
        recalculate: bool
            If True, do not extract coex score from subnetwork, perform 
            the calculation directly.

        Returns
        -------
            A pandas series containing: Coexpression Z-Score, significant, distance

        """
        if locus_a.name == locus_b.name:
            raise ValueError("locus_a and locus_b have the same locus.name")
        if recalculate: 
            score = self._coex_concordance(locus_a, locus_b)
            significant = True if score >= 3.0 else False
            distance = locus_a.distance(locus_b)
            result = pd.Series(
                [score, significant, distance],
                name=(locus_a.name, locus_b.name),
                index=["score", "significant", "distance"],
            )
        else:
            # Pull results from subnetwork
            result = self.subnetwork([locus_a, locus_b], sig_only=False).iloc[0]
        return result

    def density(self, loci, min_distance=None, by_locus=False, filter_missing_loci=False):
        """
        Calculates the density of the non-thresholded network edges
        amongst loci within locus_list. Includes parameters to perform
        measurements for loci within a certain distance of each other.
        This corrects for cis regulatory elements increasing noise
        in coexpression network.

        Parameters
        ----------
        locus_list : iter of Loci
            List of loci from which to calculate density.
        min_distance : int (default: None)
            Ignore edges between loci less than min_distance
            in density calculation.
        by_locus : bool (default: False)
            Return a per-locus breakdown of density within the subnetwork.
        filter_missing_loci : bool (default: True)
            Filter out loci that are not in the current
            Coex object (self). 

        Returns
        -------
        A network density OR density on a locus-wise basis
        """
        # filter for only loci within network
        edges = self.subnetwork(
            loci, min_distance=min_distance, 
            sig_only=False,filter_missing_loci=filter_missing_loci
        )

        if by_locus == True:
            x = pd.DataFrame.from_records(
                chain(*[
                    ((locus_a, score), (locus_b, score))
                    for locus_a, locus_b, score, sig, dis in edges.reset_index().values
                ]),
                columns=["locus", "score"],
            )
            return x.groupby("locus").agg(np.mean)
        else:
            if len(edges) == 0:
                return np.nan
            if len(edges) == 1:
                return edges.score[0]
            return np.nanmean(edges.score) / (1 / np.sqrt(len(edges)))

    def locality(
        self,
        loci,
        iter_name=None,
        include_regression=False
    ):
        """
        Computes the merged local vs global degree table

        Parameters
        ----------
        loci : iterable of camoco.Loci
            A list or equivalent of loci
        iter_name : object (default: none)
            This will be added as a column. Useful for
            locusrating bootstraps of locality and keeping
            track of which one a row came from after catting
            multiple bootstraps together.
        include_regression : bool (default: False)
            Include the OLS regression residuals and fitted values
            on local ~ global.

        Returns
        -------
        A pandas DataFrame with local, global and residual columns
        based on linear regression of local on global degree.

        """
        global_degree = self.degree.loc[[x.name for x in loci]].fillna(0)
        local_degree = self.subnetwork_degree(loci)
        degree = global_degree.merge(local_degree, left_index=True, right_index=True)
        degree.columns = ["global", "local"]
        degree = degree.sort_values(by="global")
        if include_regression:
            # set up variables to use astype to aviod pandas sm.OLS error
            loc_deg = degree["local"]
            glob_deg = degree["global"]
            ols = sm.OLS(loc_deg.astype(float), glob_deg.astype(float)).fit()
            degree["resid"] = ols.resid
            degree["fitted"] = ols.fittedvalues
            degree = degree.sort_values(by="resid", ascending=False)
        if iter_name is not None:
            degree["iter_name"] = iter_name
        return degree

    # -----------------------------------------
    #       Internal Methods
    # -----------------------------------------

    def _coex_concordance(self, locus_a, locus_b, maxnan=10, return_dict=False):
        """
        This is a sanity method to ensure that the pcc calculated
        directly from the expr profiles matches the one stored in
        the coex vector
        """
        expr_a = self.expr.loc[locus_a.name].values
        expr_b = self.expr.loc[locus_b.name].values
        mask = np.logical_and(np.isfinite(expr_a), np.isfinite(expr_b))
        if sum(mask) < maxnan:
            # too many nans to reliably calculate pcc
            return np.nan
        r = pearsonr(expr_a[mask], expr_b[mask])[0]
        # fisher transform it
        z = np.arctanh(r - 0.0000001)
        # standard normalize it
        pcc_mean = self.metadata.search(where("name") == "pcc_mean").pop()["val"]
        pcc_std = self.metadata.search(where("name") == "pcc_std").pop()["val"]
        z = (z - pcc_mean) / pcc_std
        if return_dict:
            return {"pearsonr": r, "zscore": z}
        else:
            return z

    def _calculate_coexpression(self, significance_threshold: float = 3.0):
        """
        Generates pairwise PCCs for locus expression profiles in self.expr.
        Also calculates pairwise locus distance.
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
        # This affects the mean and std fro the locus.
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

        log.info("Calculating pairwise locus distances")
        # Calculate the pairwise locus distances
        loci_positions = self.loci.m80.db.query("SELECT chromosome,start,end FROM loci")
        # Make sure types are as expected
        distances = np.absolute(PCCUP.pairwise_locus_distances(
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
            # Translate the expr indexes to the locus names
            for i, node_degree in sigs:
                degree.loc[names[i]] = node_degree
        self.m80.col["degree"] = degree

    def _calculate_locus_hierarchy(
        self, method: str = "single"
    ) -> np.ndarray:
        """
        Calculate the hierarchical locus distance for the Expr matrix
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
        locus_link = fastcluster.linkage(dists, method=method)
        return locus_link

    def _calculate_leaves(self, method: str = "single") -> None:
        """
        This calculates the leaves of the dendrogram from the coex
        """
        locus_link = self._calculate_locus_hierarchy(method=method)
        log.info("Finding the leaves")
        leaves = leaves_list(locus_link)

        # Put them in a dataframe and stow them
        leaves = pd.DataFrame(leaves, index=self.expr.index, columns=["index"])
        self.m80.col["leaves"] = leaves

    def _calculate_clusters(self) -> None:
        """
        Calculates global clusters
        """
        terms = [lp.Term(f"MCL_{i}", loci=loci) for i,loci in enumerate(self.mcl())]

        log.info("Creating Cluster Ontology")
        self.MCL = lp.Ontology.from_terms(
            "MCL",
            terms,
            rootdir=self.m80.thawed_dir
        )
        log.info("Finished finding clusters")

    # -----------------------------------------
    #       Factory Methods
    # -----------------------------------------

    @classmethod
    def from_DataFrame(
        cls,
        name: str,
        df: pd.DataFrame,
        loci: lp.Loci,
        /,
        min_expr: float = 0.001,
        max_locus_missing_data: float = 0.2,
        max_accession_missing_data: float = 0.3,
        min_single_accession_expr: float = 0.08,
        normalize_expr_values: bool = True,
        significance: float = 3.0,
        force: bool = False,
        rootdir: str = None,
    ) -> "Coex":
        """
        Create a Coex object from a Pandas DataFrame and
        a LocusPocus Loci object.

        Required Parameters
        -------------------
        name : str
            A name for the Coex object
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

        try:
            # Create an empty class
            self = cls(name)

            # Store the metadata used to build the network
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

            # include only loci that are NAMED in the refgen
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

            # Update the database ------------------------------------------------------------------------------------------
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

        except Exception as e:
            m80.delete("Coex", name, rootdir=rootdir)
            raise e

        # Return an instance
        return self
