import logging
import warnings

import numpy as np
import pandas as pd
import minus80 as m80
import locuspocus as lp

from typing import List

from .exceptions import (
    CoexEmptyError, CoexExistsError   
)

log = logging.getLogger(__name__)

class Coex(m80.Freezable):

    def __init__(
        self, 
        name:str, 
        basedir: str = None
    ) -> None:
        # Init the minus80 freezer
        super().__init__(name, basedir=basedir)
        # cold cache
        self._loci = None
        self._expr = None
        
    #-----------------------------------------
    #       Attributes
    #-----------------------------------------

    @property
    def loci(self) -> lp.Loci:
        # lazy evaluation
        if self._loci is None:
            self._loci = lp.Loci(self.m80.name,basedir=self.m80.thawed_dir)
        return self._loci

    @property
    def _expr(self) -> pd.DataFrame:
        if self.__expr is None:
            self.__expr = self.m80.col['expr']
        return self.__expr

    @_expr.setter
    def _expr(self,val):
        self.__expr = val

    #-----------------------------------------
    #       Methods
    #-----------------------------------------

    def expr(
        loci: lp.Loci = None,
        accessions: List[str] = None,
        raw: bool = False
        # gene_normalize=False
    ) -> pd.DataFrame:
        # fetch the correct expr matrix
        if raw:
            df = self.m80.col['raw_expr'] 
        else:
            df = self._expr

        return df
        
    #-----------------------------------------
    #       Factory Methods
    #-----------------------------------------
    
    @classmethod
    def from_DataFrame(
        cls, 
        name: str,
        description: str,
        df: pd.DataFrame, 
        loci: lp.Loci,
        /,
        normalize_expr_values: bool = True,
        min_expr: float = 0.001,
        max_locus_missing_data: float = 0.2,
        max_accession_missing_data: float = 0.3,
        min_single_accession_expr: float = 0.08,
        force: bool = False
    ) -> "Coex":
        '''
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
            normalize_expr_values : bool (default: True)
                Normalize the expr values using np.arcsinh. This compresses
                larger values more than smaller values similar to a log transformation.
                Inverse Hyperbolic Sine is recommended for "count" type data from RNA-Seq.
            min_expr : int (default: 0.001)
                FPKM (or equivalent) values under this threshold will be set to
                NaN and not used during correlation calculations.
            max_gene_missing_data : float (default: 0.2)
                Maximum percentage missing data a gene can have. Genes with more
                than the amount of missing data are removed.
            max_accession_missing_data : float (default: 0.5)
                Maximum percentage missing data an accession (experiment) can
                have before it is removed.
            min_single_sample_expr : float (default: 0.08)
                Genes that do not have a single accession having an expression
                value above this threshold are removed from analysis. These are
                likely presence/absence and will not have a strong coexpression
                pattern.
            force : bool
                If true and a Coex already exists with the same name, it will
                be overwritten.
        '''
        if force:
            m80.Tools.delete('Coex',name)
        # Sanity check
        if m80.Tools.available('Coex',name):
            raise CoexExistsError(f'{name} already exists.') 
        # Create an empty class
        self = cls(name)
        # get the list of loci that are in the df
        log.info("Creating filtered Loci object")
        # we need to do this using the minus80 API so its fast
        lids = []
        cur = loci.m80.db.cursor()
        
        for (LID,) in cur.executemany(
            'SELECT LID from loci WHERE name = ?',((name,) for name in df.index)        
        ): 
            

        # Calculate QC for loci and accessions
        qc_loci = pd.DataFrame({"has_id": True}, index=df.index)
        qc_accession = pd.DataFrame({"has_id": True}, index=df.columns)
        qc_loci['pass_in_loci'] = [x in loci for x in df.index]
        log.info(
            f"Found {sum(qc_loci['pass_in_loci'] == False)} "
            f"df row IDs not in {loci.m80.name}"
        )

        # Set the raw data to be the original df
        self.m80.col['raw_expr'] = df
        if normalize_expr_values:
            # Apply inverse hyperbolic sine 
            # per https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3634062/
            df = df.apply(np.arcsinh, axis=0)
        # Remove DF entries that are not in Loci

        # Return an instance
        return self

    #-----------------------------------------
    #       Static Methods
    #-----------------------------------------
    
