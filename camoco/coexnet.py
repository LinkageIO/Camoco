import warnings

import pandas as pd
import minus80 as m80
import locuspocus as lp

from typing import List

from .exceptions import (
    CoexNetEmptyError,        
)


class CoexNet(m80.Freezable):

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
    ):
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
        loci: lp.Loci
        normalize_expr_values: bool = True
    ) -> "CoexNet":
        '''
            Create a CoexNet from a Pandas DataFrame and 
            a LocusPocus Loci object.
        '''
        self = cls(name)
        # Set the raw data to be the original df
        self.m80.col['raw_expr'] = df
        if normalize_expr_values:
            # Apply inverse hyperbolic sine 
            # per https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3634062/
            df = df.apply(np.arcsinh, axis=0)

        # Return an instance
        return self

    #-----------------------------------------
    #       Static Methods
    #-----------------------------------------
    
