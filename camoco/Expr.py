#! /usr/bin/python3
from .Camoco import Camoco
from .RefGen import RefGen
from .Tools import memoize
from .Locus import Locus
from .Exceptions import CamocoGeneNameError,CamocoAccessionNameError,CamocoGeneAbsentError

from scipy.spatial.distance import pdist, squareform, euclidean
from scipy.stats import hypergeom, pearsonr
from scipy.stats.mstats import rankdata as mrankdata
from scipy.cluster.hierarchy import linkage, dendrogram
from collections import defaultdict,Counter

import matplotlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import io
import re
import string

pd.set_option('display.width', 100)

class Expr(Camoco):
    '''
        A gene expression dataset. Build, normalize, filter and 
        easily access different parts of the gene expression matrix.
    '''
    def __init__(self, name):
        # Create a camoco object
        super().__init__(name=name, type='Expr')
        # Part I: Load the Expression dataset
        self.log('Loading Expr table')
        self._expr = self._bcolz('expr')
        self._gene_qc_status = self._bcolz('gene_qc_status')
        if (self._expr is None) or (self._gene_qc_status is None):
            self._expr = pd.DataFrame()
        
        self.log('Building Expr Index')
        self._expr_index = defaultdict(
            lambda: None,
            {gene:index for index, gene in enumerate(self._expr.index)}
        )
        # Part II: Load the Reference Genome
        try:
            self.log('Loading RefGen')
            self.refgen = RefGen(self.refgen)
        except TypeError as e:
            self.log('RefGen for {} not set!', self.name)
        except NameError as e:
            self.log.warn('Refgen for {} not available, must be reset!', self.name)

    def __contains__(self, obj):
        if obj in self._expr.index:
            return True
        if obj in self._expr.columns:
            return True
        try:
            if obj.id in self._expr.index:
                return True
        except AttributeError as e:
            pass
        return False

    def __repr__(self):
        return ""

    def __str__(self):
        pass

    def num_genes(self,raw=False):
        return len(self.expr(raw=raw))

    def num_accessions(self,raw=False):
        return len(self.expr(raw=raw).columns)

    def shape(self):
        return self._expr.shape

    def zscore(self):
        pass

    def accessions(self):
        return self._expr.columns

    def genes(self, raw=False):
        # Returns a list of distinct genes
        if raw is False:
            return self.refgen.from_ids(self._expr.index)
        else:
            return self.refgen.from_ids(self._bcolz('raw_expr').index)

    def expr_profile(self, gene):
        '''
            return the expression profile for a gene
        '''
        return self._expr.loc[gene.id]

    def is_normalized(self, max_val=None, raw=False):
        if max_val is not None:
            max_val = max_val # Use the user defined max val
        elif self.rawtype.upper() == 'RNASEQ':
            max_val = 1100
        elif self.rawtype.upper() == 'MICROARRAY':
            max_val = 100
        else:
            max_val = 0
        return self._expr.apply(
            lambda col: np.nanmax(col.values) < max_val, axis=0
        )

    def max_values(self, axis=0):
        return np.nanmax(self._expr, axis=axis)

    def anynancol(self):
        '''
            A gut check method to make sure none of the expression columns
            got turned into all nans. Because apparently that is a problem.
        '''
        return any(self._expr.apply(lambda col: all(np.isnan(col)), axis=0))

    def expr(self, genes=None, accessions=None,
             raw=False,gene_normalize=False):
        '''
            Access raw and QC'd expression data.

            Parameters
            ----------
            genes : iterable of camoco.Locus objects (default: None)
                If not None, this will retrieve the expression values for
                the loci specified within the iterable, otherwise it will
                include ALL loci in the expr dataset
            accessions : iterable of str (default: None)
                If not None, will retrieve expression values for the
                accessions (experiments) specified, otherwise will
                retrieve ALL accessions.
            raw : bool (default: False)
                Flag to indicate on using the raw table versus the current
                expr table. See the transformation_log for more details on
                the difference.
            gene_normalize : bool (default: False)
                Perform standard normalization on gene-wise data
            zscore : bool (default: False)

        '''
        if raw is True:
            self.log('Extracting raw expression values')
            df = self._bcolz('raw_expr')
        else:
            df = self._expr
        if genes is not None:
            df = df.loc[[x.id for x in genes], :]
        if accessions is not None:
            df = df[accessions]
        if gene_normalize:
            df = df.apply(
                # Axis: 1 applies to ROWS!
                lambda row: (row-row.mean())/row.std(), axis=1
            )
        return df

    def plot_accession_histograms(self, bins=50, figsize=(16, 8)):
        '''
            Plot histogram of accession expression values.
        '''
        raw = self._bcolz('raw_expr')
        qcd = self._expr

        for name, values in qcd.iteritems():
            raw_values = raw[name]
            # Shorten name
            if len(name) > 20:
                name = name[0:20] + '...' + name[-11:-1]
            self.log('Plotting values for {}', name)
            # Extract out the raw values
            raw_valid = np.ma.masked_invalid(raw_values)
            # Extract out the normalized values
            valid = np.ma.masked_invalid(values)
            # Plot histograms
            f = plt.figure(figsize=figsize)
            plt.subplot(121)
            plt.hist(raw_valid[~raw_valid.mask], bins=bins)
            plt.xlim(-15, 15)
            plt.title('{}:{}'.format(self.name, name))
            plt.ylabel('Frequency')
            plt.subplot(122)
            plt.hist(valid[~valid.mask], bins=bins)
            plt.xlabel('Expression')
            plt.xlim(-15, 15)

            plt.savefig("ACC_HIST_{}:{}.png".format(self.name, name))
            plt.close(f)


    '''
        Internal Methods ------------------------------------------------------
    '''

    def _update_values(self, df, transform_name, raw=False):
        '''
            updates the 'expression' table values with values from df.
            Requires a transformation name for the log.
            Option to overwrite raw table or working table.

            Parameters
            ----------
            df : DataFrame
                Updates the internal values for the Expr object
                with values in the data frame.
            transform_name : str
                A short justification for what was done to the
                updated values.
            raw : bool (default: False)
                A flag to update the raw values. This also resets
                the current values to what is in df.

            Returns
            -------
            self : Expr Object

            Raises:
            ------
            CamocoGeneNamesError
            CamocoAccessNamesError
        '''
        # update the transformation log
        if len(set(df.columns)) != len(df.columns):
            raise CamocoAccessionNameError('Accession names must be unique')
        if len(set(df.index)) != len(df.index):
            raise CamocoGeneNameError('Gene names must be unique.')
        self._transformation_log(transform_name)
        if raw == True:
            table = 'raw_expr'
            # If we are updating the raw table, remove the
            # normal table since it assumes it came from
            # the raw table.
            self._reset(raw=False)
        else:
            table = 'expr'
            # Keep full names in raw, but compress the
            # names in the normed network
            def shorten(x):
                if len(x) > 100:
                    return x[0:89] + '...' + x[-10:-1]
                else:
                    return x
            df.columns = [shorten(x) for x in df.columns]
        # Sort the table by genes
        df = df.sort_index()
        # ensure that column names are alphanumeric
        colP = re.compile('[^A-Za-z0-9_]')
        begP = re.compile('^\d')
        df.columns = [colP.sub('_', x).strip('_') for x in df.columns.values]
        df.columns = [x if not begP.match(x[0]) else 'Exp_'+x for x in df.columns.values]
        # Also, make sure gene names are uppercase
        idxP = re.compile('[^A-Za-z0-9_, ;:().]')
        df.index = [idxP.sub('', str(x)).upper() for x in df.index.values]
        try:
            self._bcolz(table, df=df)
            self._expr = df
        except Exception as e:
            self.log('Unable to update expression table values: {}', e)
            raise e
        # Set the index
        self._expr_index = defaultdict(
            lambda: None,
            {gene:index for index, gene in enumerate(self._expr.index)}
        )
        return self

    def _get_gene_index(self,gene):
        '''
            Retrieve the row index for a gene.
            
            Parameters
            ----------
            gene : co.Locus object
                The gene object the get the index for

            Returns
            -------
            an integer containing the expr dataframe index

            Raises
            ------
            CamocoGeneAbsentError
                If the gene requested is not in the Expr dataframe
        '''
        if isinstance(gene,Locus):
            id = gene.id
        else:
            id = gene
        index = self._expr_index[id]
        if index == None:
            raise CamocoGeneAbsentError('{} not in {}'.format(id,self.name))
        return index

    def _transformation_log(self, transform=None):
        if transform is None:
            return self._global('transformation_log')
        elif transform == 'reset' or self._global('transformation_log') is None:
            self._global('transformation_log', 'raw')
        else:
            self._global(
                'transformation_log',
                self._global('transformation_log') + '->' + str(transform)
            )
            self.log('Trans. Log: {}', self._global('transformation_log'))

    def _reset(self, raw=False):
        '''
            resets the expression values to their raw
            state undoing any normalizations
        '''
        if raw:
            # kill the raw table too
            self.log('Resetting raw expression data')
            self._bcolz('raw_expr', df=pd.DataFrame())
        self.log('Resetting expression data')
        self._expr = self.expr(raw=True)
        self._bcolz('expr', df=self._expr)
        self._transformation_log('reset')


    def _normalize(self, norm_method=None, max_val=None, **kwargs):
        ''' 
            Evaluates QC expression data and re-enters 
            normalized data into database 

            Parameters
            ----------
            norm_method : The normalization method to use. This can be inferred
                from the raw data type. By default RNASeq uses np.arcsinh and
                microarray data uses np.log2. A different normalization function
                can be passed directly in. 
                Default: None (inferred from Expr.rawtype)
            max_val : This value is used to determine if any columns of the 
                dataset have already been normalized. If any 'normailzed' 
                values in an Accession column is larger than max_val, an
                exception is thown. max_val is determined by Expr.raw_type
                (default 100 for MicroArray and 1100 for RNASeq) but a 
                max_val can be passed in to override these defaults.

        '''
        self.log('------------ Normalizing')
        if all(self.is_normalized(max_val=max_val)):
            self.log("Dataset already normalized")
            self._transformation_log('DetectedPreNormalized')
        elif any(self.is_normalized(max_val=max_val)):
            raise TypeError(
                ('Attempting normalization on already normalized'
                ' dataset. Consider passing a --max_val '
                '< {} if Im wrong.').format( min(self.max_values())))
        else:
            df = self._expr
            if norm_method is not None:
                method = norm_method
            elif self.rawtype.upper() == 'RNASEQ':
                method = np.arcsinh
            elif self.rawtype.upper() == 'MICROARRAY':
                method = np.log2
            else:
                raise ValueError(
                    ('Could not guess correct normalization for {}'
                    ' pass in function through method argument.'
                    ).format(self.rawtype))
            # apply the normalization to each column (accession)
            df = df.apply(lambda col: method(col), axis=0)
            # update values
            self._update_values(df, method.__name__)

    def _quality_control(self, min_expr=0.01, max_gene_missing_data=0.2, \
        min_single_sample_expr=5, max_accession_missing_data=0.3, \
        membership=None, dry_run=False, presence_absence=False, **kwargs):
        '''
            Perform Quality Control on raw expression data. This method filters
            genes based on membership to some RefGen instance, filters based on
            a minimum FPKM or equivalent expression value, filters out genes
            and accessions with too much missing data, filters out genes which
            are lowly expressed (do not have at least one accession that meets
            an FPKM threshold, i.e. likely presence absense). See parameters
            for more details.

            Parameters
            ----------
            min_expr : int (default: 0.01)
                FPKM (or equivalent) values under this threshold will be set to
                NaN and not used during correlation calculations.
            max_gene_missing_data : float (default: 0.2)
                Maximum percentage missing data a gene can have. Genes under
                this are removed from dataset.
            min_single_sample_expr : int (default: 5)
                Genes that do not have a single accession having an expression
                value above this threshold are removed from analysis. These are
                likely presence/absence and will not have a strong coexpression
                pattern.
            max_accession_missing_data : float (default: 0.5)
                maximum percentage missing data an accession (experiment) can
                have before it is removed.
            membership : RefGen
                Genes which are not contained within this RefGen will be
                removed. Note: this could also be another object that will
                implement an interface that will check to see if gene ids are
                contained within it i.e. a set of gene ids.
            dry_run : bool (default: False)
                Used in testing to speed up calculations. Limits the QC
                dataframe to only have 100 genes.
            presence_absence : bool (default: False)
                Used to convert 0's within the data to a 0.001 after min 
                expression values are filtered out to allow for presence
                absence variation
        '''
        self.log('------------Quality Control')
        df = self.expr()
        # remember how we set the flags
        self._global('qc_min_expr', min_expr)
        self._global('qc_max_gene_missing_data', max_gene_missing_data)
        self._global('qc_min_single_sample_expr', min_single_sample_expr)
        self._global('qc_max_accession_missing_data', max_accession_missing_data)
        # Retrieve raw data as a data frame
        self.log('Raw Starting set: {} genes {} accessions'.format(
            len(df.index), len(df.columns))
        )
        # Remember why we remove certain genes
        # If TRUE it passes, if FALSE it fails!!!
        qc_gene = pd.DataFrame({'has_id':True}, index=df.index)
        qc_accession = pd.DataFrame({'has_id':True}, index=df.columns)
        # -----------------------------------------
        # Gene Membership test
        if not membership:
            membership = self.refgen
        self._global('qc_membership', str(membership))
        qc_gene['pass_membership'] = [x in membership for x in df.index]
        self.log(
            "Found out {} genes not in {}", 
            sum(qc_gene['pass_membership'] == False), 
            membership
        )

        # -----------------------------------------
        # Set minimum FPKM threshold
        self.log("Filtering expression values lower than {}", min_expr)
        df_flt = df.copy()
        # Presence absence variable et
        if presence_absence == True:
           self.log("Allowing for presence absence variation")
           #find out which values equal 0
           zero_index = df_flt == 0
        # Filter the min expression genes
        df_flt[df < min_expr] = np.nan
        if presence_absence == True:
            #change out original 0's index to a small value
            df_flt[zero_index] = 0.001
        df = df_flt
        # -----------------------------------------
        # Gene Missing Data Test
        qc_gene['pass_missing_data'] = df.apply(
            lambda x : ((sum(np.isnan(x))) < len(x)*max_gene_missing_data),
            axis=1
        )
        self.log(
            "Found {} genes with > {} missing data",
            sum(qc_gene['pass_missing_data'] == False),
            max_gene_missing_data
        )
        # -----------------------------------------
        # Gene Min Expression Test
        # filter out genes which do not meet a minimum expr
        # threshold in at least one sample
        qc_gene['pass_min_expression'] = df.apply(
            lambda x: any(x >= min_single_sample_expr),
            axis=1 # 1 is column
        )
        self.log(
            ("Found {} genes which "
            "do not have one sample above {}"), 
            sum(qc_gene['pass_min_expression'] == False),
            min_single_sample_expr
        )
        qc_gene['PASS_ALL'] = qc_gene.apply(
            lambda row: np.all(row), axis=1
        )
        df = df.loc[qc_gene['PASS_ALL'], :]
        # -----------------------------------------
        # Filter out ACCESSIONS with too much missing data
        qc_accession['pass_missing_data'] = df.apply(
            lambda col : (
                ((sum(np.isnan(col)) / len(col)) <= max_accession_missing_data)
            ),
            axis=0 # 0 is columns
        )
        self.log(
            "Found {} accessions with > {} missing data", 
            sum(qc_accession['pass_missing_data'] == False),
            max_accession_missing_data
        )
        # Update the total QC passing column
        qc_accession['PASS_ALL'] = qc_accession.apply(
            lambda row: np.all(row), axis=1
        )
        df = df.loc[:, qc_accession['PASS_ALL']]
        # Update the database
        self._bcolz('qc_accession', df=qc_accession)
        self._bcolz('qc_gene', df=qc_gene)
        # Report your findings
        self.log('Genes passing QC:\n{}', str(qc_gene.apply(sum, axis=0)))
        self.log('Accessions passing QC:\n{}', str(qc_accession.apply(sum, axis=0)))
        # Also report a breakdown by chromosome
        qc_gene = qc_gene[qc_gene['pass_membership']]
        qc_gene['chrom'] = [self.refgen[x].chrom for x in qc_gene.index]
        self.log('Genes passing QC by chromosome:\n{}',
            str(qc_gene.groupby('chrom').aggregate(sum, axis=0))
        )
        # update the df to reflect only genes/accession passing QC
        self.log('Kept: {} genes {} accessions'.format(len(df.index), len(df.columns)))
        if dry_run:
            # If dry run, take first 100 rows of QC
            self.log.warn("Dry Run")
            df = df.iloc[0:100,:]
        self._update_values(df, 'quality_control')

    @staticmethod
    def inplace_nansort(col):
        # mask invalid data
        masked_col = np.ma.masked_invalid(col)
        masked_sorted = np.sort(col[~masked_col.mask].data)
        # get ranked values
        col_sorted = np.copy(col)
        non_nan = 0
        for i, x in enumerate(~masked_col.mask):
            if x == True:
                col_sorted[i] = masked_sorted[non_nan]
                non_nan += 1
            else:
                col_sorted[i] = np.nan
        return col_sorted

    def _quantile(self):
        '''
            Perform quantile normalization across each accession.
            Each accessions gene expression values are replaced with
            ranked gene averages.
        '''
        self.log('------------ Quantile ')
        if 'quantile' in self._transformation_log():
            raise ValueError('Quantile already performed on {}', self.name)
        # Retrieve current expression DataFrame
        expr = self.expr()
        self.log('Ranking data')
        for accession_name,values in expr.iteritems():
            rank_ties = max(Counter(values).values())
            if rank_ties > len(values) * 0.20:
                raise ValueError(
                    f'{self.name}:{accession_name} has {rank_ties} '
                    f'({rank_times/len(values)}%) rank ties'
                )
        # assign ranks by accession (column)
        expr_ranks = expr.rank(axis=0, method='first', na_option='keep')
        assert np.all(np.isnan(expr) == np.isnan(expr_ranks))
        # normalize rank to be percentage
        expr_ranks = expr_ranks.apply(
            lambda col: col/np.nanmax(col.values),
            axis=0
        )
        # we need to know the number of non-nans so we can correct for their ranks later
        self.log('Sorting ranked data')
        # Sort values by accession/column, lowest to highest
        expr_sort = expr.apply(lambda col: self.inplace_nansort(col), axis=0)
        # make sure the nans weren't included in the sort or the rank
        assert np.all(np.isnan(expr) == np.isnan(expr_ranks))
        assert np.all(np.isnan(expr) == np.isnan(expr_sort))
        # calculate ranked averages
        self.log('Calculating averages')
        rank_average = expr_sort.apply(np.nanmean, axis=1)
        # we need to apply the percentages to the lenght of the
        rankmax = len(rank_average)
        self.log('Range of normalized values:{}..{} (n = {})'.format(
            min(rank_average), max(rank_average), len(rank_average))
        )
        self.log('Asserting that no Genes are nan...')
        assert sum(np.isnan(rank_average)) == 0
        self.log('Applying non-floating normalization')
        quan_expr = expr_ranks.applymap(
            lambda x : rank_average[int(x*rankmax)-1] if not np.isnan(x) else np.nan
        )
        self.log('Updating values')
        assert np.all(np.isnan(expr) == np.isnan(quan_expr))
        self._update_values(quan_expr, 'quantile')


    @property
    def _parent_refgen(self):
        return RefGen(self._global('parent_refgen'))

    def _set_refgen(self, refgen, filter=True):
        '''
            Sets the current refgen. Its complicated.
        '''
        # Keep a record of parent refgen
        self._global('parent_refgen', refgen.name)
        # Filter down to only genes in
        if filter:
            refgen = refgen.filtered_refgen(
                'Filtered{}'.format(self.name),
                'Filtered Refgen',
                refgen,
                self.genes(),
            )
        # remember to set for current instance
        self._global('refgen', refgen.name)
        self.refgen = refgen

    @property
    def _cmap(self):
        '''
            Used for the heatmap function. Retruns a matplotlib cmap which is yellow/blue
        '''
        white_middle_heatmapdict = {
            'red': ((0.0, 1.0, 1.0),
                    (0.5, 1.0, 1.0),
                    (1.0, 0.0, 0.0)),
            'green':((0.0, 1.0, 1.0),
                    (0.5, 1.0, 1.0),
                    (1.0, 0.0, 0.0)),
            'blue': ((0.0, 0.0, 0.0),
                    (0.5, 1.0, 1.0),
                    (1.0, 1.0, 1.0))
        }
        heatmapdict = {
            'red': ((0.0, 1.0, 1.0),
                    (0.5, 0.0, 0.0),
                    (1.0, 0.0, 0.0)),
            'green':((0.0, 1.0, 1.0),
                    (0.5, 0.0, 0.0),
                    (1.0, 0.0, 0.0)),
            'blue': ((0.0, 0.0, 0.0),
                    (0.5, 0.0, 0.0),
                    (1.0, 1.0, 1.0))}
        heatmap_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap', heatmapdict, 256)
        return heatmap_cmap


    ''' ------------------------------------------------------------------------------------------
            Class Methods
    '''

    @classmethod
    def create(cls, name, description, refgen, type='Expr'):
        '''
            Create an empty Expr instance. Overloads the Camoco
            create method. See Camoco.create(...)

            Parameters
            ----------
            name : str
                A name for the Expr object to reference in the Camoco database
            description : str
                A short description for the dataset
            refgen : camoco.RefGen
                A Camoco refgen object which describes the reference
                genome referred to by the genes in the dataset. This
                is cross references during import so we can pull information
                about genes we are interested in during analysis.

            Returns
            -------
                An empty Expr instance

        '''
        # Piggy back on the super create method
        self = super().create(name, description,type=type)
        # Create appropriate bcolz tables
        self._bcolz('expr', df=pd.DataFrame())
        self._bcolz('raw_expr', df=pd.DataFrame())
        # Delete existing datasets
        self._set_refgen(refgen, filter=False)
        return self

    @classmethod
    def from_table(cls, filename, name, description, refgen, rawtype=None,
            sep='\t', normalize=True, quality_control=True, **kwargs):
        '''
            Create a Expr instance from a file containing raw expression data.
            For instance FPKM or results from a microarray experiment. This is
            a convenience method which reads the table in to a pandas DataFrame
            object and passes the object the Expr.from_DataFrame(...). See the
            doc on Expr.from_DataFrame(...) for more options.

            Parameters
            ----------
            filename : str (path)
                a path the the table containing the raw expression data.
            name : str
                A short name to refer to from the camoco dataset API.
            description : str
                A short description for the dataset
            refgen : camoco.RefGen
                A Camoco refgen object which describes the reference
                genome referred to by the genes in the dataset. This
                is cross references during import so we can pull information
                about genes we are interested in during analysis.
            rawtype : str (default: None)
                This is noted here to reinforce the impotance of the rawtype
                passed to camoco.Expr.from_DataFrame. See docs there for more
                information.
            sep : str (default: \t)
                Column delimiter for the data in filename path
            normalize : bool (Default: True)
                Specifies whether or not to normalize the data so raw
                expression values lie within a log space. This is best
                practices for generating interpretable expression analyses. See
                Expr._normalize method for more information.  info.
            quality_control : bool (Default: True)
                A flag which specifies whether or not to perform QC. Parameters
                for QC are passed in using the **kwargs arguments. For default
                parameters and options see Expr._quality_control.
            **kwargs : key value pairs
                additional parameters passed to subsequent methods. (see
                Expr.from_DataFrame)

            Returns
            -------
            An Expr instance

        '''
        tbl = pd.read_table(filename, sep=sep)
        return cls.from_DataFrame(
                tbl, name, description, refgen, 
                rawtype=rawtype, **kwargs
            )

    @classmethod
    def from_DataFrame(cls, df, name, description, refgen, rawtype=None,
            normalize=True, norm_method=None, quantile=False, 
            quality_control=True, **kwargs
        ):
        ''' 
            Creates an Expr instance from a pandas DataFrame. Expects that the
            DataFrame index is gene names and the column names are accessions
            (i.e. experiments).  This is the preferred method for creating an
            Expr instance, in other words, other classmethods transform their
            data so they can call this method.

            Parameters
            ----------
            df : pandas.DataFrame
                a DataFrame containing expression data. Assumes index is the
                genes and columns is the accessions (experiment names)
            name : str
                A short name to refer to from the camoco dataset API.
            description : str
                A short description for the dataset
            refgen : camoco.RefGen
                A Camoco refgen object which describes the reference
                genome referred to by the genes in the dataset. This
                is cross references during import so we can pull information
                about genes we are interested in during analysis.
            rawtype : str (one of: 'RNASEQ' or 'MICROARRAY')
                Specifies the fundamental datatype used to measure expression.
                During importation of the raw expression data, this value is
                used to make decisions in converting data to log-space.
            normalize : bool (Default: True)
                Specifies whether or not to normalize the data so raw
                expression values lie within a log space. This is best
                practices for generating interpretable expression analyses. See
                Expr._normalize method for more information.  info.
            norm_method : None OR python function
                If rawtype is NOT RNASEQ or MICROARRY AND normalize is still
                True, the normalization method for the raw expression values
                needs to be passed in. This is for extreme customization
                situations.
            quantile : bool (Default : False)
                Specifies whether or not to perform quantile normalization on
                import.
            quality_control : bool (Default: True)
                A flag which specifies whether or not to perform QC. Parameters
                for QC are passed in using the **kwargs arguments. For default
                parameters and options see Expr._quality_control.
            **kwargs : key value pairs
                additional parameters passed to subsequent methods.
                See arguments for Expr._normalize(), Expr._quality_control()

            Returns
            -------
            An Expr instance

        '''
        # we are all pandas on the inside O.O
        self = cls.create(name, description, refgen)
        self._reset(raw=True)
        if rawtype is None:
            raise TypeError("raw_type must be one of ['RNASEQ', 'MICROARRAY']")
        self._global('rawtype', rawtype)
        # put raw values into the database
        self.log('Importing Raw Expression Values')
        self._update_values(df, 'Raw'+rawtype, raw=True)
        if quality_control:
            self.log('Performing Quality Control on genes')
            self._quality_control(**kwargs)
            assert self.anynancol() == False
        else:
            self.log('Skipping Quality Control!')
        if normalize:
            self.log('Performing Raw Expression Normalization')
            self._normalize(**kwargs)
            assert self.anynancol() == False
        if quantile:
            self.log('Performing Quantile Gene Normalization')
            self._quantile()
            assert self.anynancol() == False
        self.log('Filtering refgen: {}', refgen.name)
        self._set_refgen(refgen, filter=True)
        return self
