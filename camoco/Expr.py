#! /usr/bin/env python3

from camoco.Camoco import Camoco
from camoco.RefGen import RefGen
from camoco.Tools import memoize
from scipy.spatial.distance import pdist, squareform, euclidean
from scipy.stats import hypergeom,pearsonr
from scipy.cluster.hierarchy import linkage, dendrogram
from collections import defaultdict

import matplotlib
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import io               
import re
import string

class Expr(Camoco):
    ''' 
        A representation of gene expression. The base data structure for this 
        guy is a 
    '''
    def __init__(self,name):
        super().__init__(name=name,type='Expr') 
        # Part I: Load the Expression dataset
        try:
            self.log('Loading Expr table')
            # Open the HDF5 store
            self.hdf5 = self._hdf5(name)
            # Load the expression Data
            self._expr = self.hdf5['expr']
            self._gene_qc_status = self.hdf5['gene_qc_status']
        except KeyError as e:
            self.log('{} is empty: ({})',name,e)
            self._expr = pd.DataFrame()
        self.log('Building Expr Index')
        self._expr_index = defaultdict(
            lambda: None,
            {gene:index for index,gene in enumerate(self._expr.index)}
        ) 
        # Part II: Load the Reference Genome
        try:
            self.log('Loading RefGen')
            self.refgen = RefGen(self.refgen)
        except TypeError as e:
            self.log.warn('RefGen for {} not set!',self.name)
        except NameError as e:
            self.log.warn('Refgen for {} not available, must be reset!',self.name)

    def __contains__(self,obj):
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

    def num_genes(self):
        return len(self._expr.index)

    def num_accessions(self):
        return len(self._expr.columns)

    def shape(self):
        return self._expr.shape

    def zscore(self):
        pass

    def accessions(self):
        return self._expr.columns

    def genes(self,raw=False):
        # Returns a list of distinct genes 
        if raw is False:
            return self.refgen.from_ids(self._expr.index)
        else:
            return self.refgen.from_ids(self.hdf5['raw_expr'].index)

    def expr_profile(self,gene):
        '''
            return the expression profile for a gene
        '''
        return self._expr.loc[gene.id]

    def is_normalized(self,max_val=None,raw=False):
        if max_val is not None:
            max_val = max_val # Use the user defined max val
        elif self.rawtype.upper() == 'RNASEQ':
            max_val = 1100 
        elif self.rawtype.upper() == 'MICROARRAY':
            max_val = 100 
        return self._expr.apply(
            lambda col: np.nanmax(col.values) < max_val, axis=0
        )
    
    def anynancol(self):
        ''' 
            A gut check method to make sure none of the expression columns
            got turned into all nans. Because apparently that is a problem.
        '''
        return any(self._expr.apply(lambda col: all(np.isnan(col)),axis=0))

    def expr(self,genes=None,accessions=None,raw=False):
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
            zscore : bool (default: False)
                
        '''
        df = self.hdf5['raw_expr'] if raw is True else self._expr
        if genes is not None:
            df = df.loc[[x.id for x in genes],:]
        if accessions is not None:
            df = df[accessions]
        return df

    '''
        Internal Methods ------------------------------------------------------
    '''

    def _update_values(self,df,transform_name,raw=False):
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

        '''
        # update the transformation log
        self._transformation_log(transform_name)
        table = 'raw_expr' if raw else 'expr'
        # Sort the table by genes
        df = df.sort()
        # ensure that column names are alphanumeric
        # hdf5 doesn't like unicode characters
        pattern = re.compile('[^\w ,;]+')
        df.columns = [pattern.sub('',x) for x in df.columns.values ]
        # Also, make sure gene names are uppercase
        df.index = [pattern.sub('',x).upper() for x in df.index.values ]
        try:
            self._expr = df
            self.hdf5[table] = df
            self.hdf5.flush(fsync=True)
        except Exception as e:
            self.log('Unable to update expression table values: {}',e)
            raise
        if raw == True:
            # Get rid of the current values
            self._reset(raw=False)

    def _transformation_log(self,transform=None):
        if transform is None:
            return self._global('transformation_log')
        elif transform == 'reset' or self._global('transformation_log') is None:
            self._global('transformation_log','raw')
        else:
            self._global(
                'transformation_log',
                self._global('transformation_log') + '->' + str(transform)
            )
            self.log('Trans. Log: {}',self._global('transformation_log'))
 
    def _reset(self,raw=False):
        ''' resets the expression values to their raw state undoing any normalizations '''
        if raw:
            # kill the raw table too
            self.log('Resetting raw expression data')
            self.hdf5['raw_expr'] = pd.DataFrame()
        self.log('Resetting expression data')
        self.hdf5['expr'] = self._expr = pd.DataFrame()
        self._transformation_log('reset')


    def _normalize(self,norm_method=None,is_raw=None,max_val=None,**kwargs):
        ''' evaluates qc expression data and re-enters 
            normaized data into database '''
        self.log('------------ Normalizing')
        if all(self.is_normalized(max_val=max_val)):
            self.log("Dataset already normalized")
            self._transformation_log('DetectedPreNormalized')
        elif any(self.is_normalized(max_val=max_val)):
            # Something fucked up is happending
            raise TypeError(
                ('Attempting normalization on already normalized'
                ' dataset. Consider passing a max_val '
                '< {} if Im wrong.').format(min(self.max_values())))
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
            df = df.apply(lambda col: method(col),axis=0)
            # update values
            self._update_values(df,method.__name__)

    def _quality_control(self,min_expr=1,max_gene_missing_data=0.2,\
        min_single_sample_expr=5, max_accession_missing_data=0.5,\
        membership=None, dry_run=False, **kwargs):
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
            min_expr : int (default: 1)
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
                dataframe to only have 5000 genes.
        '''        
        self.log('------------Quality Control')
        df = self.expr(raw=True)
        # remember how we set the flags
        self._global('qc_min_expr',min_expr)
        self._global('qc_max_gene_missing_data',max_gene_missing_data)
        self._global('qc_min_single_sample_expr',min_single_sample_expr)
        self._global('qc_max_accession_missing_data',max_accession_missing_data)
        # retrieve raw data as a data frame
        self.log('Raw Starting set: {} genes {} accessions'.format(
            len(df.index),len(df.columns))
        )
        # Remember why we remove certain genes
        # If TRUE it passes, if FALSE it fails!!!
        qc_gene = pd.DataFrame({'has_id':True},index=df.index)
        qc_accession = pd.DataFrame({'has_id':True},index=df.columns)

        # -----------------------------------------
        # Gene Membership test
        if not membership:
            membership = self.refgen
        self._global('qc_membership',str(membership))
        self.log("Filtering out genes not in {}",membership)
        qc_gene['pass_membership'] = [x in membership for x in df.index]

        # -----------------------------------------
        # Set minimum FPKM threshold
        self.log("Filtering out expression values lower than {}",min_expr)
        df_flt = df.copy()
        df_flt[df < min_expr] = np.nan
        df = df_flt
        # -----------------------------------------
        # Gene Missing Data Test
        self.log(
            "Filtering out genes with > {} missing data",
            max_gene_missing_data
        )
        qc_gene['pass_missing_data'] = df.apply(
            lambda x : ((sum(np.isnan(x))) < len(x)*max_gene_missing_data),
            axis=1
        )
        # -----------------------------------------
        # Gene Min Expression Test
        # filter out genes which do not meet a minimum expr 
        # threshold in at least one sample
        self.log(("Filtering out genes which "
                  "do not have one sample above {}"),min_single_sample_expr)
        qc_gene['pass_min_expression'] = df.apply(
            lambda x: any(x >= min_single_sample_expr),
            axis=1 # 1 is column
        )
        qc_gene['PASS_ALL'] = qc_gene.apply(
            lambda row: np.all(row),axis=1
        )
        df = df.loc[qc_gene['PASS_ALL'],:] 
        # -----------------------------------------
        # Filter out ACCESSIONS with too much missing data
        self.log("Filtering out accessions with > {} missing data",max_accession_missing_data)
        qc_accession['pass_missing_data'] = df.apply(
            lambda col : (
                ((sum(np.isnan(col)) / len(col)) <= max_accession_missing_data)
            ),
            axis=0 # 0 is columns
        )
        # Update the total QC passing column
        qc_accession['PASS_ALL'] = qc_accession.apply(
            lambda row: np.all(row),axis=1
        )
        df = df.loc[:,qc_accession['PASS_ALL']] 
        # Update the database
        self.hdf5['qc_accession'] = qc_accession
        self.hdf5['qc_gene'] = qc_gene
        # Report your findings
        self.log('Gene Removal:\n{}',str(qc_gene.apply(sum,axis=0)))
        self.log('Accession Removal:\n{}',str(qc_accession.apply(sum,axis=0)))
        # Also report a breakdown by chromosome
        qc_gene = qc_gene[qc_gene['pass_membership']]
        qc_gene['chrom'] = [self.refgen[x].chrom for x in qc_gene.index] 
        self.log('Genes passing QC by chromosome:\n{}',
            str(qc_gene.groupby('chrom').aggregate(sum,axis=0))
        )
        # update the df to reflect only genes/accession passing QC
        self.log('Kept: {} genes {} accessions'.format(len(df.index),len(df.columns)))
        if dry_run:
            # If dry run, take first 100 rows of QC
            self.log.warn("Dry Run")
            df = df.iloc[0:5000,:]
        self._update_values(df,'quality_control')

    def _quantile(self):
        ''' 
            Perform quantile normalization across each accession. 
            Each accessions gene expression values are replaced with 
            ranked gene averages.
        '''
        raise NotImplementedError('This method is BROKEN!')
        # get gene by accession matrix
        self.log('------------ Quantile ')
        expr = self._expr
        self.log('Ranking data')
        # assign ranks by accession (column)
        expr_ranks = expr.rank(axis=0,method='dense',na_option='keep')
        # normalize rank to be percentage
        expr_ranks = expr_ranks.apply(lambda col: col/np.nanmax(col.values), axis=0)
        # we need to know the number of non-nans so we can correct for their ranks later
        self.log('Sorting ranked data')
        # assign accession values by order
        # NOTE this currently keeps nans where they are. It COULD change with newer versions of pandas. 
        expr_sort = expr.sort(axis=0)
        # make sure the nans weren't included in the sort or the rank
        assert np.all(np.isnan(expr) == np.isnan(expr_ranks))
        assert np.all(np.isnan(expr) == np.isnan(expr_sort))
        # calculate ranked averages
        self.log('Calculating averages')
        rank_average = expr_sort.apply(lambda row: np.mean(row[np.logical_not(np.isnan(row))]),axis=1)
        rankmax = len(rank_average)
        self.log('Range of normalized values:{}..{}'.format(min(rank_average),max(rank_average)))
        self.log('Applying non-floating normalization')
        expr = expr_ranks.applymap(
            lambda x : rank_average[int(x*rankmax)-1] if not np.isnan(x) else np.nan
        )
        self.log('Updating values')
        self._update_values(expr,'quantile')

    @property
    def _parent_refgen(self):
        return RefGen(self._global['parent_refgen'])

    def _set_refgen(self,refgen,filter=True):
        '''
            Sets the current refgen. Its complicated.
        '''
        # Keep a record of parent refgen
        self._global('parent_refgen',refgen.name)
        # Filter down to only genes in 
        if filter:
            refgen = refgen.filtered_refgen(
                'Filtered{}'.format(self.name),
                'Filtered Refgen',
                refgen,
                self.genes(),
            )
        # remember to set for current instance
        self._global('refgen',refgen.name)
        self.refgen = refgen

    @property 
    def _cmap(self):
        '''
            Used for the heatmap function. Retruns a matplotlib cmap which is yellow/blue
        '''
        heatmapdict = {
            'red': ((0.0, 1.0, 1.0),
                    (0.5, 1.0, 1.0),
                    (1.0, 0.0, 0.0)),
            'green':((0.0, 1.0, 1.0),
                    (0.5, 1.0, 1.0),
                    (1.0, 0.0, 0.0)),
            'blue': ((0.0, 0.0, 0.0),
                    (0.5, 1.0, 1.0),
                    (1.0, 1.0, 1.0))}
        heatmap_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',heatmapdict,256)
        return heatmap_cmap


    ''' ------------------------------------------------------------------------------------------
            Class Methods
    '''

    @classmethod
    def create(cls,name,description,refgen):
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
        self = super().create(name,description,type='Expr')
        # Create appropriate HDF5 tables
        self.hdf5['expr'] = pd.DataFrame()
        self.hdf5['raw_expr'] = pd.DataFrame()
        self.hdf5['gene_qc_status'] = pd.DataFrame()
        self.hdf5['accession_qc_status'] = pd.DataFrame()
        self.hdf5['coex'] = pd.DataFrame()
        self.hdf5['degree'] = pd.DataFrame()
        # Delete existing datasets 
        self._set_refgen(refgen,filter=False)
        return self

    @classmethod
    def from_table(cls,filename,name,description,refgen,rawtype=None,
            sep='\t',normalize=True,quality_control=True,**kwargs):
        ''' 
            Create a Expr instance from a file containing raw expression data.
            For instance FPKM or results from a microarray experiment. This
            is a convenience method which reads the table in to a pandas DataFrame
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
                This is noted here to reinforce the impotance of the rawtype passed to 
                camoco.Expr.from_DataFrame. See docs there for more information.
            sep : str (default: \t)
                Column delimiter for the data in filename path
            normalize : bool (Default: True)
                Specifies whether or not to normalize the data so raw expression values
                lie within a log space. This is best practices for generating interpretable 
                expression analyses. See Expr._normalize method for more information.
                info.
            quality_control : bool (Default: True)
                A flag which specifies whether or not to perform QC. Parameters for QC are passed
                in using the **kwargs arguments. For default parameters and options 
                see Expr._quality_control.
            **kwargs : key value pairs
                additional parameters passed to subsequent methods. (see Expr.from_DataFrame)

            Returns
            -------
            An Expr instance

        '''
        tbl = pd.read_table(filename,sep=sep)
        return cls.from_DataFrame(tbl,name,description,refgen,rawtype=rawtype,**kwargs)

    @classmethod
    def from_DataFrame(cls,tbl,name,description,refgen,rawtype=None,
        normalize=True,norm_method=None,quantile=False,quality_control=True,**kwargs):
        ''' 
            Creates an Expr instance from a pandas DataFrame. Expects that the DataFrame
            index is gene names and the column names are accessions (i.e. experiments).
            This is the preferred method for creating an Expr instance, in other words, 
            other classmethods transform their data so they can call this method.

            Parameters
            ----------
            tbl : pandas.DataFrame
                a DataFrame containing expression data. Assumes index is the genes and
                columns is the accessions (experiment names)
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
                Specifies the fundamental datatype used to measure expression. During importation
                of the raw expression data, this value is used to make decisions in converting data
                to log-space. 
            normalize : bool (Default: True)
                Specifies whether or not to normalize the data so raw expression values
                lie within a log space. This is best practices for generating interpretable 
                expression analyses. See Expr._normalize method for more information.
                info.
            norm_method : None OR python function
                If rawtype is NOT RNASEQ or MICROARRY AND normalize is still True, the normalization
                method for the raw expression values needs to be passed in. This is for extreme customization
                situations. 
            quantile : bool (Default : False)
                Specifies whether or not to perform quantile normalization on import. 
            quality_control : bool (Default: True)
                A flag which specifies whether or not to perform QC. Parameters for QC are passed
                in using the **kwargs arguments. For default parameters and options 
                see Expr._quality_control.
            **kwargs : key value pairs
                additional parameters passed to subsequent methods. (see Expr.from_DataFrame)

            Returns
            -------
            An Expr instance

        '''
        # we are all pandas on the inside O.O
        self = cls.create(name,description,refgen)
        self._reset(raw=True)
        if rawtype is None:
            self.log(('WARNING: not passing in a rawtype makes'
                      ' downstream normalization hard...'))
            rawtype = ''
        self._global('rawtype',rawtype)
        # put raw values into the database
        self.log('Importing Raw Expression Values')
        self._update_values(tbl,'Raw'+rawtype,raw=True)
        if quality_control:
            self.log('Performing Quality Control on genes')
            self._quality_control(**kwargs)
        assert self.anynancol() == False
        if normalize:
            self.log('Performing Raw Expression Normalization')
            self._normalize(**kwargs)
            assert self.anynancol() == False
        if quantile:
            self.log('Performing Quantile Gene Normalization')
            self._quantile()
            assert self.anynancol() == False
        self.log('Filtering refgen: {}',refgen.name)
        self._set_refgen(refgen,filter=True)
        return self


    '''
        Unimplemented ---------------------------------------------------------------------------------
    '''
       
    def heatmap(self,genes=None,accessions=None,filename=None,figsize=(16,16), maskNaNs=True, 
        cluster_x=True, cluster_y=True,cluster_method="euclidian", title=None, zscore=True,raw=False, 
        heatmap_unit_label='Expression Z Score',png_encode=False):
        ''' 
            Draw clustered heatmaps of an expression matrix
        '''
        from matplotlib import rcParams
        rcParams.update({'figure.autolayout': True})
        dm = self._expr(genes=genes,accessions=accessions,zscore=zscore,raw=raw).T
        D = np.array(dm)
        row_labels = dm.index
        col_labels = dm.columns
        f = plt.figure(figsize=figsize,facecolor='white')
        # add matrix plot
        axmatrix = f.add_axes([0.3, 0.1, 0.5, 0.6])
        def masked_corr(x,y):
            mask = np.logical_and(np.isfinite(x),np.isfinite(y)) 
            if cluster_method == "euclidean":
                return euclidean(x[mask],y[mask])
            else:
                return pearsonr(x[mask],y[mask])[1]
        # add first dendrogram
        if cluster_y and len(dm.index) > 1:
            # calculate the squareform of the distance matrix
            D1 = squareform(pdist(D, masked_corr))
            ax1 = f.add_axes([0.09, 0.1, 0.2, 0.6])
            ax1.set_frame_on(False)
            Y = linkage(D1, method='complete')
            Z1 = dendrogram(Y, orientation='right')
            row_labels = row_labels[Z1['leaves']]
            D = D[Z1['leaves'], :]
            ax1.set_xticks([])
            ax1.set_yticks([])
        # add second dendrogram
        if cluster_x and len(dm.columns) > 1:
            D2 = squareform(pdist(D.T, masked_corr))
            ax2 = f.add_axes([0.3, 0.71, 0.5, 0.2])
            ax2.set_frame_on(False)
            Y = linkage(D2, method='complete')
            Z2 = dendrogram(Y)
            D = D[:, Z2['leaves']]
            col_labels = col_labels[Z2['leaves']]
            ax2.set_xticks([])
            ax2.set_yticks([])
        if title:
            plt.title(title)
        vmax = max(np.nanmin(abs(D)),np.nanmax(abs(D)))
        vmin = vmax*-1
        self.log("Extremem Values: {}",vmax)
        # Handle NaNs
        if maskNaNs:
            nan_mask = np.ma.array(D,mask=np.isnan(D))
            cmap = self._cmap
            cmap.set_bad('grey',1.0)
        else:
            cmap = self.__cmap
        im = axmatrix.matshow(D,aspect='auto',cmap=cmap,vmax=vmax,vmin=vmin)
        # Handle Axis Labels
        axmatrix.set_xticks(np.arange(D.shape[1]))
        axmatrix.xaxis.tick_bottom()
        axmatrix.tick_params(axis='x',labelsize='xx-small')
        axmatrix.set_xticklabels(col_labels,rotation=90,ha='center')
        axmatrix.yaxis.tick_right()
        axmatrix.set_yticks(np.arange(D.shape[0]))
        axmatrix.set_yticklabels(row_labels)
        axmatrix.tick_params(axis='y',labelsize='x-small')
        plt.gcf().subplots_adjust(right=0.15)
        # Add color bar
        axColorBar = f.add_axes([0.09,0.75,0.2,0.05])
        f.colorbar(im,orientation='horizontal',cax=axColorBar,
            ticks=np.arange(np.ceil(vmin),np.ceil(vmax),int((vmax-vmin)/2))
        )
        plt.title(heatmap_unit_label)
        if filename:
            plt.savefig(filename)
            plt.close()
        if png_encode is True:
            imgdata = io.BytesIO()
            plt.savefig(imgdata)
            return base64.encodebytes(imgdata.getvalue()).decode()
            
        return pd.DataFrame(
            data=D,
            index=row_labels,
            columns=col_labels
        )

    def plot_value_hist(self,groupby='accession',raw=False,bins=50,figsize=(16,16),title='',log=False):
        ''' Plots Value histograms on one of the expression matrix axis'''
        for group,df in self._expr(long=True,raw=raw).groupby(groupby):
            self.log('Plotting values for {}',group)
            plt.clf()
            plt.hist(
                list(filter(lambda x: not np.isnan(x),df.value)),
                bins=bins,
                log=log
            )
            plt.title(group+title)
            plt.xlabel('Expression')
            plt.ylabel('Frequency')
            plt.savefig(
                filename="{}{}_VALUES.png".format(group,title),
                figsize=figsize
            )
