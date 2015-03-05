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

class Expr(Camoco):
    ''' A representation of gene expression. '''
    def __init__(self,name,description=None,basedir='~/.camoco'):
        super().__init__(name=name,description=description,type='Expr',basedir=basedir) 
        try:
            self.log('Loading RefGen')
            self.log('Loading Expr table')
            self.hdf5 = self._hdf5(name)
            try:
                self.expr = self.hdf5['expr']
            except KeyError as e:
                self.log('{} is empty: ({})',name,e)
                self.expr = pd.DataFrame()
            self.log('Building Expr Index')
            self._expr_index = defaultdict(lambda: None,{gene:index for index,gene in enumerate(self.expr.index)}) 
            self.refgen = RefGen(self.refgen)
        except NameError as e:
            self.log.warn('Refgen for {} not available, must be reset!',self.name)

    def __contains__(self,obj):
        if obj in self.expr.index:
            return True
        if obj in self.expr.columns:
            return True
        try:
            if obj.id in self.expr.index:
                return True
        except AttributeError as e:
            pass
        return False

    def __repr__(self):
        return ""

    def __str__(self):
        pass

    def num_genes(self):
        return len(self.expr.index)

    def num_accessions(self):
        return len(self.expr.columns)

    def shape(self):
        return self.expr.shape

    def zscore(self):
        pass

    def accessions(self):
        return self.expr().columns

    def genes(self,enumerated=True,raw=False):
        # Returns a list of distinct genes 
        if raw is False:
            return self.refgen.from_ids(self.expr.index,enumerated=enumerated)
        else:
            return self.refgen.from_ids(self.hdf5['raw_expr'].index,enumerated=enumerated)

    def expr_profile(self,gene):
        '''
            return the expression profile for a gene
        '''
        return self.expr.loc[gene.id]

    def is_normalized(self,max_val=None,raw=False):
        if max_val is not None:
            max_val = max_val # Use the user defined max val
        elif self.rawtype.upper() == 'RNASEQ':
            max_val = 1100 
        elif self.rawtype.upper() == 'MICROARRAY':
            max_val = 100 
        return self.expr.apply(lambda col: np.nanmax(col.values) < max_val ,axis=0)
    
    def anynancol(self):
        ''' A gut check method to make sure none of the expression columns
            got turned into all nans. Because apparently that is a problem.'''
        return any(self.expr.apply(lambda col: all(np.isnan(col)),axis=0))

    '''
        Internal Methods --------------------------------------------------------------------------------
    '''

    def _update_values(self,df,transform_name,raw=False):
        '''
            updates the 'expression' table values with values from df.
            Requires a transformation name for the log. 
            Option to overwrite raw table or working table.
        '''
        # update the transformation log
        self._transformation_log(transform_name)
        table = 'raw_expr' if raw else 'expr'
        # Sort the table by genes
        df = df.sort()
        try:
            self.hdf5[table] = self.expr = df
            self.hdf5.flush(fsync=True)
        except Exception as e:
            self.log('Unable to update expression table values: {}',e)
            raise

    def _transformation_log(self,transform=None):
        if transform is None:
            return self._global('transformation_log')
        elif transform == 'reset' or self._global('transformation_log') is None:
            self._global('transformation_log','raw')
        else:
            self._global('transformation_log',self._global('transformation_log') + '->' + str(transform))
            self.log('Transformation Log: {}',self._global('transformation_log'))
 
    def _reset(self,raw=False):
        ''' resets the expression values to their raw state undoing any normalizations '''
        if raw:
            # kill the raw table too
            self.log('Resetting raw expression data')
            self.hdf5['raw_expr'] = pd.DataFrame()
        self.log('Resetting expression data')
        self.hdf5['expr'] = self.expr = pd.DataFrame()
        self._transformation_log('reset')


    def _normalize(self,method=None,is_raw=None,max_val=None,**kwargs):
        ''' evaluates qc expression data and re-enters 
            normaized data into database '''
        self.log('------------ Normalizing')
        if all(self.is_normalized(max_val=max_val)):
            self.log("Dataset already normalized")
            self._transformation_log('DetectedPreNormalized')
        elif any(self.is_normalized(max_val=max_val)):
            # Something fucked up is happending
            raise TypeError('Attempting normalization on already normalized dataset. Consider passing a max_val < {} if Im wrong.'.format(min(self.max_values())))
        else:
            df = self.expr
            if method is not None:
                method = method
            elif self.rawtype.upper() == 'RNASEQ':
                method = np.arcsinh
            elif self.rawtype.upper() == 'MICROARRAY':
                method = np.log2
            else:
                raise ValueError('Could not guess correct normalization for {}, pass in function through method argument.'.format(self.rawtype))
            # apply the normalization to each column (accession)
            df = df.apply(lambda col: method(col),axis=0)
            # update values
            self._update_values(df,method.__name__)

    def _quality_control(self,min_expr=1,max_gene_missing_data=0.2,min_single_sample_expr=5, \
        max_accession_missing_data=0.5,membership=None,dry_run=False,**kwargs):
        ''' Sets quality control flag for all values in expression table '''        
        self.log('------------Quality Control')
        df = self.hdf5['raw_expr']
        # remember how we set the flags
        self._global('qc_min_expr',min_expr)
        self._global('qc_max_gene_missing_data',max_gene_missing_data)
        self._global('qc_min_single_sample_expr',min_single_sample_expr)
        self._global('qc_max_accession_missing_data',max_accession_missing_data)
        self._global('qc_membership',str(membership))
        # retrieve raw data as a data frame
        self.log('Raw Starting set: {} genes {} accessions'.format(len(df.index),len(df.columns)))
        # ------ Membership test
        if not membership:
            membership = self.refgen
        self.log("Filtering out genes not in {}",membership)
        df = df[[x in membership for x in df.index]]
        self.log('Kept: {} genes {} accessions'.format(len(df.index),len(df.columns)))
        # -----------------------------------------
        # Set minimum FPKM threshold
        self.log("Filtering out expression values lower than {}",min_expr)
        df_flt = df.copy()
        df_flt[df < min_expr] = np.nan
        df = df_flt
        self.log('Kept: {} genes {} accessions'.format(len(df.index),len(df.columns)))
        # -----------------------------------------
        # filter out genes with too much missing data
        self.log("Filtering out genes with > {} missing data",max_gene_missing_data)
        df = df.loc[df.apply(lambda x : ((sum(np.isnan(x))) < len(x)*max_gene_missing_data),axis=1),:] 
        self.log('Kept: {} genes {} accessions'.format(len(df.index),len(df.columns)))
        # -----------------------------------------
        # filter out genes which do not meet a minimum expr threshold in at least one sample
        self.log("Filtering out genes which do not have one sample above {}",min_single_sample_expr)
        df = df[df.apply(lambda x: any(x >= min_single_sample_expr),axis=1)]
        self.log('Kept: {} genes {} accessions'.format(len(df.index),len(df.columns)))
        # -----------------------------------------
        # Filter out ACCESSIONS with too much missing data
        self.log("Filtering out accessions with > {} missing data",max_accession_missing_data)
        accession_mask = df.apply(lambda x : ((sum(np.isnan(x))) / len(x)),axis=0)
        for i,percent in enumerate(accession_mask):
            if percent > max_accession_missing_data:
                self.log("\tRemoved: {}: missing {} data",df.columns[i],percent*100)
        df = df.loc[:,np.logical_not(accession_mask > max_accession_missing_data)] 
        self.log('Kept: {} genes {} accessions'.format(len(df.index),len(df.columns)))
        if dry_run:
            # If dry run, take first 100 rows of QC
            self.log.warn("Dry Run")
            df = df.iloc[0:5000,:]
        self._update_values(df,'quality_control')

    def _quantile(self):
        ''' Perform quantile normalization across each accession. 
            Each accessions gene expression values are replaced with 
            ranked gene averages.'''
        # get gene by accession matrix
        self.log('------------Quantile')
        expr = self.expr
        self.log('Ranking data')
        # assign ranks by accession
        expr_ranks = expr.rank(axis=0,method='dense')
        # normalize rank to be percentage
        import pdb; pdb.set_trace()
        expr_ranks = expr_ranks.apply(lambda col: col/np.nanmax(col.values), axis=0)
        # we need to know the number of non-nans so we can correct for their ranks later
        self.log('Sorting ranked data')
        # assign accession values by order
        # NOTE all NANs get sorted here ...
        expr_sort = expr.apply(np.sort,axis=0)
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
        return co.RefGen(self._global['parent_refgen'])

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


    ''' ------------------------------------------------------------------------------------------
            Class Methods
    '''
    @classmethod
    def from_table(cls,filename,name,description,refgen,rawtype=None,basedir='~/.camoco'
        ,sep='\t',normalize=True,quality_control=True,**kwargs):
        ''' Returns a Expr instance read in from a table file '''
        tbl = pd.read_table(filename,sep=sep)
        return cls.from_DataFrame(tbl,name,description,refgen,rawtype=rawtype,**kwargs)

    @classmethod
    def from_DataFrame(cls,tbl,name,description,refgen,rawtype=None,basedir='~/.camoco',
        normalize=True,quantile=True,quality_control=True,**kwargs):
        ''' Imports a Expr instance based on a pandas table (genes as rows and accessions as cols)'''
        # we are all pandas on the inside O.O
        self = cls(name=name,description=description,basedir=basedir)
        self._set_refgen(refgen,filter=False)
        self._reset(raw=True)
        if rawtype is None:
            self.log('WARNING: not passing in a rawtype makes downstream normalization hard...')
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

    @property 
    def cmap(self):
        '''
            Used for the heatmap function. Retruns a cmap which is yellow/blue
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



    '''
        Unimplemented ---------------------------------------------------------------------------------
    '''
       
    def heatmap(self,genes=None,accessions=None,filename=None,figsize=(16,16), maskNaNs=True, 
        cluster_x=True, cluster_y=True,cluster_method="euclidian", title=None, zscore=True,raw=False, 
        heatmap_unit_label='Expression Z Score',png_encode=False):
        ''' Draw clustered heatmaps of an expression matrix'''
        from matplotlib import rcParams
        rcParams.update({'figure.autolayout': True})
        dm = self.expr(genes=genes,accessions=accessions,zscore=zscore,raw=raw).T
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
            cmap = self.cmap
            cmap.set_bad('grey',1.0)
        else:
            cmap = self.cmap
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
        for group,df in self.expr(long=True,raw=raw).groupby(groupby):
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


