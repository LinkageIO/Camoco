#! /usr/bin/env python3
from camoco.Camoco import Camoco
from camoco.RefGen import RefGen
from camoco.Tools import memoize
from scipy.spatial.distance import pdist, squareform, euclidean
from scipy.stats import hypergeom,pearsonr
from scipy.cluster.hierarchy import linkage, dendrogram

import matplotlib
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import io               

class Expr(Camoco):
    ''' A representation of gene expression. '''
    def __init__(self,name,description=None,basedir='~/.camoco'):
        super().__init__(name=name,description=description,type='Expr',basedir=basedir) 
        self._create_tables()
        try:
            self.refgen = RefGen(self.refgen)
        except NameError as e:
            self.log.warn('Refgen for {} not available, must be reset!',self.name)

    def __repr__(self):
        return ""

    def __str__(self):
        pass

    def __contains__(self,obj):
        try:
            # Object is a gene
            return True if obj in self.genes() else False
        except Exception as e:
            pass
        try:
            # Can be a string object
            return True if obj in (x.id for x in self.genes()) else False
        except Exception as e:
            pass
        self.log("object '{}' is not correct type to test for membership in {}",obj,self.name)

    def summary(self):
        pass

    @memoize
    def accessions(self):
        return [x[0] for x in self.db.cursor().execute('SELECT DISTINCT(accession) FROM expression')]

    @property
    @memoize
    def num_accessions(self):
        return self.db.cursor().execute("SELECT COUNT(DISTINCT(accession)) FROM expression").fetchone()[0]

    # Dont remove this! it is required for filter_refgen and is rather slow if not memoized
    @memoize
    def genes(self,enumerated=True,raw=False):
        # Returns a list of distinct genes 
        table = 'raw_expression' if raw else 'expression'
        return self.refgen.from_ids([
            x[0] for x in self.db.cursor().execute('SELECT DISTINCT(gene) FROM {}'.format(table))
        ],enumerated=enumerated)

    @memoize
    def num_genes(self):
        return self.db.cursor().execute('SELECT COUNT(DISTINCT(gene)) FROM expression').fetchone()[0]

    def del_accession(self,name):
        cur = self.db.cursor()
        try:
            cur.exectue('START TRANSACTION') 
            cur.execute('DELETE FROM expression WHERE accession = ?',(name))
            cur.exectue('END TRANSACTION') 
        except Exception as e:
            cur.execute('ROLLBACK') 

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
        self.reset(raw=True)
        self._set_refgen(refgen)
        if rawtype is None:
            self.log('WARNING: not passing in a rawtype makes downstream normalization hard...')
            rawtype = ''
        self._global('rawtype',rawtype)
        # put raw values into the database
        self.log('Importing Raw Expression Values')
        self.update_values(tbl,'Raw'+rawtype,raw=True)
        if quality_control:
            self.log('Performing Quality Control on genes')
            self.quality_control(**kwargs)
        assert self.anynancol() == False
        if normalize:
            self.log('Performing Raw Expression Normalization')
            self.normalize(**kwargs)
            assert self.anynancol() == False
        if quantile:
            self.log('Performing Quantile Gene Normalization')
            self.quantile()
            assert self.anynancol() == False
        return self

    def reset(self,raw=False):
        ''' resets the expression values to their raw state undoing any normalizations '''
        cur = self.db.cursor()
        cur.execute('BEGIN TRANSACTION')
        if raw:
            # kill the raw table too
            self.log('Resetting raw expression data')
            cur.execute('DELETE FROM raw_expression;') 
        self.log('Resetting expression data')
        cur.execute('DELETE FROM expression;') 
        cur.execute('END TRANSACTION')
        if raw and cur.execute('SELECT COUNT(*) FROM raw_expression').fetchone()[0] != 0:
            raise ValueError("Raw Expression Table NOT EMPTY!")
        if cur.execute('SELECT COUNT(*) FROM expression').fetchone()[0] != 0:
            raise ValueError("Expression Table NOT EMPTY!")
        self._set_refgen(None)
        self.transformation_log('reset')

    def update_values(self,df,transform_name,raw=False):
        ''' updates the 'expression' table values with values from df.
            Requires a transformation name for the logs. 
            Option to overwrite raw table or working table.
            '''
        # update the transformation log
        self.transformation_log(transform_name)
        table = 'raw_expression' if raw else 'expression'
        cur = self.db.cursor()
        # Sort the table by genes
        df = df.sort()
        try:
            # lets be safe about this
            cur.execute('BEGIN TRANSACTION')
            self.log('Removing old values from {}...',table)
            cur.execute('DELETE FROM {};'.format(table)) 
            cur.execute('END TRANSACTION')
            cur.execute('BEGIN TRANSACTION')
            self.log('Importing values, shape: {} {}',df.shape[0],df.shape[1])
            cur.executemany('''
                INSERT OR REPLACE INTO {} (gene, accession, value) VALUES (?,?,?)'''.format(table),
                [(gene.upper(),accession, float(value)) for (accession,gene),value in df.unstack().iteritems()] 
                #          ^             ^ All indices are uppercase in database 
            )
            self.log('Imported {} values',df.shape[0]*df.shape[1])
            cur.execute('END TRANSACTION')
        except Exception as e:
            cur.execute('ROLLBACK') 
            self.log('Unable to update expression table values: {}',e)
            raise

    def max_values(self,group_by='accession',raw=False):
        return self.expr(raw=raw,long=True).groupby(group_by).apply(
            lambda x: np.nanmax(x.value.values # I named my column value just like pandas ....
        ))


    def is_normalized(self,max_val=None,raw=False):
        if max_val is not None:
            max_val = max_val # Use the user defined max val
        elif self.rawtype.upper() == 'RNASEQ':
            max_val = 1100 
        elif self.rawtype.upper() == 'MICROARRAY':
            max_val = 100 
        return self.expr(raw=raw).apply(lambda col: np.nanmax(col.values) < max_val ,axis=0)
    
    def anynancol(self):
        ''' A gut check method to make sure none of the expression columns
            got turned into all nans. Because apparently that is a problem.'''
        return any(self.expr().apply(lambda col: all(np.isnan(col)),axis=0))

    def normalize(self,method=None,is_raw=None,max_val=None,**kwargs):
        ''' evaluates qc expression data and re-enters 
            normaized data into database '''
        self.log('------------ Normalizing')
        if all(self.is_normalized(max_val=max_val)):
            self.log("Dataset already normalized")
            self.transformation_log('DetectedPreNormalized')
        elif any(self.is_normalized(max_val=max_val)):
            # Something fucked up is happending
            raise TypeError('Attempting normalization on already normalized dataset. Consider passing a max_val < {} if Im wrong.'.format(min(self.max_values())))
        else:
            df = self.expr(raw=False,long=False)
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
            self.update_values(df,method.__name__)

    def quality_control(self,min_expr=1,max_gene_missing_data=0.2,min_single_sample_expr=5, \
        max_accession_missing_data=0.5,membership=None,dry_run=False,**kwargs):
        ''' Sets quality control flag for all values in expression table '''        
        self.log('------------Quality Control')
        df = self.expr(raw=True,long=False)
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
        self.update_values(df,'quality_control')

    def expr(self,genes=None,accessions=None,long=False,raw=False,zscore=False):
        ''' returns expression data '''
        # Set up dynamic query
        tbl = 'raw_expression' if raw else 'expression'
        gene_filter = " WHERE gene in ('{}')".format("','".join([x.id.upper() for x in genes])) if genes is not None else ""
        accession_filter = " AND accession in ('{}')".format("','".join(accessions)) if accessions is not None else ""
        query = 'SELECT * FROM {} {} {};'.format(tbl,gene_filter,accession_filter)
        # pull and create data frame
        df = pd.DataFrame(
            self.db.cursor().execute(query).fetchall(),
            columns = ['gene','accession','value']
        )
        if long is False:
            df = df.pivot('gene','accession','value')
        if zscore is True:
            if long is True:
                raise ValueError('If zscore is true, long must be false')
            else:
                df = df.apply(lambda x: (x-x.mean())/(x.std()),axis=1)
        return df

    def quantile(self):
        ''' Perform quantile normalization across each accession. 
            Each accessions gene expression values are replaced with 
            ranked gene averages.'''
        # get gene by accession matrix
        self.log('------------Quantile')
        expr = self.expr()
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
        self.update_values(expr,'quantile')

    def transformation_log(self,transform=None):
        if transform is None:
            return self._global('transformation_log')
        elif transform == 'reset' or self._global('transformation_log') is None:
            self._global('transformation_log','raw')
        else:
            self._global('transformation_log',self._global('transformation_log') + '->' + str(transform))
            self.log('Transformation Log: {}',self._global('transformation_log'))
        
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

    @property 
    def cmap(self):
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

    def _set_refgen(self,refgen):
        # update database for future
        if refgen is not None:
            self._global('refgen',refgen.name)
        else:
            self._global('refgen',None)
        # remember to set for current instance
        self.refgen = refgen

    def _filter_refgen(self):
        ''' Filter the refgen to only contain genes available in COB. Only do this after the expression table
            has been populated!!'''
        self.log("Filtering custom refgen for {} genes in {}",self.num_genes(),self.name)
        filtered_refgen = RefGen.from_RefGen(
            "{}RefGen".format(self.name),
            "RefGen for {} filtered from {}".format(self.name,self.refgen.name),
            self.refgen,
            gene_filter=self,
            basedir = self.basedir
        )
        self._set_refgen(filtered_refgen)


    def _create_tables(self):
        cur = self.db.cursor()
        cur.execute(''' 
            CREATE TABLE IF NOT EXISTS raw_expression (
                gene TEXT,
                accession TEXT,
                value REAL,
                UNIQUE (gene,accession) ON CONFLICT FAIL
            );
            CREATE TABLE IF NOT EXISTS expression (
                gene TEXT,
                accession TEXT,
                value REAL,
                UNIQUE(gene,accession) ON CONFLICT FAIL
            );
            CREATE INDEX IF NOT EXISTS genes_ind ON expression(gene);
            CREATE INDEX IF NOT EXISTS accessions_ind ON expression(accession);
        ''')
