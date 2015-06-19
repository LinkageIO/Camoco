#!/usr/bin/python3
import pyximport; pyximport.install()
import camoco.PCCUP as PCCUP

from camoco.Camoco import Camoco
from camoco.RefGen import RefGen
from camoco.Locus import Locus,Gene
from camoco.Expr import Expr
from camoco.Tools import memoize

from numpy import matrix,arcsinh,tanh
from collections import defaultdict,Counter
from itertools import chain
from subprocess import Popen, PIPE
from scipy.spatial.distance import squareform
from scipy.misc import comb
from scipy.stats import norm
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std

import pandas as pd
import numpy as np
import itertools
import matplotlib.pylab as plt

from scipy.stats import pearsonr

class COB(Expr):
    def __init__(self,name):
        super().__init__(name=name)
        #self.hdf5 = self._hdf5(name)
        try:
            self.log('Loading coex table')
            self.coex = self.hdf5['coex']
        except KeyError as e:
            self.log("{} is empty ({})",name,e)
        try:
            self.log('Loading Global Degree')
            self.degree = self.hdf5['degree']
        except KeyError as e:
            self.log("{} is empty ({})",name,e)

    def __repr__(self):
        return '<COB: {}>'.format(self.name)

    def __str__(self):
        return self.__repr__()

    def summary(self):
        print( '''
            COB Dataset: {} - {} - {}
                Desc: {}
                RawType: {}
                TransformationLog: {}
                Num Genes: {}
                Num Accessions: {}
                Edge FDR: {}
        '''.format(
                self.name, self.organism, self.build,
                self.description,
                self.rawtype,
                self._transformation_log(),
                self.num_genes(),
                self.num_accessions(),
                self.edge_FDR
        ))

    @property
    @memoize
    def edge_FDR(self):
        # get the percent of significant edges
        num_sig = sum(self.coex['significant'])/len(self.coex)
        # calulate the number expected
        num_exp = 1-norm.cdf(int(self._global('significance_threshold')))
        # FDR is the percentage expected over the percentage found
        return num_exp/num_sig

    def neighbors(self,gene,sig_only=True):
        '''
            Returns a DataFrame containing the neighbors for gene.
    
            Parameters
            ----------
            gene : co.Locus
                The gene for which to extract neighbors
    
            Returns
            -------
            A DataFrame containing edges
        '''
        gene_id = self._expr_index[gene.id]
        neighbor_indices = PCCUP.coex_neighbors(gene_id,self.num_genes())
        edges = self.coex.iloc[neighbor_indices]
        if sig_only:
            return edges[edges.significant == 1]
        else:
            return edges

    def coexpression(self,gene_a,gene_b):
        ''' 
            Returns a coexpression z-score between two genes. This
            is the pearson correlation coefficient of the two genes'
            expression profiles across the accessions (experiments).
            This value is pulled from the 

            Parameters
            ----------
            gene_a : camoco.Locus
                The first gene
            gene_b : camoco.Locus
                The second gene
        
            Returns
            -------
            Coexpression Z-Score 

        '''
        # Grab the indices in the original expression matrix
        ids = np.array([self._expr_index[gene_a.id],self._expr_index[gene_b.id]])
        # We need the number of genes
        num_genes = self.num_genes()
        index = PCCUP.coex_index(ids,num_genes)[0]
        return self.coex.iloc[index]

    def subnetwork(self,gene_list=None,sig_only=True,min_distance=100000,
        filter_missing_gene_ids=True):
        '''
            Input: a gene list (passing None gives you all genes)
            Output: a dataframe containing all edges EXCLUSIVELY between genes
                within list
        '''
        if gene_list is None:
            df = self.coex
        else:
            ids = np.array([self._expr_index[x.id] for x in gene_list])
            if filter_missing_gene_ids:
                # filter out the Nones 
                ids = np.array(list(filter(None,ids)))
            num_genes = self.num_genes()
            # Grab the coexpression indices for the genes
            indices = PCCUP.coex_index(ids,num_genes)
            df = self.coex.iloc[indices]
        if min_distance:
            df = df[df.distance >= min_distance]
        if sig_only:
            df = df[df.significant == 1]
        return df

    def trans_locus_density(self,locus_list,return_mean=True):
        '''
            Calculates the 
        '''
        pass

    def density(self,gene_list,return_mean=True,min_distance=50000):
        ''' 
            calculates the denisty of the non-thresholded network amongst genes
            not within a certain distance of each other. This corrects for
            cis regulatory elements increasing noise in coexpression network 
        '''
        # filter for only genes within network
        edges = self.subnetwork(gene_list,
            min_distance=min_distance,sig_only=False
        )
        if len(edges) == 0:
            return np.nan
        if len(edges) == 1:
            return edges.score[0]
        if return_mean:
            return np.nanmean(edges.score)/(1/np.sqrt(len(edges)))
            # old code worth a look
            # return ((np.nanmean(edges.score)/((np.nanstd(edges.score))/np.sqrt(len(edges)))),
            #        (np.nanmean(edges.score)/(1/np.sqrt(len(edges)))),
            #        (np.nanmedian(edges.score)/((np.nanstd(edges.score))/np.sqrt(len(edges)))))
        else:
            return edges

    def plot_scores(self,filename=None,pcc=True,bins=50):
        ''' 
            Plot the histogram of PCCs.
        '''
        if filename is None:
            filename = self.name+'.png'
        plt.clf()
        # grab the scores only and put in a np array to save space (pandas DF was HUGE)
        scores = (self.coex.score.values + float(self._global('pcc_mean'))) * float(self._global('pcc_std'))
        if pcc:
            self.log('Transforming scores')
            # Transform Z-scores to pcc scores (inverse fisher transform)
            scores = tanh(scores)
        plt.hist(scores,bins=bins)
        plt.xlabel('PCC')
        plt.ylabel('Freq')
        plt.savefig(filename) 

    def to_dat(self,gene_list=None,filename=None,sig_only=False,min_distance=0):
        '''
            Outputs a .DAT file (see Sleipnir library)
        '''
        if filename is None:
            filename = self.name + '.dat'
        with open(filename, 'w') as OUT:
            self.log("Creating .dat file")
            self.subnetwork(
                gene_list,sig_only=sig_only,min_distance=min_distance
            )['score'].to_csv(OUT,sep='\t')
            self.log('Done')

    def mcl(self,gene_list=None,I=2.0,scheme=7,min_distance=100000,min_cluster_size=0,max_cluster_size=10e10):
        ''' 
            A *very* thin wrapper to the MCL program. The MCL program must
            be accessible by a subprocess (i.e. by the shell).
            Returns clusters (as list) as designated by MCL. 
            
            Parameters
            ----------
            gene_list : a gene iterable
                These are the genes which will be clustered
            I : float (default: 2.0)
                This is the inflation parameter passed into mcl.
            scheme : int in 1:7
                MCL accepts parameter schemes. See mcl docs for more details
            min_distance : int (default: 100000)
                The minimum distance between genes for which to consider co-expression
                interactions. This filters out cis edges. 
            min_cluster_size : int (default: 0)
                The minimum cluster size to return. Filter out clusters smaller
                than this.
            max_cluster_size : float (default: 10e10)
                The maximum cluster size to return. Filter out clusters larger
                than this. 

            Returns
            -------
            A list clusters containing a lists of genes within each cluster
        '''
        # output dat to tmpfile
        tmp = self._tmpfile()
        self.to_dat(filename=tmp.name,gene_list=gene_list,min_distance=min_distance,sig_only=True)
        # build the mcl command 
        cmd = "mcl {} --abc -scheme {} -I {} -o -".format(tmp.name,scheme,I)
        self.log("running MCL: {}",cmd)
        try:
            p = Popen(cmd, stdout=PIPE, stderr=self.log_file, shell=True)
            self.log('waiting for MCL to finish...')
            sout = p.communicate()[0]
            p.wait()
            self.log('...Done')
            if p.returncode==0:
                # Filter out cluters who are smaller than the min size
                return list(filter(lambda x: len(x) > min_cluster_size and len(x) < max_cluster_size,
                    # Generate ids from the refgen
                    [ self.refgen.from_ids([gene.decode('utf-8') for gene in line.split()]) for line in sout.splitlines() ]
                ))
            else:
                self.log( "MCL failed: return code: {}".format(p.returncode))
        except FileNotFoundError as e:
            self.log('Could not find MCL in PATH. Make sure its installed and shell accessible as "mcl".')

    def local_degree(self,genes):
        '''
            Returns the local degree of a list of genes
        '''
        local_degree = pd.DataFrame(
            list(Counter(chain(*self.subnetwork(genes,sig_only=True).index.get_values())).items()),
            columns=['Gene','Degree']
        ).set_index('Gene')
        # We need to find genes not in the subnetwork and add them as degree 0
        degree_zero_genes = pd.DataFrame( # The code below is optimized
            [(gene.id,0) for gene in genes if gene.id not in local_degree.index],
            columns=['Gene','Degree']        
        ).set_index('Gene')

        return pd.concat([local_degree,degree_zero_genes])

    def global_degree(self,genes):
        '''
            Returns the global degree of a list of genes
        '''
        try:
            if isinstance(genes,Locus):
                return self.degree.ix[genes.id].Degree
            else:
                return self.degree.ix[[x.id for x in genes]].fillna(0)
        except KeyError as e:
            return 0

    def locality(self, gene_list,bootstrap_name=None,include_regression=False):
        '''
            Computes the merged local vs global degree table
            
            Parameters
            ----------
            gene_list : iterable of camoco.Loci
                A list or equivalent of loci
            bootstrap_name : object (default: none)
                This will be added as a column. Useful for
                generating bootstraps of locality and keeping
                track of which one a row came from after catting
                multiple bootstraps together. 

            Returns
            -------
            A pandas DataFrame with local,global and residual columns
            based on linear regression of local on global degree. 

        '''
        self.log('Fetching Degree ... ')
        degree = self.local_degree(gene_list)
        # Add on the global degree
        degree['global'] = self.global_degree(self.refgen.from_ids(degree.index.values))['Degree']
        degree.columns = ['local','global']
        degree = degree.sort('global')
        if include_regression:
            # Add the regression lines
            ols = sm.OLS(degree['local'],degree['global']).fit()
            degree['resid'] = ols.resid
            degree['fitted'] = ols.fittedvalues
        if bootstrap_name is not None:
            degree['bootstrap_name'] = bootstrap_name
        return degree


    def plot_locality(self,gene_list,bootstraps=10,num_windows=100,sd_thresh=2):
        '''
            Make a fancy locality plot.
        '''
        # Generate a blank fig
        fig,ax = plt.subplots(figsize=(8,6)) 
        fig.hold(True)
        # Y axis is local degree (what we are TRYING to predict)
        degree = self.locality(gene_list).sort('global')
        ax.set_ylim(0,max(degree['local']))
        ax.set_xlim(0,max(degree['global']))
        if bootstraps > 0:
            bs = pd.concat(
                [self.locality(
                    self.refgen.bootstrap_candidate_genes(gene_list)
                ) for x in range(10)]
            ).sort('global')
            ax.set_ylim(0,max(bs['local']))
            ax.set_xlim(0,max(bs['global']))
            plt.plot(bs['global'],bs['local'],'ro',alpha=0.05,label='Bootstraps')
        # Plot the bootstraps and the empirical
        plt.plot(degree['global'],degree['local'],'bo',label='Empirical')
        emp_ols = sm.OLS(degree['local'],degree['global']).fit()
        ax.plot(degree['global'],emp_ols.fittedvalues,'k:',label='Empirical OLS')

        if bootstraps > 0:
            # Get the OLS
            bs_ols = sm.OLS(bs['local'],bs['global']).fit()
            bs['resid'] = bs_ols.resid
            bs['fitted'] = bs_ols.fittedvalues
            ax.plot(bs['global'],bs_ols.fittedvalues,'g--',label='bootstrap OLS')
            # Do lowess on the residuals
            # We only care about windows within the empirical part
            window_tick = len(bs)/num_windows
            bs['window'] = [int(x/window_tick) for x in range(len(bs))]
            # get std for each window
            win_std = bs.groupby('window').apply(lambda df: df['resid'].std()).to_dict()
            bs['std_envelope'] = [win_std[x] for x in bs.window.values]
            # Plot confidence intervals
            prstd, iv_l, iv_u = wls_prediction_std(bs_ols)           
            ax.plot(bs['global'], iv_u, 'g--',label='conf int.')
            ax.plot(bs['global'], iv_l, 'g--')
            # plot the  
            ax.plot(
                bs['global'],bs['fitted']+(sd_thresh*bs['std_envelope']),'r--'
                ,label='{} s.d. envelope'.format(sd_thresh)
            )
            ax.plot(bs['global'],bs['fitted']-(sd_thresh*bs['std_envelope']),'r--')
        ax.set_xlabel('Number Global Interactions')
        ax.set_ylabel('Number Local Interactions')
        legend = ax.legend(loc='best')
        return plt


    ''' ------------------------------------------------------------------------------------------
            Internal Methods
    '''
    def _coex_indices(self,ids):
        return 


    def _calculate_coexpression(self,significance_thresh=3):
        ''' 
            Generates pairwise PCCs for gene expression profiles in self._expr.
            Also calculates pairwise gene distance.
        '''
        # Start off with a fresh set of genes we can pass to functions
        tbl = pd.DataFrame(
            list(itertools.combinations(self._expr.index.values,2)),
            columns=['gene_a','gene_b']
        )
        # Now add coexpression data
        self.log("Calculating Coexpression")
        # Calculate the PCCs
        pccs = 1-PCCUP.pair_correlation(np.ascontiguousarray(self._expr.as_matrix()))
        # return the long form of the 
        assert len(pccs) == len(tbl)
        tbl['score'] = pccs
        # correlations of 1 dont transform well, they cause infinities
        tbl.loc[tbl['score'] == 1,'score'] = 0.99999999
        tbl.loc[tbl['score'] == -1,'score'] = -0.99999999
        # Perform fisher transform on PCCs
        tbl['score'] = np.arctanh(tbl['score'])
        # Sometimes, with certain datasets, the NaN mask overlap completely for the
        # two genes expression data making its PCC a nan. This affects the mean and std fro the gene.
        valid_scores = np.ma.masked_array(tbl['score'],np.isnan(tbl['score']))
        # Calculate Z Scores
        pcc_mean = valid_scores.mean()
        pcc_std = valid_scores.std()
        # Remember these so we can go back to PCCs
        self._global('pcc_mean',pcc_mean)
        self._global('pcc_std',pcc_std)
        tbl['score'] = (valid_scores-pcc_mean)/pcc_std
        # Assign significance
        self._global('significance_threshold',significance_thresh)
        tbl['significant'] = pd.Series(list(tbl['score'] >= significance_thresh),dtype='int_')
        self.log("Calculating Gene Distance")
        distances = self.refgen.pairwise_distance(gene_list=self.refgen.from_ids(self._expr.index))
        assert len(distances) == len(tbl)
        tbl['distance'] = distances
        # Reindex the table to match genes
        self.log('Indexing coex table')
        tbl.set_index(['gene_a','gene_b'],inplace=True)
        # put in the hdf5 store
        self._build_tables(tbl)
        self.log("Done")
        return self

    def _calculate_degree(self):
        try:
            self.log('Building Degree')
            self.hdf5['degree'] = pd.DataFrame(
                list(Counter(chain(*self.subnetwork(sig_only=True).index.get_values())).items()),
                columns=['Gene','Degree']
            ).set_index('Gene')
            self.hdf5.flush(fsync=True)
            self.degree = self.hdf5['degree']
        except Exception as e:
            self.log("Something bad happened:{}",e)
            raise

    def _build_tables(self,tbl):
        try:
            self.log("Building Database")
            ## HDF5 Store
            self.hdf5['coex'] = tbl
            self.log("Flushing Database")
            self.hdf5.flush(fsync=True)
            self.coex = self.hdf5['coex']
            self.log("Done")
        except Exception as e:
            self.log("Something bad happened:{}",e)
            raise

    def _coex_concordance(self,gene_a,gene_b,maxnan=10):
        '''
            This is a sanity method to ensure that the pcc calculated
            directly from the expr profiles matches the one stored in 
            the database
        '''
        expr_a = self.expr_profile(gene_a).values
        expr_b = self.expr_profile(gene_b).values
        mask = np.logical_and(np.isfinite(expr_a),np.isfinite(expr_b))
        if sum(mask) < maxnan:
            # too many nans to reliably calculate pcc 
            return np.nan
        r = pearsonr(expr_a[mask],expr_b[mask])[0]
        # fisher transform it
        z = np.arctanh(r-0.0000001)
        # standard normalize it
        z = (z - float(self._global('pcc_mean')))/float(self._global('pcc_std'))
        return z
 

    ''' ------------------------------------------------------------------------------------------
            Class Methods -- Factory Methods
    '''

    @classmethod
    def from_Expr(cls,expr):
        ''' 
            Create a COB instance from an camoco.Expr (Expression) instance.
            A COB inherits all the methods of a Expr instance and implements additional
            coexpression specific methods. This method accepts an already build Expr
            instance and then performs the additional computations needed to build a 
            full fledged COB instance.

            Parameters
            ----------
            expr : camoco.Expr
                
            Returns
            -------
            camoco.COB instance
            
        '''
        # The Expr object already exists, just get a handle on it
        self = expr
        self._calculate_coexpression()
        self._calculate_degree()
        return self

    @classmethod
    def from_DataFrame(cls,df,name,description,refgen,rawtype=None,**kwargs):
        '''
            The method will read the table in (as a pandas dataframe), 
            build the Expr object passing all keyword arguments in **kwargs
            to the classmethod Expr.from_DataFrame(...). See additional 
            **kwargs in COB.from_Expr(...)

            Parameters
            ----------
            df : pandas.DataFrame
                A Pandas dataframe containing the expression information.
                Assumes gene names are in the index while accessions (experiments)
                are stored in the columns.
            name : str
                Name of the dataset stored in camoco database
            description : str
                Short string describing the dataset
            refgen : camoco.RefGen
                A Camoco refgen object which describes the reference
                genome referred to by the genes in the dataset. This 
                is cross references during import so we can pull information
                about genes we are interested in during analysis. 
            rawtype : str (default: None)
                This is noted here to reinforce the impotance of the rawtype passed to 
                camoco.Expr.from_DataFrame. See docs there for more information.
            **kwargs : key value pairs
                additional parameters passed to subsequent methods. (see Expr.from_DataFrame)

        '''
        # Create a new Expr object from a data frame
        expr = super().from_DataFrame(
            df,name,description,refgen,rawtype,**kwargs
        )
        return cls.from_Expr(expr)

    @classmethod
    def from_table(cls,filename,name,description,refgen,rawtype=None,sep='\t',**kwargs):
        '''
            Build a COB Object from an FPKM or Micrarray CSV. This is a
            convenience method which handles reading in of tables.
            Files need to have gene names as the first column and 
            accession (i.e. experiment) names as the first row. All
            kwargs will be passed to COB.from_DataFrame(...). See
            docstring there for option descriptions.

            Parameters
            ----------
            filename : str (path)
                the path to the FPKM table in csv or tsv
            name : str
                Name of the dataset stored in camoco database
            description : str
                Short string describing the dataset
            refgen : camoco.RefGen
                A Camoco refgen object which describes the reference
                genome referred to by the genes in the dataset. This 
                is cross references during import so we can pull information
                about genes we are interested in during analysis. 
            rawtype : str (default: None)
                This is noted here to reinforce the impotance of the rawtype passed to 
                camoco.Expr.from_DataFrame. See docs there for more information.
            sep : str (default: \\t)
                Specifies the delimiter of the file referenced by the filename parameter.
            **kwargs : key value pairs
                additional parameters passed to subsequent methods.

            Returns
            -------
                a COB object
        '''
        return cls.from_DataFrame(pd.read_table(filename,sep=sep),name,description,refgen,rawtype=rawtype,**kwargs)


    '''
        Unimplemented ---------------------------------------------------------------------------------
    '''


    def next_neighbors(self,gene_list):
        ''' returns a list containing the strongest connected neighbors '''
        pass

    def neighborhood(self,gene_list):
        ''' Input: A gene List
            Output: a Dataframe containing gene ids which have at least one edge
                    with another gene in the input list. Also returns global degree '''
        pass

    def lcc(self,gene_list,min_distance=None):
        ''' returns an igraph of the largest connected component in graph '''
        pass

    def seed(self, gene_list, limit=65): 
        ''' Input: given a set of nodes, add on the next X strongest connected nodes ''' 
        pass

    def graph(self,gene_list,min_distance=None):
        ''' Input: a gene list
            Output: a iGraph object '''
        pass

    def coordinates(self,gene_list,layout=None):
        ''' returns the static layout, you can change the stored layout by passing 
            in a new layout object. If no layout has been stored or a gene does not have
            coordinates, returns (0,0) for each mystery gene'''
        pass

    def plot(self,gene_list,filename=None,width=3000,height=3000,layout=None,**kwargs):
        pass

    def compare_to_COB(self,COB_list,filename=None,gridsize=100,extent=[-10,10,-10,10]):
        ''' Compare the edge weights in this COB to another COB. Prints out edge weights to file'''
        for oCOB in COB_list:
            self.log("Comparing {} to {}",self.name,oCOB.name)
            filename = "{}_to_{}".format(self.name,oCOB.name)
            # Print out the edge comparisons for each common gene
            self.log("Printing out common gene edges")
            if not os.path.exists(filename+'.tsv'):
                with open(filename+'.tsv','w') as OUT:
                # merge the tables of edges on common genes
                    common_edges = self.coex.join(oCOB.coex, how='inner')

    def compare_to_dat(self,filename,sep="\t",score_cutoff=3):
        ''' Compare the number of genes with significant edges as well as degree with a DAT file '''
        pass

    def compare_degree(self,obj,score_cutoff=3):
        ''' Compares the degree of one COB to another '''
        pass
