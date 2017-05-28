#!/usr/bin/env python3

import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib.pylab as plt
matplotlib.style.use('ggplot')
#matplotlib.use('Agg')

import argparse
import sys
import os
import copy
import re
import sqlite3

import numpy as np
import scipy as sp
from scipy.stats import hypergeom
from scipy import tril_indices,triu_indices

from itertools import chain

import camoco as co
import pandas as pd
pd.set_option('display.width',300)
import glob as glob

from .Camoco import Camoco
from .Tools import log

class Overlap(Camoco):
    '''
        The Overlap class represents the statistical enrichment of co-expression
        among a set of loci. 

    '''
    def __init__(self,gwas):
        if not isinstance(gwas, str):
            gwas = gwas.name
        super().__init__(gwas,type='Overlap')
        self._global('gwas',gwas)
        self._create_tables()
        # Read the overlap results in from the database
        # It is unfortunate that we need to use the sqlite3 package here
        self.results = pd.read_sql(
            'SELECT * FROM overlap;',
            sqlite3.connect(self.db.filename)
        )

    '''
        Private Methods
    '''

    def _create_tables(self):
        cur = self.db.cursor()
        cur.execute('''
                CREATE TABLE IF NOT EXISTS "overlap" (
                  "COB" TEXT,
                  "FlankLimit" INTEGER,
                  "Method" TEXT,
                  "NumBootstraps" INTEGER,
                  "Ontology" TEXT,
                  "Term" TEXT,
                  "TermCollapsedLoci" INTEGER,
                  "TermLoci" INTEGER,
                  "TermPValue" REAL,
                  "WindowSize" INTEGER,
                  "bs_mean" REAL,
                  "bs_std" REAL,
                  "fdr" REAL,
                  "fitted" REAL,
                  "gene" TEXT,
                  "global" REAL,
                  "local" REAL,
                  "num_random" REAL,
                  "num_real" REAL,
                  "score" REAL,
                  "zscore" REAL,
                  UNIQUE (COB, Ontology, Term, gene, WindowSize, FlankLimit) ON CONFLICT REPLACE
                );
        ''')

    def _build_indices(self):
        cur = self.db.cursor()
        cur.execute('CREATE INDEX IF NOT EXISTS gene ON overlap(gene)')

    '''
        Methods
    '''
    def get_data(self, gene=None, cob=None, term=None,
        windowSize=None, flankLimit=None):
        '''
            Function to get data using any of the parameters it is normaly queried by
        '''
        # Base section of query
        query = "SELECT * FROM overlap WHERE"
        
        # Throw arguments into a dictionary for effective looping
        args = {
            'gene':gene, 'COB':cob, 'Term':term,
            'WindowSize':windowSize, 'FlankLimit':flankLimit}
        
        # For each argument, add a clase to the SQL query
        for k,v in args.items():
            if not(v is None):
                if isinstance(v,(set,list)):
                    ls = "{}".format("','".join([str(x) for x in v]))
                else:
                    ls = v
                query += " {} IN ('{}') AND".format(k,ls)
        
        # Peel off unneeded things at the end of the query, depending args
        if query.endswith(' AND'):
            query = query.rstrip(' AND')
        if query.endswith(' WHERE'):
            query = query.rstrip(' WHERE')
        query += ';'
        
        # Run the query
        cur = self.db.cursor()
        print('Executing: {}'.format(query))
        cur.execute(query)
        
        # Throw the results into a DataFrame for conviniece
        return pd.read_sql(query,sqlite3.connect(self.db.filename)).set_index('gene')


    def generate_output_name(self):
        # Handle the different output schemes
        if self.args.out is None:
            self.args.out = '{}_{}_{}_{}_{}_{}.tsv'.format(
                self.cob.name,
                self.ont.name,
                self.args.candidate_window_size,
                self.args.candidate_flank_limit,
                self.args.method,
                ':'.join(self.args.terms)
            )
        if not self.args.out.endswith('.overlap.tsv'):
            self.args.out = self.args.out + '.overlap.tsv'
        if os.path.dirname(self.args.out) != '':
            os.makedirs(os.path.dirname(self.args.out),exist_ok=True)
        if os.path.exists(self.args.out) and self.args.force != True:
            print(
                "Output for {} exists! Skipping!".format(
                    self.args.out
                ),file=sys.stderr
            )
            return

    def snp2gene(self,term):
        if 'effective' in self.args.snp2gene:
            # Map to effective
            return term.effective_loci(
                window_size=self.args.candidate_window_size
            )
        elif 'strongest' in self.args.snp2gene:
            return term.strongest_loci(
                window_size=self.args.candidate_window_size,
                attr=self.args.strongest_attr,
                lowest=self.args.strongest_higher
            )

    def generate_bootstraps(self,loci,overlap):
        '''
            bootstrapping procedure. Our target here is to provide enough bootstraps
            to identify loci that are significant at n==1000 bootstraps. The auto 
            procedure will continue untill we meet n==1000 OR we find 50 bootstraps
            that are have higher score such that we will never be significant at 1000
            boostraps (0.05 * 1000 = 50).
        '''
        target_score = overlap.score.mean()
        max_bs = 1000
        num_bs = 0
        bs = []
        if self.args.num_bootstraps == 'auto':
            # Create a bullshit generator... err bootstraps
            bs_generator = (co.Overlap.overlap(self,loci,bootstrap=True,iter_name=x) \
                for x in range(max_bs)
            )
            while num_bs <= 50 and len(bs) < 1000: 
                # Add 50 bootstraps to current bootstraps
                bs = [next(bs_generator) for x in range(50)] + bs
                # Find the number of bs more extreme than the empirical score
                num_bs = sum([df.score.mean() >= target_score for df in bs])
                self.cob.log(
                    "Iteration: {} -- current pval: {} {}% complete",
                    len(bs),num_bs/len(bs),num_bs/50
                )
        else:
            # Be a lil broke back noodle and explicitly bootstrap
            bs = [co.Overlap.overlap(self, loci,bootstrap=True,iter_name=x) \
                for x in range(int(self.args.num_bootstraps))
            ]
        return pd.concat(bs)


    def overlap(self,loci,bootstrap=False,iter_name=None):
        '''
            Calculate Network-Term Overlap based on the method in CLI args

            This method must return a DataFrame 
        '''
        # generate emirical and bootstrap values 
        if self.args.method == 'density':
            return self.cob.trans_locus_density(
                loci,
                flank_limit=self.args.candidate_flank_limit,
                by_gene = True,
                bootstrap = bootstrap,
                iter_name = iter_name
            ) 
        elif self.args.method == 'locality':
            return self.cob.trans_locus_locality(
                loci,
                flank_limit=self.args.candidate_flank_limit,
                by_gene=True,
                bootstrap=bootstrap,
                iter_name=iter_name,
                include_regression=True,
            ).rename(columns={'resid':'score'})


    def num_below_fdr(self,fdr_cutoff=0.3):
        return pd.pivot_table(
            self.results,
            index=['Method','COB','WindowSize','FlankLimit'],
            columns='Term',
            values='fdr',
            aggfunc=lambda x: sum(x<=fdr_cutoff)
        )

    def SNP2Gene_breakdown(self,COB=None):
        '''
        Provides a breakdown of SNP to gene mapping parameters for each term in the Overlap.
        Includes the number of initial Loci, the number of collapsed Loci (within a window)
        and the number of candidate genes (within a window and up to a flank limit)

        Parameters
        ----------
        COB : str (default: 'average')
            If specfified, the results will be composed only of SNP to gene
            mappings from a single COB network. If 'average' is specified,
            the results will be the SET of genes across all COB networks.
        '''
        # Get some help
        def bp_to_kb(bp):
            return "{}KB".format(int(bp/1000))
        def get_level(df,level):
            ''' Returns the level values by name '''
            level_index = df.columns.names.index(level)
            return df.columns.levels[level_index]
        # Prepare the data frame results
        if COB == None:
            results = self.results
        else:
            results = self.results.query('COB=="{}"'.format(COB))
        # Total for the Ionome
        ont = co.GWAS(self.results.Ontology.unique()[0])
        ref = co.COB(self.results.COB.unique()[0])._parent_refgen
        # Make an aggregate term
        total = co.Term('total',loci=set(chain(* [x.loci for x in ont.terms()])))
        # Calculate number of SNPs
        snps = pd.DataFrame(pd.pivot_table(results,index="Term",values='TermLoci'))
        snps.columns = pd.MultiIndex.from_product([['GWAS SNPs'],['-'],['-']],names=['Name','WindowSize','FlankLimit'])
        snps.ix['Total'] = len(total.loci)
        # Calculate number of Candidate Loci
        loci = pd.pivot_table(results,index="Term",columns=['WindowSize'],values='TermCollapsedLoci')
        for window_size in loci.columns:
            loci.ix['Total',window_size] = len(total.effective_loci(window_size))
        loci.columns = pd.MultiIndex.from_product([['Collapsed Loci'],list(map(bp_to_kb,loci.columns)),['-']],names=['Name','WindowSize','FlankLimit'])
        # Calculate number of Candidate Genes
        genes = pd.pivot_table(results,index='Term',columns=['WindowSize','FlankLimit'],values='gene',aggfunc=lambda x: len(set(x)))
        for window_size in get_level(genes,'WindowSize'):
            for flank_limit in get_level(genes,'FlankLimit'):
                genes.ix['Total',(window_size,flank_limit)] = len(ref.candidate_genes(total.effective_loci(window_size=window_size),flank_limit=flank_limit))
        genes.columns = pd.MultiIndex.from_product(
            [['Candidate Genes'],
                list(map(bp_to_kb,get_level(genes,"WindowSize"))),
                get_level(genes,'FlankLimit')
            ],
            names=['Name','WindowSize','FlankLimit']
        )
        results = snps.join(loci).join(genes)
        #ionome_eff_loci = [len()]
        return results.astype(int)
     
    def plot_pval_heatmap(self,filename=None,pval_cutoff=0.05,
                          collapse_snp2gene=False,figsize=(15,10),
                          skip_terms=None):
        '''
            Generate a heatmap based on TermPVal

            Parameters
            ----------
            filename : str
                output file name
            pval_cutoff : float (default:0.05)
                The term p-value cutoff for shading the heatmap
            collapse_snp2gene : bool (default: False)
                If true, candidates will be collapsed around
                snp to gene mapping parameters
            figsize : tuple(int,int) (default:(15,10))
                Control the size of the figure
            skip_terms : iter (default:None)
                If provided, terms in the iterable will
                not be plotted
        '''
        methods = sorted(self.results.Method.unique())
        cobs = sorted(self.results.COB.unique())
        fig,axes = plt.subplots(len(cobs),len(methods))
        fig.set_size_inches(figsize)
        def get_axis(i,j,axes):
            if len(axes.shape) == 1:
                return axes[i]
            else:
                return axes[i,j]
        for i,cob in enumerate(cobs):
            for j,method in enumerate(methods):
                axis = get_axis(i,j,axes)
                data = pd.pivot_table(
                    self.results.ix[(self.results.COB==cob) & (self.results.Method==method)],
                    index=['WindowSize','FlankLimit'],
                    columns=['Term'],
                    values='TermPValue',
                    aggfunc=np.mean
                )
                if collapse_snp2gene:
                    data = pd.DataFrame(data.apply(min)).T
                if skip_terms:
                    for term in skip_terms:
                        del data[term]
                #data[data > pval_cutoff] = np.nan
                #data[data < pval_cutoff] = 0
                axis.set_frame_on(False)
                cmap = plt.cm.Greens_r
                cmap.set_bad('white',1.0)
                im = axis.pcolormesh(
                    np.ma.masked_invalid(data),
                    cmap=cmap,
                    #linewidth=0.1,
                    edgecolor='lightgrey',
                    vmin=0,
                    vmax=0.05
                )
                # Make the layout more natural
                if j == len(methods) - 1 :
                    axis.set_ylabel(cob,fontsize=10)
                    axis.set_yticklabels(data.index.values,rotation='45')
                    axis.yaxis.set_label_position("right")
                if j == 0:
                    axis.set_yticklabels(
                        axis.get_yticklabels(),size=7
                    )
                    axis.set(
                        yticks=np.arange(len(data.index))+0.5
                    )
                else:
                    axis.set_yticks([])
                    axis.set_yticklabels([])
                axis.yaxis.set_ticks_position('left')
                axis.invert_yaxis()
                if i == 0:
                    axis.set_title(method,y=1.1,size=15)
                    axis.xaxis.tick_top()
                    axis.set_xticks(np.arange(len(data.columns))+0.5)
                    axis.set_xticklabels(
                        [re.sub('\d','',x) for x in data.columns.values], 
                        rotation='45',
                        ha='left',
                        size=7
                    )
                else:
                    axis.set_xticks([])
                    axis.set_xticklabels([])
                axis.set_yticklabels(data.index.values)
        # Create a new axis to append the color bar to
        #fig.subplots_adjust(right=0.8)
        #cax = fig.add_axes([.95,0.5,0.01,0.4])
        #divider = make_axes_locatable(get_axes(0,0,))
        #cax = divider.append_axes("bottom", size="5%", pad=0.05)
        #cbar = fig.colorbar(
        #    im, cax=cax, orientation='vertical', ticks=[0,0.025,0.05]
        #)
        #cbar.set_ticklabels(['0','0.025','≥0.05'])
        #plt.tight_layout(pad=0.4,w_pad=0.5, h_pad=1.0)
        if filename is not None:
            plt.savefig(filename,dpi=300)
            plt.close()
        return (fig,axes)

    def high_priority_candidates(self,fdr_cutoff=0.3,min_snp2gene_obs=2):
        '''
            Return the number of candidate genes, seen at multiple SNP2Gene
            mappings and in multiple networks (both)
        '''
        df = self.results[np.isfinite(self.results.fdr)]
        df = df[df.fdr <= fdr_cutoff]
        # Filter out genes that do not occur at 2+ SNP to gene mappings, 
        # they will never be included
        #df['fdr'] = fdr_cutoff
        df = df.groupby(
            ['COB','Term','Method','gene']
        ).filter(lambda df: len(df)>=min_snp2gene_obs).copy()
        # Find genes that were significant in either density or locality
        either = df.copy()
        either.loc[:,'Method'] = 'either'
        either.loc[:,'COB'] = 'Any'
        either.drop_duplicates(inplace=True)
        # Find genes that are significant in both density and locality 
        # in the same element and the same network
        both_same_net = df.groupby(
            ['COB','Term','gene']
        ).filter(lambda df: len(df.Method.unique()) == 2).copy()
        both_same_net.loc[:,'Method'] = 'both_same_net'
        # Find genes that are significant in both density and locality 
        # in the same element but any network
        both_any_net = df.groupby(
            ['Term','gene']
        ).filter(lambda df: len(df.Method.unique()) == 2).copy()
        both_any_net.loc[:,'Method'] = 'both_any_net'
        both_any_net['COB'] = 'Any'
        # Concatenate
        return pd.concat([df,both_same_net,both_any_net,either])

    def adjacency(self, min_snp2gene_obs=2,fdr_cutoff=0.3):
        '''
            Return a matrix showing the number of shared HPO genes by Term.
            The diagonal of the matrix is the number of genes discoverd by that 
            term. The upper diagonal shows the overlap between the row and column
            and the lower diagonal shows the hypergeomitric pval for the overlap
            between the two terms. The universe used is the number of unique genes
            in the overlap results.
        '''
        df = self.high_priority_candidates(fdr_cutoff=fdr_cutoff,min_snp2gene_obs=min_snp2gene_obs)
        # 
        x={df[0]:set(df[1].gene) for df in df.groupby('Term')}                     
        adj = []                                                                        
        #num_universe = len(set(chain(*x.values())))
        num_universe = len(self.results.gene.unique())
        for a in x.keys():                                                              
            for b in x.keys():                                                          
                num_common = len(set(x[a]).intersection(x[b]))
                if a != b:
                    pval = hypergeom.sf(num_common-1,num_universe,len(x[a]),len(x[b]))
                else:
                    # This will make the diagonal of the matrix be the number HPO genes
                    # for the element
                    pval = len(x[a])
                adj.append((a,b,num_common,pval)) 
        adj = pd.DataFrame(adj)                                                         
        overlap = pd.pivot_table(adj,index=0,columns=1,values=2)
        # Mask out the lower diagonal on the overalp matrix
        overlap.values[tril_indices(len(overlap))] = 0
        pvals = pd.pivot_table(adj,index=0,columns=1,values=3)
        # Mask out the upper tringular on the pvals matrix
        pvals.values[triu_indices(len(pvals),1)] = 0
        return overlap+pvals

    def num_candidate_genes(self,fdr_cutoff=0.3,min_snp2gene_obs=2):
        candidates = self.high_priority_candidates(fdr_cutoff=fdr_cutoff,min_snp2gene_obs=min_snp2gene_obs)
        candidates['fdr'] = fdr_cutoff
        # Pivot and aggregate
        return pd.pivot_table(
            candidates,
            columns=['fdr','Method','COB'],
            index=['Term'],
            values='gene',
            aggfunc=lambda x: len(set(x))
        ).fillna(0).astype(int)

    
    ''' ----------------------------------------------------------------------------------
        Class Methods
    '''
    @classmethod
    def create(cls, gwas, description=None):
        if not isinstance(gwas, str):
            gwas = gwas.name
        self = super().create(gwas, description, type='Overlap')
        return self
 
    @classmethod
    def from_csv(cls,dir='./',sep='\t'):
        # Read in Gene specific data
        results = pd.concat(
            [pd.read_table(x,sep=sep) \
                for x in glob.glob(dir+"/*.overlap.tsv") ]
        )
        gwas = results.Ontology.unique()[0]
        self = cls.create(gwas)
        self.results = results
        # Add the results to the sqlite table
        #self._bcolz('overlap',df=self.results)
        self.results.to_sql('overlap',sqlite3.connect(self.db.filename),if_exists='replace') 
        return self

    @classmethod
    def from_CLI(cls,args):
        '''
            Implements an interface for the CLI to perform overlap
            Analysis
        '''
        self = lambda:None 
        # Build base camoco objects
        self.args = args
        self.cob = co.COB(args.cob)
        self.ont = co.GWAS(args.gwas)
        self.generate_output_name = co.Overlap.generate_output_name(self)

        # Generate a terms iterable
        if 'all' in self.args.terms:
            terms = self.ont.iter_terms()
        else:
            terms = [self.ont[term] for term in self.args.terms]
        all_results = list()

        results = []
        # Iterate through terms and calculate
        for term in terms:
            if term.id in self.args.skip_terms:
                self.cob.log('Skipping {} since it was in --skip-terms',term.id)
            self.cob.log(
                "Calculating Overlap for {} of {} in {} with window:{} and flank:{}",
                term.id,
                self.ont.name,
                self.cob.name,
                self.args.candidate_window_size,
                self.args.candidate_flank_limit
            )
            if args.dry_run:
                continue
            # Generate SNP2Gene Loci
            loci = co.Overlap.snp2gene(self,term)
            if len(loci) < 2:
                self.cob.log('Not enough genes to perform overlap')
                continue
            self.overlap = co.Overlap.overlap(self, loci)
            overlap = co.Overlap.overlap(self, loci)
            bootstraps = co.Overlap.generate_bootstraps(self,loci,overlap)
            bs_mean = bootstraps.groupby('iter').score.apply(np.mean).mean()
            bs_std  = bootstraps.groupby('iter').score.apply(np.std).mean()
            # Calculate z scores for density
            if bs_std != 0:
                overlap['zscore'] = (overlap.score-bs_mean)/bs_std
                bootstraps['zscore'] = (bootstraps.score-bs_mean)/bs_std
            else:
                # If there is no variation, make all Z-scores 0
                overlap['zscore'] = bootstraps['zscore'] = 0
            # Calculate FDR
            overlap['fdr'] = np.nan
            max_zscore = int(overlap.zscore.max()) + 1
            for zscore in np.arange(0, max_zscore,0.25):
                num_random = bootstraps\
                        .groupby('iter')\
                        .apply(lambda df: sum(df.zscore >= zscore))\
                        .mean()
                num_real = sum(overlap.zscore >= zscore)
                # Calculate FDR
                if num_real != 0 and num_random != 0:
                    fdr = num_random / num_real
                elif num_real != 0 and num_random == 0:
                    fdr = 0
                else:
                    fdr = 1
                overlap.loc[overlap.zscore >= zscore,'fdr'] = fdr
                overlap.loc[overlap.zscore >= zscore,'num_real'] = num_real
                overlap.loc[overlap.zscore >= zscore,'num_random'] = num_random
                overlap.loc[overlap.zscore >= zscore,'bs_mean'] = bs_mean
                overlap.loc[overlap.zscore >= zscore,'bs_std'] = bs_std
                overlap.sort_values(by=['zscore'],ascending=False,inplace=True)
            overlap_pval = (
                (sum(bootstraps.groupby('iter').apply(lambda x: x.score.mean()) >= overlap.score.mean()))\
                / len(bootstraps.iter.unique())
            )
            # This gets collated into all_results below
            overlap['COB'] = self.cob.name
            overlap['Ontology'] = self.ont.name
            overlap['Term'] = term.id
            overlap['WindowSize'] = self.args.candidate_window_size
            overlap['FlankLimit'] = self.args.candidate_flank_limit
            overlap['TermLoci'] = len(term.loci)
            overlap['TermCollapsedLoci'] = len(loci)
            overlap['TermPValue'] = overlap_pval
            overlap['NumBootstraps'] = len(bootstraps.iter.unique())
            overlap['Method'] = self.args.method
            results.append(overlap.reset_index())
        if not args.dry_run:
            self.results = pd.concat(results)
            self.results.to_csv(self.args.out,sep='\t',index=None)
            overlap_object = cls.create(self.ont)
            overlap_object.results = results
            self.results.to_sql('overlap',sqlite3.connect(overlap_object.db.filename),if_exists='append',index=False)
            


