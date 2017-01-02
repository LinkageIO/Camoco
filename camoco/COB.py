#!/usr/bin/python3

import camoco.PCCUP as PCCUP

from .Camoco import Camoco
from .RefGen import RefGen
from .Locus import Locus,Gene
from .Expr import Expr
from .Tools import memoize,available_datasets
from .Term import Term
from .Ontology import Ontology

from math import isinf
from numpy import matrix, arcsinh, tanh
from collections import defaultdict, Counter
from itertools import chain
from subprocess import Popen, PIPE
from scipy.spatial.distance import squareform
from scipy.misc import comb
from scipy.stats import norm
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
from statsmodels.sandbox.regression.predstd import wls_prediction_std

import matplotlib.pyplot as plt
import statsmodels.api as sm
import networkx as nx
import pandas as pd
import numpy as np
import itertools
from scipy.misc import comb 

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

import statsmodels.api as sm
import sys
import json
import gc

from scipy.stats import pearsonr

class COB(Expr):
    def __init__(self, name):
        super().__init__(name=name)
        self.log('Loading coex table')
        self.coex = self._bcolz('coex')
        if self.coex is None:
            self.log("{} is empty", name)
        self.log('Loading Global Degree')
        self.degree = self._bcolz('degree')
        if self.degree is None:
            self.log("{} is empty", name)
        if not available_datasets('Ontology','{}MCL'.format(name))\
            and self.coex is not None:
            self._calculate_clusters()
        self.log('Loading Clusters')
        self.clusters = self._bcolz('clusters')
        if self.clusters is None:
            self.log('Clusters not loaded for: {} ()', name)
            self.MCL = None
        else:
            self.MCL = Ontology('{}MCL'.format(self.name))


    def __repr__(self):
        return '<COB: {}>'.format(self.name)

    def __str__(self):
        return self.__repr__()

    def summary(self,file=sys.stdout):
        print( '''
            COB Dataset: {}
                Desc: {}
                RawType: {}
                TransformationLog: {}
                Num Genes: {:,}({:.2g}%)
                Num Accessions: {}
                Num Edges: {:,}

            Raw
            ------------------
            Num Raw Genes: {:,}
            Num Raw Accessions: {}

            QC Parameters
            ------------------
            min expr level: {} 
                - expression below this is set to NaN
            max gene missing data: {} 
                - genes missing more than this percent are removed
            max accession missing data: {}
                - Accession missing more than this percent are removed
            min single sample expr: {} 
                - genes must have at least this amount of expression in 
                  on accession

            Clusters
            ------------------
            Num clusters (size >= 10): {}


        '''.format(
            # Dataset
            self.name,
            self.description,
            self.rawtype,
            self._transformation_log(),
            self.num_genes(),
            (self.num_genes()/self.num_genes(raw=True))*100,
            self.num_accessions(),
            len(self.coex),
            # Raw
            len(self.expr(raw=True)),
            len(self.expr(raw=True).columns),
            # QC
            self._global('qc_min_expr'),
            self._global('qc_max_gene_missing_data'),
            self._global('qc_min_single_sample_expr'),
            self._global('qc_max_accession_missing_data'),
            # Clusters
            sum(self.clusters.groupby('cluster').apply(len) >= 10)
        ), file=file)

    def qc_gene(self):
        qc_gene = self._bcolz('qc_gene')
        qc_gene['chrom'] = [self.refgen[x].chrom for x in qc_gene.index]
        return qc_gene.groupby('chrom').aggregate(sum, axis=0)

    @property
    @memoize
    def edge_FDR(self):
        # get the percent of significant edges
        num_sig = np.sum(self.coex['significant'])/len(self.coex)
        # calulate the number expected
        num_exp = 1-norm.cdf(int(self._global('significance_threshold')))
        # FDR is the percentage expected over the percentage found
        return num_exp/num_sig

    def set_sig_edge_zscore(self,zscore):
        self.coex.significant = self.coex.score >= zscore

    def neighbors(self, gene, sig_only=True, names_as_index=True, names_as_cols=False):
        '''
            Returns a DataFrame containing the neighbors for gene.

            Parameters
            ----------
            gene : co.Locus
                The gene for which to extract neighbors
            sig_only : bool
                A flag to include only significant interactions.
            names_as_index : bool (default: True)
                Include gene names as the index.
            names_as_cols : bool (default: False)
                Include gene names as two columns named 'gene_a' and 'gene_b'.

            Returns
            -------
            A DataFrame containing edges
        '''
        # Find the neighbors
        gene_id = self._expr_index[gene.id]
        neighbor_indices = PCCUP.coex_neighbors(gene_id, self.num_genes())
        neighbor_indices.sort()
        edges = self.coex.iloc[neighbor_indices]
        del neighbor_indices
        if sig_only:
            edges = edges[edges.significant == 1]
        if len(edges) == 0:
            return pd.DataFrame(columns=['score','distance','significant'])
        
        # Find the indexes if necessary
        if names_as_index or names_as_cols:
            names = self._expr.index.values
            ids = edges.index.values
            ids = PCCUP.coex_expr_index(ids, self.num_genes())
            edges.insert(0,'gene_a', names[ids[:,0]])
            edges.insert(1,'gene_b', names[ids[:,1]])
            del ids; del names;
        if names_as_index:
            edges = edges.set_index(['gene_a','gene_b'])
        
        return edges
        
    def coexpression(self, gene_a, gene_b):
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
        if gene_a.id == gene_b.id:
            # We don't cache these results
            score = self._coex_concordance(gene_a,gene_b)
            significant = 0
            distance = 0
            return pd.Series(
                [score,significant,distance],
                name=(gene_a.id,gene_b.id),
                index = ['score','significant','distance']
            )
        return self.subnetwork([gene_a,gene_b],sig_only=False).iloc[0]

    def subnetwork(self, gene_list=None, sig_only=True, min_distance=None,
        filter_missing_gene_ids=True, trans_locus_only=False,
        names_as_index=True, names_as_cols=False):
        '''
            Extract a subnetwork of edges exclusively between genes
            within the gene_list. Also includes various options for
            what information to report, see Parameters.

            Parameters
            ---------
            gene_list : iter of Loci
                The genes from which to extract a subnetwork.
                If gene_list is None, the function will assume
                gene_list is all genes in COB object (self).
            sig_only : bool
                A flag to include only significant interactions.
            min_distance : bool (default: None)
                If not None, only include interactions that are
                between genes that are a `min_distance` away from
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
        '''
        num_genes = self.num_genes()
        if gene_list is None:
            # Return the entire DataFrame
            df = self.coex.copy()
        else:
            # Extract the ids for each Gene
            gene_list = set(sorted(gene_list))
            ids = np.array([self._expr_index[x.id] for x in gene_list])
            if filter_missing_gene_ids:
                # filter out the Nones
                ids = np.array([x for x in ids if x is not None])
            if len(ids) == 0:
                df = pd.DataFrame(columns=['score','significant','distance'])
            else:
                # Grab the coexpression indices for the genes
                ids = PCCUP.coex_index(ids, num_genes)
                ids.sort()
                df = self.coex.iloc[ids]
                del ids
        if min_distance is not None:
            df = df[df.distance >= min_distance]
        if names_as_index or names_as_cols or trans_locus_only:
            names = self._expr.index.values
            ids = df.index.values
            if len(ids) > 0:
                ids = PCCUP.coex_expr_index(ids, num_genes)
                df.insert(0,'gene_a', names[ids[:,0]])
                df.insert(1,'gene_b', names[ids[:,1]])
                del ids; del names;
            else:
                df.insert(0,'gene_a',[])
                df.insert(0,'gene_b',[])
        if names_as_index:
            df = df.set_index(['gene_a','gene_b'])
        if trans_locus_only:
            try:
                parents = {x.id:x.attr['parent_locus'] for x in gene_list}
            except KeyError as e:
                raise KeyError(
                    "Each locus must have 'parent_locus'"
                    " attr set to calculate trans only"
                )
            df['trans'] = [
                parents[gene_a] != parents[gene_b] for gene_a,gene_b in \
                zip(df.index.get_level_values(0),df.index.get_level_values(1))
            ]
        if sig_only:
            df = df.ix[df.significant]
        return df

    def cluster_coefficient(self, locus_list, flank_limit,
        trans_locus=True, bootstrap=False, by_gene=True, iter_name=None):
        '''
            Calculates the clustering coefficient for genes which span loci.

            Parameters
            ----------
            locus_list : iter of Loci
                an iterable of loci
            flank_limit : int
                The number of flanking genes passed to be pulled out
                for each locus (passed onto the refgen.candidate_genes method)
            return_mean : bool (default: True)
                If false, raw edges will be returned
            bootstrap : bool (default: False)
                If true, candidate genes will be bootstrapped from the COB
                reference genome
            by_gene : bool (default: False)
                Return a per-gene breakdown of density within the subnetwork.
            iter_name : str (default: None)
                Optional string which will be added as a column. Useful for
                keeping track of bootstraps in an aggregated data frame.

            Returns
            -------
            Clustering coefficient of interactions if return_mean is True
            otherwise a dataframe of trans edges

        '''
        raise NotImplementedError()


    def trans_locus_density(self, locus_list,flank_limit,
        return_mean=True, bootstrap=False, by_gene=False,
        iter_name=None):
        '''
            Calculates the density of edges which span loci. Must take in a locus
            list so we can exlude cis-locus interactions.

            Parameters
            ----------
            locus_list : iter of Loci
                an iterable of loci
            flank_limit : int
                The number of flanking genes passed to be pulled out
                for each locus (passed onto the refgen.candidate_genes method)
            return_mean : bool (default: True)
                If false, raw edges will be returned
            bootstrap : bool (default: False)
                If true, candidate genes will be bootstrapped from the COB
                reference genome
            by_gene : bool (default: False)
                Return a per-gene breakdown of density within the subnetwork.
            iter_name : str (default: None)
                Optional string which will be added as a column. Useful for
                keeping track of bootstraps in an aggregated data frame.

            Returns
            -------
            Z-score of interactions if return_mean is True
            otherwise a dataframe of trans edges

        '''
        # convert to list of loci to lists of genes
        if not bootstrap:
            genes_list = self.refgen.candidate_genes(
                locus_list, flank_limit=flank_limit, chain=True,
                include_parent_locus=True
            )
        else:
            genes_list = self.refgen.bootstrap_candidate_genes(
                locus_list, flank_limit=flank_limit, chain=True,
                include_parent_locus=True
            )
        # Extract the edges for the full set of genes
        edges = self.subnetwork(
            genes_list,
            min_distance=0,
            sig_only=False,
            trans_locus_only=True,
            names_as_index=True
        )
        if by_gene == True:
            # Filter out trans edges
            gene_split = pd.DataFrame.from_records(
                chain(
                    *[((gene_a,score),(gene_b,score)) \
                     for gene_a,gene_b,score,*junk \
                     in edges[edges.trans==True].reset_index().values]
                ),columns=['gene','score']
            )
            gene_split = gene_split.groupby('gene').agg(np.mean)
            if iter_name is not None:
                gene_split['iter'] = iter_name
            gene_split.index.name = 'gene'
            return gene_split
        else:
            if return_mean:
                scores = edges.loc[edges['trans']==True, 'score']
                return np.nanmean(scores)/(1/np.sqrt(len(scores)))
            else:
                return edges.loc[edges['trans']==True,]

    def trans_locus_locality(self, locus_list, flank_limit,
        bootstrap=False, by_gene=False, iter_name=None,
        include_regression=False):
        '''
            Computes a table comparing local degree to global degree
            of genes COMPUTED from a set of loci.
            NOTE: interactions from genes originating from the same
            locus are not counted for global or local degree.

            Parameters
            ----------
            locus_list : iterable of camoco.Loci
                A list or equivalent of loci
            flank_limit : int
                The number of flanking genes passed to be pulled out
                for each locus (passed onto the refgen.candidate_genes method)
            bootstrap : bool (default: False)
                If true, candidate genes will be bootstrapped from the COB
                reference genome
            iter_name : object (default: none)
                This will be added as a column. Useful for
                generating bootstraps of locality and keeping
                track of which one a row came from after catting
                multiple bootstraps together.
            by_gene : bool (default: False)
                Return a per-gene breakdown of density within the subnetwork.
            include_regression : bool (default: False)
                Include the OLS regression residuals and fitted values
                on local ~ global.

            Returns
            -------
            A pandas DataFrame with local, global and residual columns
            based on linear regression of local on global degree.
        '''
        # convert to list of loci to lists of genes
        if not bootstrap:
            genes_list = self.refgen.candidate_genes(
                locus_list, flank_limit=flank_limit, chain=True,
                include_parent_locus=True
            )
        else:
            genes_list = self.refgen.bootstrap_candidate_genes(
                locus_list, flank_limit=flank_limit, chain=True,
                include_parent_locus=True
            )
        #self.log("Found {} candidate genes", len(genes_list))
        # Get global and local degree for candidates
        gdegree = self.global_degree(genes_list, trans_locus_only=True)
        ldegree = self.local_degree(genes_list, trans_locus_only=True)
        # Merge the columns
        degree = ldegree.merge(gdegree,left_index=True,right_index=True)
        degree.columns = ['local', 'global']
        degree = degree.sort_values(by='global')
        degree.index.name = 'gene'
        if include_regression:
            # Add the regression lines
            ols = sm.OLS(degree['local'], degree['global']).fit()
            degree['resid'] = ols.resid
            degree['fitted'] = ols.fittedvalues
            degree = degree.sort_values(by='resid',ascending=False)
        if iter_name is not None:
            degree['iter'] = iter_name
        return degree

    def density(self, gene_list, min_distance=None, by_gene=False):
        '''
            Calculates the denisty of the non-thresholded network edges
            amongst genes within gene_list. Includes parameters to perform
            measurements for genes within a certain distance of each other.
            This corrects for cis regulatory elements increasing noise
            in coexpression network.

            Parameters
            ----------
            gene_list : iter of Loci
                List of genes from which to calculate density.
            min_distance : int (default: None)
                Ignore edges between genes less than min_distance
                in density calculation.
            by_gene : bool (default: False)
                Return a per-gene breakdown of density within the subnetwork.

            Returns
            -------
            A network density OR density on a gene-wise basis
        '''
        # filter for only genes within network
        edges = self.subnetwork(gene_list,
            min_distance=min_distance, sig_only=False
        )

        if by_gene == True:
            x = pd.DataFrame.from_records(
                chain(
                    *[((gene_a,score),(gene_b,score)) \
                     for gene_a,gene_b,score,sig,dis \
                     in edges.reset_index().values]
                ),columns=['gene','score']
            )
            return x.groupby('gene').agg(np.mean)
        else:
            if len(edges) == 0:
                return np.nan
            if len(edges) == 1:
                return edges.score[0]
            return np.nanmean(edges.score)/(1/np.sqrt(len(edges)))

    def to_dat(self, gene_list=None, filename=None, sig_only=True, min_distance=0):
        '''
            Outputs a .DAT file (see Sleipnir library)
        '''
        if filename is None:
            filename = self.name + '.dat'
        with open(filename, 'w') as OUT:
            # Get the score table
            self.log('Pulling the scores for the .dat')
            score = self.subnetwork(gene_list, sig_only=sig_only, 
            min_distance=min_distance, names_as_index=False,
            names_as_cols=False)
            del self.coex
            
            # Drop unecessary columns
            score.drop(['distance','significant'], axis=1, inplace=True)
            
            # Find the ids from those
            self.log("Finding the IDs")
            names = self._expr.index.values
            ids = PCCUP.coex_expr_index(score.index.values, self.num_genes())
            score.insert(0,'gene_a', names[ids[:,0]])
            score.insert(1,'gene_b', names[ids[:,1]])
            del ids; del names;
            
            # Print it out!
            self.log('Writing the .dat')
            score.to_csv(
                OUT, columns=['gene_a','gene_b','score'],index=False, sep='\t')
            del score
            self.coex = self._bcolz('coex')
            self.log('Done')

    def to_graphml(self,file, gene_list=None, sig_only=True, min_distance=0):
        # Get the edge indexes
        self.log('Getting the network.')
        edges = self.subnetwork(gene_list=gene_list, sig_only=sig_only,
        min_distance=min_distance, names_as_index=False,
        names_as_cols=False).index.values
        del self.coex
        
        # Find the ids from those
        names = self._expr.index.values
        edges = PCCUP.coex_expr_index(edges, self.num_genes())
        df = pd.DataFrame(index=np.arange(edges.shape[0]))
        df['gene_a'] = names[edges[:,0]]
        df['gene_b'] = names[edges[:,1]]
        del edges; del names;

        # Build the NetworkX network
        self.log('Building the graph.')
        net = nx.from_pandas_dataframe(df,'gene_a','gene_b')
        del df;

        # Print the file
        self.log('Writing the file.')
        nx.write_graphml(net,file)
        del net;
        self.coex = self._bcolz('coex')
        return

    def to_treeview(self, filename, cluster_method='mcl', gene_normalize=True):
        dm = self.expr(gene_normalize=gene_normalize)
        if cluster_method == 'leaf':
            order = self._bcolz('leaves').sort('index').index.values
        elif cluster_method == 'mcl':
            order = self._bcolz('clusters').loc[dm.index].\
                    fillna(np.inf).sort('cluster').index.values
        else:
            order = dm.index
        dm = dm.loc[order, :]

        aliases_raw = self.refgen.db.cursor().execute('SELECT alias,id FROM aliases').fetchall()
        aliases = dict()
        for alias,id in aliases_raw:
            if id in aliases:
                aliases[id] += alias + ' '
            else:
                aliases[id] = alias + ' '
        als = pd.Series(aliases)
        dm.insert(0,'Aliases',als)
        dm.to_csv(filename)


    def to_json(self, gene_list=None, filename=None, sig_only=True, min_distance=None):
        '''
            Produce a JSON network object that can be loaded in cytoscape.js
            or Cytoscape v3+.
        '''
        net = {
            'nodes' : [],
            'edges' : []
        }
        parents = defaultdict(list)
        # generate the subnetwork for the genes
        for gene in gene_list:
            net['nodes'].append(
                    {'data':{
                        'id':str(gene.id),
                        'parent':str(gene.attr['parent_locus']),
                        'classes':'gene'
                    }}
            )
            parents[str(gene.attr['parent_locus'])].append(gene.id)
        # Add parents first
        for parent,children in parents.items():
            net['nodes'].insert(0,
                {'data':{
                    'id':parent,
                    'parent':None,
                    'classes':'snp'
                }}
            )
            for child in children:
                net['edges'].append(
                    {'data':{
                        'source': child,
                        'target' : parent,
                        'score' : 50,
                        'distance' : 0
                    }}
                )
        # Get the edge indexes
        self.log('Getting the network.')
        edges = self.subnetwork(gene_list=gene_list, sig_only=sig_only,
        min_distance=min_distance, names_as_index=False,
        names_as_cols=True)
        del self.coex
    
        for source,target,score,significant,distance in subnet.itertuples(index=False):
            net['edges'].append(
                {'data':{
                    'source': source,
                    'target' : target,
                    'score' : score,
                    'distance' : fix_val(distance)
                }}
            )
        # Return the correct output
        net = {'elements' : net}
        if filename:
            with open(filename,'w') as OUT:
                print(json.dumps(net),file=OUT)
                del net
                self.coex = self._bcolz('coex')
        else:
            net = json.dumps(net)
            self.coex = self._bcolz('coex')
            return net


    def mcl(self, gene_list=None, I=2.0, scheme=7, min_distance=None,
            min_cluster_size=0, max_cluster_size=10e10):
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
            min_distance : int (default: None)
                The minimum distance between genes for which to consider
                co-expression interactions. This filters out cis edges.
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
        self.to_dat(
            filename=tmp.name, gene_list=gene_list,
            min_distance=min_distance, sig_only=True
        )
        # build the mcl command
        cmd = "mcl {} --abc -scheme {} -I {} -o -".format(tmp.name, scheme, I)
        self.log("running MCL: {}", cmd)
        try:
            p = Popen(cmd, stdout=PIPE, stderr=sys.stderr, shell=True)
            self.log('waiting for MCL to finish...')
            sout = p.communicate()[0]
            p.wait()
            self.log('MCL done, Reading results.')
            if p.returncode==0:
                # Filter out cluters who are smaller than the min size
                return list(
                    filter(
                        lambda x: len(x) > min_cluster_size and len(x) < max_cluster_size,
                        # Generate ids from the refgen
                        [ self.refgen.from_ids([gene.decode('utf-8') \
                                for gene in line.split()]) \
                                for line in sout.splitlines()
                        ]

                    )
                )
            else:
                raise ValueError( "MCL failed: return code: {}".format(p.returncode))
        except FileNotFoundError as e:
            self.log('Could not find MCL in PATH. Make sure its installed and shell accessible as "mcl".')

    def local_degree(self, gene_list, trans_locus_only=False):
        '''
            Returns the local degree of a list of genes

            gene_list : iterable (co.Locus object)
                a list of genes for which to retrieve local degree for. The
                genes must be in the COB object (of course)
            trans_locus_only : bool (default: False)
                only count edges if they are from genes originating from
                different loci. Each gene MUST have 'parent_locus' set in
                its attr object.

        '''
        subnetwork = self.subnetwork(
            gene_list, sig_only=True, trans_locus_only=trans_locus_only
        )
        if trans_locus_only:
            subnetwork = subnetwork.ix[subnetwork.trans]
        local_degree = pd.DataFrame(
            list(Counter(
                chain(*subnetwork.index.get_values())
            ).items()),
            columns=['Gene', 'Degree']
        ).set_index('Gene')
        # We need to find genes not in the subnetwork and add them as degree 0
        # The code below is ~optimized~
        # DO NOT alter unless you know what you're doing :)
        degree_zero_genes = pd.DataFrame(
            [(gene.id, 0) for gene in gene_list if gene.id not in local_degree.index],
            columns=['Gene', 'Degree']
        ).set_index('Gene')
        return pd.concat([local_degree, degree_zero_genes])

    def global_degree(self, gene_list, trans_locus_only=False):
        '''
            Returns the global degree of a list of genes
        '''
        try:
            if isinstance(gene_list, Locus):
                if trans_locus_only:
                    raise ValueError('Cannot calculate cis degree on one gene.')
                return self.degree.ix[gene_list.id].Degree
            else:
                degree = self.degree.ix[[x.id for x in gene_list]].fillna(0)
                if trans_locus_only:
                    degree = degree - self.cis_degree(gene_list)
                return degree
        except KeyError as e:
            return 0

    def cis_degree(self, gene_list):
        '''
            Returns the number of cis interactions for each gene in the gene
            list. **each gene object MUST have its 'parent_locus' attr set!!

            gene_list : iterable of co.Gene
        '''
        subnetwork = self.subnetwork(
            gene_list, sig_only=True, trans_locus_only=True
        )
        # Invert the trans column
        subnetwork['cis'] = np.logical_not(subnetwork.trans)
        subnetwork = subnetwork.ix[subnetwork.cis]
        local_degree = pd.DataFrame(
            list(Counter(
                chain(*subnetwork.index.get_values())
            ).items()),
            columns=['Gene', 'Degree']
        ).set_index('Gene')
        # We need to find genes not in the subnetwork and add them as degree 0
        # The code below is ~optimized~
        # DO NOT alter unless you know what you're doing :)
        degree_zero_genes = pd.DataFrame(
            [(gene.id, 0) for gene in gene_list if gene.id not in local_degree.index],
            columns=['Gene', 'Degree']
        ).set_index('Gene')
        return pd.concat([local_degree, degree_zero_genes])


    def locality(self, gene_list, iter_name=None, include_regression=False):
        '''
            Computes the merged local vs global degree table

            Parameters
            ----------
            gene_list : iterable of camoco.Loci
                A list or equivalent of loci
            iter_name : object (default: none)
                This will be added as a column. Useful for
                generating bootstraps of locality and keeping
                track of which one a row came from after catting
                multiple bootstraps together.
            include_regression : bool (default: False)
                Include the OLS regression residuals and fitted values
                on local ~ global.

            Returns
            -------
            A pandas DataFrame with local, global and residual columns
            based on linear regression of local on global degree.

        '''
        global_degree = self.global_degree(gene_list)
        local_degree = self.local_degree(gene_list)
        degree = global_degree.merge(
            local_degree,left_index=True,right_index=True
        )
        degree.columns = ['global', 'local']
        degree = degree.sort_values(by='global')
        if include_regression:
            # Add the regression lines
            ols = sm.OLS(degree['local'], degree['global']).fit()
            degree['resid'] = ols.resid
            degree['fitted'] = ols.fittedvalues
            degree = degree.sort_values(by='resid',ascending=False)
        if iter_name is not None:
            degree['iter_name'] = iter_name
        return degree



    ''' ----------------------------------------------------------------------
        Plotting Methods
    '''

    def plot(self, filename=None, genes=None,accessions=None,
             gene_normalize=True, raw=False,
             cluster_method='mcl', include_accession_labels=None,
             include_gene_labels=None, avg_by_cluster=False,
             min_cluster_size=10, cluster_accessions=True):
        '''
            Plots a heatmap of genes x expression.

            Parameters
            ----------
            filename : str 
                If specified, figure will be written to output filename
            genes : co.Locus iterable (default: None)
                An iterable of genes to plot expression for
            accessions : iterable of str
                An iterable of strings to extract for expression values.
                Values must be a subset of column values in expression matrix
            gene_normalize: bool (default: True)
                normalize gene values in heatmap to show expression patterns.
            raw : bool (default: False)
                If true, raw expression data will be used. Default is to use
                the normailzed, QC'd data.
            cluster_method : str (default: mcl)
                Specifies how to organize the gene axis in the heatmap. If
                'mcl', genes will be organized by MCL cluster. If 'leaf', 
                genes will be organized based on hierarchical clustering,
                any other value will result in genes to be in sorted order.
            include_accession_labels : bool (default: None)
                Force the rendering of accession labels. If None, accession 
                lables will be included as long as there are less than 30.
            include_gene_lables : bool (default: None)
                Force rendering of gene labels in heatmap. If None, gene
                labels will be rendered as long as there are less than 100.
            avg_by_cluster : bool (default: False)
                If True, gene expression values will be averaged by cluster
                showing a single row per cluster.
            min_cluster_size : int ( default: 10)
                If avg_by_cluster, only cluster sizes larger than min_cluster_size
                will be included.
            cluster_accessions : bool (default: True)
                If true, accessions will be clustered
            plot_dendrogram : bool (default: True)
                If true, dendrograms will be plotted

            Returns
            -------
            a populated matplotlib figure object

        '''
        # Get leaves of genes
        dm = self.expr(genes=genes,accessions=accessions,
                raw=raw,gene_normalize=gene_normalize)
        if cluster_method == 'leaf':
            order = self._bcolz('leaves').sort_values(by='index').index.values
        elif cluster_method == 'mcl':
            order = self._bcolz('clusters').loc[dm.index].\
                    fillna(np.inf).sort_values(by='cluster').index.values
        else:
            # No cluster order
            order = dm.index
        # rearrange expression by leaf order
        dm = dm.loc[order, :]

        # Optional Average by cluster
        if avg_by_cluster == True:
            # Extract clusters
            dm = self.clusters.groupby('cluster').\
                    filter(lambda x: len(x) >= min_cluster_size).\
                    groupby('cluster').\
                    apply(lambda x: self.expr(genes=self.refgen[x.index]).mean()).\
                    apply(lambda x: (x-x.mean())/x.std() ,axis=1)
            if len(dm) == 0:
                self.log.warn('No clusters larger than {} ... skipping',min_cluster_size)
                return None
        # Get leaves of accessions
        if cluster_accessions:
            accession_pccs = (1 - PCCUP.pair_correlation(
                np.ascontiguousarray(
                    # PCCUP expects floats
                    self._expr.as_matrix().T.astype('float')
                )
            ))
            link = linkage(1-accession_pccs, method='complete')
            accession_dists = leaves_list(link)
            # Order by accession distance
            dm = dm.loc[:,dm.columns[accession_dists]]


        # Save plot if provided filename
        fig = plt.figure(
            facecolor='white',
            figsize=(20,20)
        )
        ax = fig.add_subplot(111)
        nan_mask = np.ma.array(dm, mask=np.isnan(dm))
        cmap = self._cmap
        cmap.set_bad('grey', 1.0)
        vmax = max(np.nanmin(abs(dm)), np.nanmax(abs(dm)))
        vmin = vmax*-1
        im = ax.matshow(dm, aspect='auto', cmap=cmap, vmax=vmax, vmin=vmin)
        # Intelligently add labels
        if (include_accession_labels is None and len(dm.columns) < 30) \
            or include_accession_labels == True:
                ax.set(
                    xticks=np.arange(len(dm.columns)),
                    xticklabels=dm.columns.values
                )
                ax.set_xticklabels(dm.columns, rotation=90)
        if (include_gene_labels is None and len(dm.index) < 100) \
            or include_gene_labels == True:
                ax.set(
                    yticks=np.arange(len(dm.index)),
                    yticklabels=dm.index.values
                )
        fig.colorbar(im)
        # Save if you wish
        if filename is not None:
            plt.savefig(filename,dpi=300)
            plt.close()
        return fig



    def plot_scores(self, filename=None, pcc=True, bins=50):
        '''
            Plot the histogram of PCCs.

            Parameters
            ----------
            filename : str (default: None)
                The output filename, if none will return the matplotlib object
            pcc : bool (default:True)
                flag to convert scores to pccs
            bins : int (default: 50)
                the number of bins in the histogram
        '''
        fig = plt.figure(figsize=(8, 6))
        # grab the scores only and put in a
        # np array to save space (pandas DF was HUGE)
        scores = self.coex.score.dropna().values
        if pcc:
            self.log('Transforming scores')
            scores = (scores * float(self._global('pcc_std'))) \
                + float(self._global('pcc_mean'))
            # Transform Z-scores to pcc scores (inverse fisher transform)
            scores = np.tanh(scores)
        plt.hist(scores, bins=bins)
        plt.xlabel('PCC') if pcc else plt.xlabel('Z-Score')
        plt.ylabel('Freq')
        if filename is not None:
            plt.savefig(filename)
            plt.close()
        else:
            return fig


    def plot_locality(self, gene_list, bootstraps=10,
                      num_windows=100, sd_thresh=2):
        '''
            Make a fancy locality plot.
        '''
        # Generate a blank fig
        fig, ax = plt.subplots(figsize=(8, 6))
        fig.hold(True)
        # Y axis is local degree (what we are TRYING to predict)
        degree = self.locality(gene_list).sort('global')
        ax.set_ylim(0, max(degree['local']))
        ax.set_xlim(0, max(degree['global']))
        if bootstraps > 0:
            bs = pd.concat(
                [self.locality(
                    self.refgen.bootstrap_candidate_genes(gene_list)
                ) for x in range(10)]
            ).sort('global')
            ax.set_ylim(0, max(bs['local']))
            ax.set_xlim(0, max(bs['global']))
            plt.plot(bs['global'], bs['local'], 'ro', alpha=0.05, label='Bootstraps')
        # Plot the bootstraps and the empirical
        plt.plot(degree['global'], degree['local'], 'bo', label='Empirical')
        emp_ols = sm.OLS(degree['local'], degree['global']).fit()
        ax.plot(degree['global'], emp_ols.fittedvalues, 'k:', label='Empirical OLS')

        if bootstraps > 0:
            # Get the OLS
            bs_ols = sm.OLS(bs['local'], bs['global']).fit()
            bs['resid'] = bs_ols.resid
            bs['fitted'] = bs_ols.fittedvalues
            ax.plot(bs['global'], bs_ols.fittedvalues, 'g--', label='bootstrap OLS')
            # Do lowess on the residuals
            # We only care about windows within the empirical part
            window_tick = len(bs)/num_windows
            bs['window'] = [int(x/window_tick) for x in range(len(bs))]
            # get std for each window
            win_std = bs.groupby('window').apply(lambda df: df['resid'].std()).to_dict()
            bs['std_envelope'] = [win_std[x] for x in bs.window.values]
            # Plot confidence intervals
            prstd, iv_l, iv_u = wls_prediction_std(bs_ols)
            ax.plot(bs['global'], iv_u, 'g--', label='conf int.')
            ax.plot(bs['global'], iv_l, 'g--')
            # plot the
            ax.plot(
                bs['global'], bs['fitted']+(sd_thresh*bs['std_envelope']), 'r--'
                , label='{} s.d. envelope'.format(sd_thresh)
            )
            ax.plot(bs['global'], bs['fitted']-(sd_thresh*bs['std_envelope']), 'r--')
        ax.set_xlabel('Number Global Interactions')
        ax.set_ylabel('Number Local Interactions')
        legend = ax.legend(loc='best')
        return plt

    def compare_degree(self, obj, diff_genes=10, score_cutoff=3):
        '''
            Compares the degree of one COB to another.

            Parameters
            ----------
            obj : COB instance
                The object you are comparing the degree to.
            diff_genes : int (default: 10)
                The number of highest and lowest different
                genes to report
            score_cutoff : int (default: 3)
                The edge score cutoff used to called
                significant.
        '''
        self.log("Comparing degrees of {} and {}", self.name, obj.name)

        # Put the two degree tables in the same table
        lis = pd.concat(
            [self.degree.copy(), obj.degree.copy()],
            axis=1, ignore_index=True
        )

        # Filter the table of entries to ones where both entries exist
        lis = lis[(lis[0] > 0) & (lis[1] > 0)]
        delta = lis[0] - lis[1]

        # Find the stats beteween the two sets,
        # and the genes with the biggest differences
        delta.sort(ascending=False)
        highest = sorted(
            list(dict(delta[:diff_genes]).items()),
            key=lambda x: x[1], reverse=True
        )
        lowest = sorted(
            list(dict(delta[-diff_genes:]).items()),
            key=lambda x: x[1], reverse=False
        )
        ans = {
            'correlation_between_cobs':lis[0].corr(lis[1]),
            'mean_of_difference':delta.mean(),
            'std_of_difference':delta.std(),
            ('bigger_in_'+self.name):highest,
            ('bigger_in_'+obj.name):lowest
        }

        return ans


    ''' ----------------------------------------------------------------------
            Internal Methods
    '''

    def _calculate_coexpression(self, significance_thresh=3):
        '''
            Generates pairwise PCCs for gene expression profiles in self._expr.
            Also calculates pairwise gene distance.
        '''
        # Drop the old Coex table to save memory
        try:
            del self.coex
        except AttributeError:
            pass

        # 1. Calculate the PCCs
        self.log("Calculating Coexpression")
        pccs = (1 - PCCUP.pair_correlation(
            np.ascontiguousarray(
                # PCCUP expects floats
                self._expr.as_matrix().astype('float')
            )
        ))
        
        self.log("Applying Fisher Transform")
        pccs[pccs >= 1.0] = 0.9999999
        pccs[pccs <= -1.0] = -0.9999999
        pccs = np.arctanh(pccs)
        gc.collect();
        
        self.log("Calculating Mean and STD")
        # Sometimes, with certain datasets, the NaN mask overlap
        # completely for the two genes expression data making its PCC a nan.
        # This affects the mean and std fro the gene.
        pcc_mean = np.ma.masked_array(pccs, np.isnan(pccs)).mean()
        self._global('pcc_mean', pcc_mean)
        gc.collect()
        pcc_std = np.ma.masked_array(pccs, np.isnan(pccs)).std()
        self._global('pcc_std', pcc_std)
        gc.collect()
        
        # 2. Calculate Z Scores
        self.log("Finding adjusted scores")
        pccs = (pccs-pcc_mean)/pcc_std
        gc.collect()
        
        # 3. Build the dataframe
        self.log("Build the dataframe")
        tbl = pd.DataFrame(
            pccs, 
            index=np.arange(len(pccs)), 
            columns=['score'], 
            copy=False
        )
        del pccs
        gc.collect()
        
        # 3. Calculate Gene Distance
        self.log("Calculating Gene Distance")
        tbl['distance'] = self.refgen.pairwise_distance(
            gene_list=self.refgen.from_ids(self._expr.index))
        gc.collect()
        
        # 4. Assign significance
        self.log("Thresholding Significant Network Interactions")
        self._global('significance_threshold', significance_thresh)
        tbl['significant'] = tbl['score'] >= significance_thresh
        gc.collect()
        
        # 5. Store the table
        self.log("Storing the coex table")
        self._bcolz('coex', df=tbl)
        del tbl
        gc.collect()
        
        # 6. Load the new table into the object
        self.coex = self._bcolz('coex')
        self.log("Done")
        return self


    def _calculate_degree(self):
        '''
            Calculates degrees of genes within network. Stores
            them in our feather store.
        '''
        self.log('Building Degree')
        # Get significant expressions and dump coex from memory for time being
        # Generate a df that starts all genes at 0
        names = self._expr.index.values
        self.degree = pd.DataFrame(0,index=names,columns=['Degree'])
        # Get the index and find the counts
        self.log('Calculating Gene degree')
        sigs = self.coex[self.coex.significant]
        sigs = sigs.index.values
        sigs = PCCUP.coex_expr_index(sigs, len(self._expr.index.values))
        sigs = list(Counter(chain(*sigs)).items())
        if len(sigs) > 0:
            # Translate the expr indexes to the gene names
            for i,degree in sigs:
                self.degree.ix[names[i]] = degree
        self._bcolz('degree', df=self.degree)
        # Cleanup
        del sigs 
        del names
        gc.collect()
        return self
        
    def _calculate_leaves(self):
        '''
            This calculates the leaves of the dendrogram from the coex
        '''
        # We need to recreate the original PCCs
        self.log('Calculating Leaves')
        if len(self.coex) == 0:
            raise ValueError('Cannot calculate leaves without coex')
        pcc_mean = float(self._global('pcc_mean'))
        pcc_std  = float(self._global('pcc_std'))
        
        # Get score column and dump coex from memory for time being
        dists = self.coex.score
        del self.coex
        
        # Subtract pccs from 1 so we do not get negative distances
        dists = (dists * pcc_std) + pcc_mean
        dists = np.tanh(dists)
        dists = 1 - dists
        gc.collect()
        
        # Find the leaves from hierarchical clustering
        self.log("Finding the leaves")
        dists = leaves_list(linkage(dists, method='single'))
        gc.collect()
        
        # Put them in a dataframe and stow them
        self.leaves = pd.DataFrame(dists,index=self._expr.index,columns=['index'])
        self._bcolz('leaves', df=self.leaves)
        
        # Cleanup and reinstate the coex table
        del dists; gc.collect();
        self.coex = self._bcolz('coex')
        return self
    
    def _calculate_clusters(self):
        '''
            Calculates global clusters
        '''
        clusters = self.mcl()
        self.log('Building cluster dataframe')
        names = self._expr.index.values
        self.clusters = pd.DataFrame(np.nan,index=names,columns=['cluster'])
        if len(clusters) > 0:
            self.clusters = pd.DataFrame(
                data=[(gene.id, i) for i, cluster in enumerate(clusters) \
                        for gene in cluster],
                columns=['Gene', 'cluster']
            ).set_index('Gene')
            self._bcolz('clusters', df=self.clusters)
        self.log('Creating Cluster Ontology')
        terms = []
        for i,x in enumerate(self.clusters.groupby('cluster')):
            genes = self.refgen[x[1].index.values]
            terms.append(Term(
                'MCL{}'.format(i),
                desc='{} MCL Cluster {}'.format(self.name,i),
                loci = genes
            ))
        self.MCL = Ontology.from_terms(
            terms,
            '{}MCL'.format(self.name),
            '{} MCL Clusters'.format(self.name),
            self.refgen)
        self.log('Finished finding clusters')
        return self
    
    def _coex_concordance(self, gene_a, gene_b, maxnan=10):
        '''
            This is a sanity method to ensure that the pcc calculated
            directly from the expr profiles matches the one stored in
            the database
        '''
        expr_a = self.expr_profile(gene_a).values
        expr_b = self.expr_profile(gene_b).values
        mask = np.logical_and(np.isfinite(expr_a), np.isfinite(expr_b))
        if sum(mask) < maxnan:
            # too many nans to reliably calculate pcc
            return np.nan
        r = pearsonr(expr_a[mask], expr_b[mask])[0]
        # fisher transform it
        z = np.arctanh(r-0.0000001)
        # standard normalize it
        z = (z - float(self._global('pcc_mean'))) \
            / float(self._global('pcc_std'))
        return z

    ''' -----------------------------------------------------------------------
            Class Methods -- Factory Methods
    '''
    @classmethod
    def create(cls, name, description, refgen):
        self = super().create(name, description, refgen)
        self._bcolz('gene_qc_status', df=pd.DataFrame())
        self._bcolz('accession_qc_status', df=pd.DataFrame())
        self._bcolz('coex', df=pd.DataFrame())
        self._bcolz('degree', df=pd.DataFrame())
        self._bcolz('mcl_cluster', df=pd.DataFrame())
        self._bcolz('leaves', df=pd.DataFrame())
        self._expr_index = defaultdict(
            lambda: None,
            {gene:index for index, gene in enumerate(self._expr.index)}
        )
        return self

    @classmethod
    def from_Expr(cls, expr):
        '''
            Create a COB instance from an camoco.Expr (Expression) instance.
            A COB inherits all the methods of a Expr instance and implements
            additional coexpression specific methods. This method accepts an
            already build Expr instance and then performs the additional
            computations needed to build a full fledged COB instance.

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
        self._calculate_leaves()
        self._calculate_clusters()
        return self

    @classmethod
    def from_DataFrame(cls, df, name, description,
                       refgen, rawtype=None, **kwargs):
        '''
            The method will read the table in (as a pandas dataframe),
            build the Expr object passing all keyword arguments in **kwargs
            to the classmethod Expr.from_DataFrame(...). See additional
            **kwargs in COB.from_Expr(...)

            Parameters
            ----------
            df : pandas.DataFrame
                A Pandas dataframe containing the expression information.
                Assumes gene names are in the index while accessions
                (experiments) are stored in the columns.
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
                This is noted here to reinforce the impotance of the rawtype
                passed to camoco.Expr.from_DataFrame. See docs there
                for more information.
            **kwargs : key,value pairs
                additional parameters passed to subsequent methods.
                (see Expr.from_DataFrame)

        '''
        # Create a new Expr object from a data frame
        expr = super().from_DataFrame(
            df, name, description, refgen, rawtype, **kwargs
        )
        return cls.from_Expr(expr)

    @classmethod
    def from_table(cls, filename, name, description,
                   refgen, rawtype=None, sep='\t', index_col=None, **kwargs):
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
                This is noted here to reinforce the importance of the rawtype
                passed to camoco.Expr.from_DataFrame. See docs there for
                more information.
            sep : str (default: \\t)
                Specifies the delimiter of the file referenced by the
                filename parameter.
            index_col : str (default: None)
                If not None, this column will be set as the gene index
                column. Useful if there is a column name in the text file
                for gene names.
            **kwargs : key value pairs
                additional parameters passed to subsequent methods.

            Returns
            -------
                a COB object
        '''
        df = pd.read_table(
            filename,
            sep=sep,
            compression='infer',
            index_col=index_col
        )
        return cls.from_DataFrame(
            df, name, description, refgen,
            rawtype=rawtype,**kwargs
        )


    '''
        Unimplemented ---------------------------------------------------------------------------------
    '''

    def next_neighbors(self, gene_list):
        ''' returns a list containing the strongest connected neighbors '''
        raise NotImplementedError()

    def neighborhood(self, gene_list):
        ''' Input: A gene List
            Output: a Dataframe containing gene ids which have at least
            one edge with another gene in the input list. Also returns
            global degree
        '''
        raise NotImplementedError()

    def lcc(self, gene_list, min_distance=None):
        ''' returns an igraph of the largest connected component in graph '''
        raise NotImplementedError()

    def seed(self, gene_list, limit=65):
        ''' Input: given a set of nodes, add on the next X strongest connected
            nodes '''
        raise NotImplementedError()

    def graph(self, gene_list, min_distance=None):
        ''' Input: a gene list
            Output: a iGraph object '''
        raise NotImplementedError()

    def coordinates(self, gene_list, layout=None):
        ''' returns the static layout, you can change the stored layout by
            passing in a new layout object. If no layout has been stored or a gene
            does not have coordinates, returns (0, 0) for each mystery gene'''
        raise NotImplementedError()

def fix_val(val):
    if isinf(val):
        return -1
    if np.isnan(val):
        # because Fuck JSON
        return "null"
    else:
        return val
