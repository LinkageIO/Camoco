#!/usr/bin/python3

import camoco.PCCUP as PCCUP

from .Camoco import Camoco
from .RefGen import RefGen
from .Locus import Locus,Gene
from .Expr import Expr
from .Tools import memoize

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

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

import statsmodels.api as sm
import sys
import json

from scipy.stats import pearsonr

class COB(Expr):
    def __init__(self, name):
        super().__init__(name=name)
        try:
            self.log('Loading coex table')
            self.coex = self.hdf5['coex']
        except KeyError as e:
            self.log("{} is empty ({})", name, e)
        try:
            self.log('Loading Global Degree')
            self.degree = self.hdf5['degree']
        except KeyError as e:
            self.log("{} is empty ({})", name, e)
        try:
            self.log('Loading Clusters')
            self.clusters = self.hdf5['clusters']
        except KeyError as e:
            self.log('Clusters not loaded for: {} ()', name, e)

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

            QC
            ------------------
            min expr level: {}
            max gene missing data: {}
            min single sample expr: {}
            max accession missing data: {}

            Clusters
            ------------------
            Num clusters (size > 3): {}


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
            sum(self.clusters.groupby('cluster').apply(len) > 3)
        ), file=file)

    def qc_gene(self):
        qc_gene = self.hdf5['qc_gene']
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

    def neighbors(self, gene, sig_only=True):
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
        neighbor_indices = PCCUP.coex_neighbors(gene_id, self.num_genes())
        edges = self.coex.iloc[neighbor_indices]
        if sig_only:
            return edges[edges.significant == 1]
        else:
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
        # Grab the indices in the original expression matrix
        ids = np.array([self._expr_index[gene_a.id], self._expr_index[gene_b.id]])
        # We need the number of genes
        num_genes = self.num_genes()
        index = PCCUP.coex_index(ids, num_genes)[0]
        return self.coex.iloc[index]

    def subnetwork(self, gene_list=None, sig_only=True, min_distance=None,
        filter_missing_gene_ids=True, trans_locus_only=False):
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

            Returns
            -------
            A pandas.DataFrame containing the edges. Columns
            include score, significant (bool), inter-genic distance,
            and

        '''
        if gene_list is None:
            df = self.coex
        else:
            ids = np.array([self._expr_index[x.id] for x in set(gene_list)])
            if filter_missing_gene_ids:
                # filter out the Nones
                ids = np.array(list(filter(None, ids)))
            num_genes = self.num_genes()
            # Grab the coexpression indices for the genes
            indices = PCCUP.coex_index(ids, num_genes)
            df = self.coex.iloc[indices]
        if min_distance is not None:
            df = df.loc[df.distance >= min_distance, :]
        if sig_only:
            df = df.loc[df.significant == 1, :]
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
        return df.copy()

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
        self.log("Found {} candidate genes", len(genes_list))
        # Extract the edges for the full set of genes
        edges = self.subnetwork(
            genes_list,
            min_distance=0,
            sig_only=False,
            trans_locus_only=True
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
        self.log("Found {} candidate genes", len(genes_list))
        # Get global and local degree for candidates
        gdegree = self.global_degree(genes_list, trans_locus_only=True)
        ldegree = self.local_degree(genes_list, trans_locus_only=True)
        # Merge the columns
        degree = ldegree.merge(gdegree,left_index=True,right_index=True)
        degree.columns = ['local', 'global']
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

    def to_dat(self, gene_list=None, filename=None, sig_only=False, min_distance=0):
        '''
            Outputs a .DAT file (see Sleipnir library)
        '''
        if filename is None:
            filename = self.name + '.dat'
        with open(filename, 'w') as OUT:
            self.log("Creating .dat file")
            self.subnetwork(
                gene_list, sig_only=sig_only, min_distance=min_distance
            )['score'].to_csv(OUT, sep='\t')
            self.log('Done')

    def to_graphml(self,file, gene_list=None,sig_only=True,min_distance=0):
        # Get all the graph edges (and the nodes implicitly)
        self.log('Getting the network.')
        edges = self.subnetwork(gene_list=gene_list, sig_only=sig_only, min_distance=min_distance).index.values

        # Build the NetworkX network
        self.log('Building the graph.')
        net = nx.Graph()
        net.add_edges_from(edges)

        # Print the file
        self.log('Writing the file.')
        nx.write_graphml(net,file)
        return

    def to_treeview(self, filename, cluster_method='mcl', gene_normalize=True):
        dm = self.expr(gene_normalize=gene_normalize)
        if cluster_method == 'leaf':
            order = self.hdf5['leaves'].sort('index').index.values
        elif cluster_method == 'mcl':
            order = self.hdf5['clusters'].loc[dm.index].\
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


    def to_json(self, gene_list=None, filename=None):
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
        # Now do the edges
        subnet = self.subnetwork(gene_list)
        subnet.reset_index(inplace=True)
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
        else:
            return json.dumps(net)


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
            self.log('...Done')
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
             include_gene_labels=None):
        '''
            Plots a heatmap of genes x expression.
        '''
        # Get leaves
        dm = self.expr(genes=genes,accessions=accessions,
                raw=raw,gene_normalize=gene_normalize)
        if cluster_method == 'leaf':
            order = self.hdf5['leaves'].sort('index').index.values
        elif cluster_method == 'mcl':
            order = self.hdf5['clusters'].loc[dm.index].\
                    fillna(np.inf).sort_values(by='cluster').index.values
        else:
            # No cluster order
            order = dm.index
        # rearrange expression by leaf order
        dm = dm.loc[order, :]
        # Save plot if provided filename
        fig = plt.figure(
            facecolor='white'
        )
        ax = fig.add_subplot(111)
        nan_mask = np.ma.array(dm, mask=np.isnan(dm))
        cmap = self._cmap
        cmap.set_bad('grey', 1.0)
        vmax = max(np.nanmin(abs(dm)), np.nanmax(abs(dm)))
        vmin = vmax*-1
        im = ax.matshow(dm, aspect='auto', cmap=cmap, vmax=vmax, vmin=vmin)
        if (include_accession_labels is None and len(dm.columns) < 30) \
            or include_accession_labels == True:
                ax.set(
                    xticks=np.arange(len(dm.columns)),
                    xticklabels=dm.columns.values
                )
                ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=90)
        if (include_gene_labels is None and len(dm.index) < 100) \
            or include_gene_labels == True:
                ax.set(
                    yticks=np.arange(len(dm.index)),
                    yticklabels=dm.index.values
                )
        fig.colorbar(im)
        # Save if you wish
        if filename is not None:
            plt.savefig(filename,figsize=(5,10))
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
        # Start off with a fresh set of genes we can pass to functions
        tbl = pd.DataFrame(
            list(itertools.combinations(self._expr.index.values, 2)),
            columns=['gene_a', 'gene_b']
        )

        # 1. Calculate the PCCs
        self.log("Calculating Coexpression")
        pccs = 1 - PCCUP.pair_correlation(
            np.ascontiguousarray(self._expr.as_matrix())
        )
        assert len(pccs) == len(tbl)
        tbl['score'] = pccs
        # correlations of 1 dont transform well, they cause infinities
        tbl.loc[tbl['score'] == 1, 'score'] = 0.99999999
        tbl.loc[tbl['score'] == -1, 'score'] = -0.99999999
        # Perform fisher transform on PCCs
        tbl['score'] = np.arctanh(tbl['score'])
        # Sometimes, with certain datasets, the NaN mask overlap
        # completely for the two genes expression data making its PCC a nan.
        # This affects the mean and std fro the gene.
        valid_scores = np.ma.masked_array(
            tbl['score'],
            np.isnan(tbl['score'])
        )

        # 2. Calculate Z Scores
        pcc_mean = valid_scores.mean()
        pcc_std = valid_scores.std()
        # Remember these so we can go back to PCCs
        self._global('pcc_mean', pcc_mean)
        self._global('pcc_std', pcc_std)
        tbl['score'] = (valid_scores-pcc_mean)/pcc_std

        # 3. Assign significance
        self._global('significance_threshold', significance_thresh)
        tbl['significant'] = pd.Series(
            list(tbl['score'] >= significance_thresh),
            dtype='int_'
        )

        # 4. Calculate Gene Distance
        self.log("Calculating Gene Distance")
        distances = self.refgen.pairwise_distance(
            gene_list=self.refgen.from_ids(self._expr.index)
        )
        assert len(distances) == len(tbl)
        tbl['distance'] = distances
        # Reindex the table to match genes
        self.log('Indexing coex table')
        tbl.set_index(['gene_a', 'gene_b'], inplace=True)

        # 5. Put table in the hdf5 store
        self._build_tables(tbl)
        self.log("Done")
        return self

    def _calculate_clusters(self):
        '''
            Calculates global clusters
        '''
        clusters = self.mcl()
        self.hdf5['clusters'] = pd.DataFrame(
            data=[(gene.id, i) for i, cluster in enumerate(clusters) \
                    for gene in cluster],
            columns=['Gene', 'cluster']
        ).set_index('Gene')
        self.hdf5.flush(fsync=True)
        self.clusters = self.hdf5['clusters']
        return self

    def _calculate_degree(self):
        '''
            Calculates degrees of genes within network. Stores
            them in our HDF5 store.
        '''
        self.log('Building Degree')
        self.hdf5['degree'] = pd.DataFrame(
            data=list(Counter(chain(
                *self.subnetwork(sig_only=True).index.get_values()
            )).items()),
            columns=['Gene', 'Degree']
        ).set_index('Gene')
        self.hdf5.flush(fsync=True)
        self.degree = self.hdf5['degree']
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
        # Subtract pccs from 1 so we do not get negative distances
        dists = 1 - np.tanh((self.coex.score * pcc_std)+pcc_mean)
        self.leaves = pd.DataFrame(
            leaves_list(linkage(dists, method='single')),
            index=self._expr.index,
            columns=['index']
        )
        # store the leaves
        self.hdf5['leaves'] = self.leaves
        self.hdf5.flush(fsync=True)
        return self


    def _build_tables(self, tbl):
        try:
            self.log("Building Database")
            ## HDF5 Store
            self.hdf5['coex'] = tbl
            self.log("Flushing Database")
            self.hdf5.flush(fsync=True)
            self.coex = self.hdf5['coex']
            self.log("Done")
        except Exception as e:
            self.log("Something bad happened:{}", e)
            raise

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
        self.hdf5['gene_qc_status'] = pd.DataFrame()
        self.hdf5['accession_qc_status'] = pd.DataFrame()
        self.hdf5['coex'] = pd.DataFrame()
        self.hdf5['degree'] = pd.DataFrame()
        self.hdf5['mcl_cluster'] = pd.DataFrame()
        self.hdf5['leaves'] = pd.DataFrame()
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
        return cls.from_DataFrame(
            pd.read_table(
                filename,sep=sep,
                compression='infer',
                index_col=index_col
            ),
            name,description,refgen,
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


