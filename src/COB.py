#!/usr/bin/python

from __future__ import print_function

import matplotlib as plt
#plt.use('Agg')
import numpy as np
import time
import igraph as ig
import sys
import pandas as pd

from scipy.stats import hypergeom
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform

from COBDatabase import *
from COBDataset import *
from COBLocus import *
from COBGene import *
from COBQTL import *
from COBSNP import *
from Genome import *
from Chrom import *

class Log(object):
    log_file = sys.stderr
    def init(self,log_file=None):
        if log_file:
            self.log_file = open(log_file,'w')
    def log(self,*args):
        ''' shared object logging '''
        print(time.ctime(),'-',*args,file=self.log_file)

def available_datasets():
        c = COBDatabase()
        return c.query('''SELECT * FROM datasets''')

def edit_dataset(dataset_name, field, new_field):
        c = COBDatabase()
        assert field in ["name","description","FPKM",'gene_build'],'nope'
        try:
            c.execute("UPDATE datasets SET {} = '{}' WHERE name = {};",field,new_field,dataset_name)
        except Exception as e:
            c.log("{} was not changeable in {}",field,dataset_name)

class COB(COBDatabase,Log):
    ''' 
        This class implements the basic interface to the COB Database. It provides routine
        analysis functions for basic network queries.
    '''
    def __init__(self,network=None):
        ''' Ititialized a COB object based on a dataset Name '''
        super(COB,self).__init__()
        # HouseKeeping
        try:
            (self.DatasetID, self.name, self.desc,
             self.FPKM, self.build,
             self.date_added
            ) = self.fetchone('''
                SELECT * FROM datasets WHERE name = '{}' ''',network
            )
        except Exception as e:
            self.log("Dataset not found.")
            

    @property        
    def geneIDS(self):
        return self.query('''
            SELECT DISTINCT(GrameneID) FROM expression 
            JOIN genes ON expression.GeneID = genes.ID
            WHERE DatasetID = {}''',self.DatasetID).GrameneID

    @property
    def genes(self):
        return [COBGene(x,self.build) for x in self.query('''
            SELECT DISTINCT(GrameneID) FROM expression 
            JOIN genes ON expression.GeneID = genes.ID 
            WHERE DatasetID = {}''',self.DatasetID).GrameneID]
    
    @property
    def accessions(self):
        return self.query('''
            SELECT name FROM accessions
            WHERE DatasetID = {}
         ''',self.DatasetID).name

    @property
    def coex_table_name(self):
        return "coex_{}".format(self.name.replace(' ','_'))

    def neighbors(self,gene_list):
        ''' 
            Input : a list of COBGene Objects
            Output: Returns the neighbors for a list of genes as a DataFrame
            Columns returned: query_name, queryID, neighbor_name, neighbor, score
        '''
        return pd.concat([
            self.query('''
            SELECT query.GrameneID as query_name, gene_b as query, 
                   neighbor.GrameneID as neighbor_name, gene_a as neighbor, score
            FROM {} coex
            JOIN genes neighbor ON coex.gene_a = neighbor.ID
            JOIN genes query    ON coex.gene_b = query.ID
            WHERE significant = 1 AND gene_b IN ({})''',
                 self.coex_table_name, ",".join(map(str,[x.ID for x in gene_list]))),
            self.query('''
            SELECT query.GrameneID as query_name, gene_a as query, 
                   neighbor.GrameneID as neighbor_name, gene_b as neighbor, score
            FROM {} coex 
            JOIN genes neighbor ON coex.gene_b = neighbor.ID
            JOIN genes query    ON coex.gene_a = query.ID
            WHERE significant = 1 AND gene_a IN ({})''', 
                self.coex_table_name, ','.join(map(str,[x.ID for x in gene_list])))
        ])

    def num_neighbors(self,gene_list):
        ''' 
            Input: a list of COBGene Object
            Output: Dataframe containing the number of significant global
                    and local neighbors for each gene in list
        '''
        neighbors = self.neighbors(gene_list)
        def get_nums(df):
            global_length =  len(df)
            local_length = len(set(neighbors['query']).intersection(set(df.neighbor)))
            return pd.Series({
                'global': global_length,
                'local' : local_length 
                })
        return neighbors.groupby('query_name').apply(lambda x : get_nums(x))

    def neighbors_score(self,gene_list):
        ''' returns a series containing the strongest connected neighbors '''
        neighbors = self.neighbors(gene_list)
        if len(neighbors) > 0:
            scores = neighbors.groupby('neighbor_name').score.sum()
            scores.sort(ascending=False)
        return scores

    def neighborhood(self,gene_list):
        ''' Input: A gene List
            Output: a Dataframe containing gene ids which have at least one edge
                    with another gene in the input list. Also returns global degree '''
        if len(gene_list) == 0:
            return []
        self.log("Analyzing {} Genes for {} Network",len(gene_list),self.name)
        self.log("Found {} genes in COB",len(gene_list))
        neighbors = self.num_neighbors(gene_list)
        local = neighbors[neighbors.local >= 1]
        self.log("Found {} genes in subnetwork",len(local))
        return(local)

    def subnetwork(self,gene_list):
        ''' Input: a gene list
            Output: a dataframe containing all edges between genes within list
        '''
        if len(gene_list) == 0:
            return pd.DataFrame()
        subnet = self.query(''' 
            SELECT source.GrameneID source, gene_a, 
                target.GrameneID as target, gene_b, score 
            FROM {} coex 
            JOIN genes source ON source.ID = gene_a
            JOIN genes target ON target.ID = gene_b
            WHERE significant = 1 
            AND gene_a IN ({})
            AND gene_b IN ({}) ''', 
                self.coex_table_name, 
                ','.join(map(str,[x.ID for x in gene_list])), 
                ",".join(map(str,[x.ID for x in gene_list]))
        )
        return subnet
 
    def density(self,gene_list):
        ''' calculate the density of the ENTIRE network between a set of genes 
            Input: A list of COBGene Objects
            Output: a scalar value representing the density of genes
        '''
        if len(gene_list) == 0:
            return 0
        scores = self.query("SELECT score FROM {} coex WHERE gene_a IN ({}) AND gene_b IN ({})",
            self.coex_table_name,
            ",".join(map(str,[x.ID for x in gene_list])),
            ",".join(map(str,[x.ID for x in gene_list]))
        ).score
        if len(scores) == 0:
            return 0
        return (scores.mean()/(1/np.sqrt(len(scores))))
 
    def seed(self, gene_list, max_show = 65): 
        ''' Input: given a set of nodes, add on the next X strongest connected nodes ''' 
        if len(gene_list) == 0:
            return []
        neighbors = self.neighbors_score(gene_list)
        seed_set =  set(self._gene_ids(gene_list)).union(neighbors.index[0:min(len(neighbors),max_show)])
        return seed_set

    def gene_expr_vals(self,gene_list):
        ''' Input: A list of COBGenes
            Output: A dataframe containing normalized expression values for gene_list
        '''
        expr_vals = self.query(
         '''SELECT 
            GrameneID, expression.GeneID, expression.AccessionID, expression.DatasetID,
                 value, meanExpr, stdExpr, name, class,description
            FROM expression
            JOIN avgexpr ON expression.GeneID = avgexpr.GeneID 
                AND expression.DatasetID = avgexpr.DatasetID
            JOIN accessions ON expression.AccessionID = accessions.ID 
                AND expression.DatasetID = accessions.DatasetID
            JOIN genes ON expression.GeneID = genes.ID
            WHERE expression.GeneID IN ({})
            AND   expression.DatasetID = {}''',
            ",".join(map(str,[x.ID for x in gene_list])), 
            self.DatasetID
        )
        def zscore(group):
            group.index = group.name
            norm = (group.value-group.value.mean())/group.value.std()
            return norm
        return pd.DataFrame({gid:zscore(group) for gid,group in expr_vals.groupby('GrameneID')})

    def graph(self,gene_list):
        ''' Input: a gene list
            Output: a iGraph object '''
        if len(gene_list) == 0:
            return []
        # retrieve all first neighbors
        edges = self.subnetwork(gene_list)
        nodes = list(set(edges.source).union(edges.target))
        graph = ig.Graph()
        for node in nodes:
            graph.add_vertex(name=node,label=node)
        for edge in edges[['source','target','score']].itertuples():
            graph.add_edge(source=graph.vs.find(name=edge[1]),target=graph.vs.find(name=edge[2]),weight=edge[3])
        return graph 
        

    def heatmap(self,dm,filename=None):

        D = np.array(dm) 
        D1 = squareform(pdist(dm, metric='euclidean'))
        D2 = squareform(pdist(dm.T, metric='euclidean'))
        
        f = plt.pylab.figure(figsize=(16, 16),facecolor='white')

        # add first dendrogram
        ax1 = f.add_axes([0.09, 0.1, 0.2, 0.6])
        ax1.set_frame_on(False)
        Y = linkage(D1, method='complete')
        Z1 = dendrogram(Y, orientation='right')
        ax1.set_xticks([])
        ax1.set_yticks([])

        # add second dendrogram
        ax2 = f.add_axes([0.3, 0.71, 0.5, 0.2])
        ax2.set_frame_on(False)
        Y = linkage(D2, method='complete')
        Z2 = dendrogram(Y)
        ax2.set_xticks([])
        ax2.set_yticks([])

        # add matrix plot
        axmatrix = f.add_axes([0.3, 0.1, 0.5, 0.6])
        idx1 = Z1['leaves']
        idx2 = Z2['leaves']
        D = D[idx1, :]
        D = D[:, idx2]
        vmax = max(abs(D).min(),abs(D).max())
        vmin = vmax*-1
        im = axmatrix.matshow(D, aspect='auto', cmap=self.cmap,vmax=vmax,vmin=vmin)
        # Handle Axis Labels
        axmatrix.set_xticks(np.arange(D.shape[1]))
        axmatrix.xaxis.tick_bottom()
        axmatrix.tick_params(axis='x',labelsize='xx-small')
        axmatrix.set_xticklabels(dm.columns[idx2],rotation=90,ha='right')
        axmatrix.yaxis.tick_right()
        axmatrix.set_yticks(np.arange(D.shape[0]))
        axmatrix.set_yticklabels(dm.index[idx1])
        axmatrix.tick_params(axis='y',labelsize='x-small')

        axColorBar = f.add_axes([0.09,0.75,0.2,0.05])
        f.colorbar(im,orientation='horizontal',cax=axColorBar,ticks=np.arange(np.ceil(vmin),np.ceil(vmax),2))

        if filename:
            plt.pyplot.savefig(filename)
            plt.pyplot.close()

        return {'ordered' : D, 'rorder' : Z1['leaves'], 'corder' : Z2['leaves']} 

# Unimplemented Below Here ------------------------------------------
    def url(self,gene_list,size=65):
        gene_list = self.ids2genes(gene_list)
        return "http://lovelace.cs.umn.edu/cob?results_source={}&action=seed&q_list={}&neighborhood_size={}".format(
            self.network,
            ",".join(map(str,gene_list)),
            size
        )   
           
    def go_enrichment(self,gene_list,pval_cutoff = 0.05):
        if self.verbose:
            self.log("Calculating go enrichment for {}",",".join(map(str,gene_list)))
        if not gene_list:
            return pd.DataFrame()
        gene_ids = self.genes2ids(gene_list)
        go_terms = self.query('''SELECT DISTINCT(go_id) FROM mzn_gene_go_terms 
            JOIN mzn_go_terms ON go_id = term_id 
            WHERE gene_id IN ({}) ''',",".join(map(str,gene_ids)))
        if go_terms.empty:
            return pd.DataFrame()
        terms = []
        for id in go_terms.go_id.values:
            annotation = self.query("SELECT * FROM mzn_go_terms WHERE term_id = {}",id).iloc[0]
            genes_in_go_term = self.query('''SELECT gene_id FROM mzn_gene_go_terms WHERE go_id = {}''',id)
            num_genes_in_go_term = len(genes_in_go_term)
            overlap = set(genes_in_go_term.gene_id).intersection(set(gene_ids))
            num_genes_total = self.num_go_genes_total
            pval = hypergeom.sf(len(overlap),num_genes_total,num_genes_in_go_term,len(gene_ids))
            terms.append({
                "term_name" : annotation.term_name,
                "pvalue" : pval,
                "num_in_go" : num_genes_in_go_term,
                "num_overlap":len(overlap),
                "num_total" : num_genes_total,
                "num_in_cluster":len(gene_ids),
                "short_desc" : annotation.term_short_desc
            })
        terms = pd.DataFrame(terms).sort("pvalue",ascending=True)
        return terms[terms.pvalue <= pval_cutoff]
            
    def gene_summary(self,gene_list, gene_in=dict(), sep = "\n"):
        ''' Returns a table summarizing features of the gene list. kwargs include gene_in
            which is a dictionary containing keys which are name ofgroups of loci and values which
            are lists of loci and each will be translated to a column in the final dataframe'''
        tbl = pd.DataFrame()
        tbl['GRAMENEID'] = [x.gramene_id for x in gene_list]
        tbl['COMMON']    = [x.common_name for x in gene_list] 
        tbl['CHR'] = [x.chrom for x in gene_list]
        tbl['START'] = [x.start for x in gene_list]
        tbl['STOP'] = [x.end for x in gene_list]
        tbl['GO_TERM'] = [sep.join(x.go_terms.term_name) for x in gene_list]
        tbl['GO_DESC'] = [sep.join(x.go_terms.term_short_desc) for x in gene_list]
        tbl['ARAB_ORTH'] = [sep.join(x.arab_ortho.arab_name) for x in gene_list]
        for loci_name,locus_list in gene_in.items():
            tbl[loci_name] = [sep.join(
                map(str,filter(
                    lambda locus: str(locus.id) if gene in locus else None, locus_list 
                ))) for gene in gene_list
            ]
        return tbl


    ############################################################################
    # Helper Functions
    ############################################################################
    def log(self, string, *args, **kwargs):
        print('[LOG]',time.ctime(), "-", string.format(*args),file=self.log_file)

    def _gene_ids(self, gene_list, attr="ID"):
        ''' return a list of ids from a gene list '''
        return [getattr(x,attr) for x in gene_list]

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
        heatmap_cmap = plt.colors.LinearSegmentedColormap('my_colormap',heatmapdict,256)
        return heatmap_cmap

    def __repr__(self):
        return '''
        Name: {}
        Description: {}
        DatasetID: {}
        FPKM: {}
        Gene Build: {}
        Date Added: {}
        Number of Genes: {}
        Number of accessions: {}
        '''.format(self.name, self.desc, self.DatasetID, self.FPKM,
            self.build, self.date_added, len(self.geneIDS), len(self.accessions)) 



    def compare_to_dat(self,filename,sep="\t",score_cutoff=3):
        ''' Compare the number of genes with significant edges as well as degree with a DAT file '''
        self.log("Reading in {}",filename)
        ref_dat = pd.read_table(filename,sep=sep,names=['gene_a','gene_b','score'])
        self.log("Retrieving network edges...")
        dat = self.query('''SELECT source.GrameneID as gene_a, target.GrameneID as gene_b, score 
            FROM {} coex 
            JOIN genes source ON coex.gene_a = source.ID
            JOIN genes target ON coex.gene_b = target.ID
            WHERE score >= {};
            ''',self.coex_table_name,score_cutoff)
        # Get a list of genes in the dataset as well as reference
        genes = set(dat.gene_a).union(set(dat.gene_b))
        self.log("Found {} genes in {}",len(genes),self.name)
        ref_genes = set(ref_dat.gene_a).union(set(ref_dat.gene_b))
        self.log("Found {} genes in {}",len(ref_genes),filename)
        self.log("Found {} genes in intersection",len(genes.intersection(ref_genes)))
        if len(genes) != len(ref_genes):
            # Draw a venn diagram you lazy sac
            self.log("Genes not in {}: {}",self.name,"\n".join(ref_genes.difference(genes.intersection(ref_genes))))
            self.log("Genes not in {}: {}",filename,"\n".join(genes.difference(genes.intersection(ref_genes))))
        def compare_edges(gene):
            # Compare the degree of each gene with the reference dat
            ref_edges = ref_dat[(ref_dat.gene_a == gene) | (ref_dat.gene_b == gene)]
            edges     = dat[(dat.gene_a == gene)|(dat.gene_b == gene)]
            ref_neighbors = set(ref_edges.gene_a).union(set(ref_edges.gene_b))
            neighbors     = set(edges.gene_a).union(set(edges.gene_b))
            return (len(neighbors),len(ref_neighbors))
        return [compare_edges(gene) for gene in genes.intersection(ref_genes)]

