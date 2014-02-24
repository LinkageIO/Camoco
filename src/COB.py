#!/usr/bin/python

from __future__ import print_function
import pandas as pd
import pandas.io.sql as psql
import numpy as np
import sys
import os.path as path
import MySQLdb
import time
import igraph as ig

from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform

import matplotlib as plt

class COBDatabase(object):
    def sensitive(sens):
        ''' retrieves sensitive variables from the 'sensitive' directory in the 
            base dir. Useful for when you want to post code online and not have
            people trash your database... '''
        with open(path.join("/heap/cobDev/pyCOB/sensitive",sens),'r') as sens_file:
            return sens_file.readline().strip()
    # The Database Connection needs to be shared between all classes which inherit COBDatabase
    # Otherwise we open too many connections
    db = MySQLdb.connect(
        host   = sensitive("MySQLHost"),
        user   = sensitive("MySQLUser"),
        passwd = sensitive("MySQLPasswd"),
        db     = sensitive("MySQLDB")
    )

    def __init__(self,network="Developmental"):
        self.network = network
    def query(self, query, *variables):
        ''' Returns a pandas dataframe '''
        query = query.format(*variables)
        return psql.frame_query(query,self.db)

    def add_gene(self,gene_name):
        '''adds the gene name to the database '''
        cur = self.db.cursor()
        cur.execute("INSERT IGNORE INTO mzn_gene_name (gene_name) VALUES ('{}')".format(gene_name))

class Chrom(object):
    def __init__(self,id,length):
        self.id = id
        self.length = length
    
    def rQTL(self,length):
        ''' returns a random QTL within the chromosome '''
        start = np.random.randint(0,self.length)
        while start+length > self.length:
            start = np.random.randint(0,self.length)
        return QTL(chrom=self.id,start=start,end=start+length,id="RQTL-{}".format(length))
    def rLocus(self):
        ''' returns a random Locus from within the chromosome '''
        pos = np.random.randint(0,self.length)
        return Locus(chrom=self.id,start=pos,end=pos,id='rLocus-chr{}:{}'.format(self.id,pos))

class Genome(object):   
    def __init__(self,id,chroms=list()):
        self.id = id
        self.chroms = chroms 
    def add_chromosone(self,chrom):
        ''' adds a chromosome to itself '''
        self.chroms.append(chrom)
    def rChrom(self):
        ''' returns a random chromosome object from genome '''
        rindex = np.random.randint(0,len(self.chroms))
        return self.chroms[rindex]
    def rQTL(self,length):
        ''' returns a random QTL of specified lengths '''
        return self.rChrom().rQTL(length)
    def rLocus(self):
        ''' returns a random 'SNP' from the genome '''
        return self.rChrom().rLocus()
        
class Locus(COBDatabase):
    def __init__(self,chrom,start,end,id=None,gene_build='4a.53'):
        super(Locus,self).__init__()
        # Make sure we arent going nutso
        if start > end:
            oldstart = start
            start = end
            end = oldstart
        self.chrom = chrom
        self.start = start
        self.end   = end
        self.id    = id
        self.gene_build = gene_build

    def genes(self):
        ''' returns genes within the locus '''
        genes = self.query('''SELECT * FROM mzn_gene_loci WHERE chromosome = {}
            AND chromo_start > {}
            AND chromo_end < {}
            AND gene_build = '{}'
            ORDER BY chromo_start''', self.chrom, self.start, self.end, self.gene_build)
        return genes
    def upstream_genes(self,limit=100):
        ''' returns the X amount of upstream genes '''
        genes = self.query('''SELECT * from mzn_gene_loci WHERE chromosome = {}
            AND chromo_end < {}
            AND gene_build = '{}'
            ORDER BY chromo_start DESC
            LIMIT {}''',self.chrom,self.start,self.gene_build,int(limit))      
        return genes   
    
    def downstream_genes(self,limit=100):
        ''' returns the X amount of upstream genes '''
        genes = self.query(''' SELECT * FROM mzn_gene_loci WHERE chromosome = {}
            AND chromo_start > {} 
            AND gene_build = '{}'
            ORDER BY chromo_start
            LIMIT {}''',self.chrom,self.end,self.gene_build,int(limit))
        return genes
    def flanking_genes(self,limit=100):
        return pd.concat([
            self.upstream_genes(limit=np.ceil(limit/2)),
            self.downstream_genes(limit=np.ceil(limit/2))
        ])

    def __contains__(self,locus):
        if (locus.chrom == self.chrom and
               (( locus.start >= self.start and locus.start <= self.end)
               or(locus.end   <= self.end   and locus.end   >= self.start)
            )):
            return True
        else:
            return False
    def __len__(self):
        return self.end - self.start

    def __str__(self):
        return self.id
    def __repr__(self):
        return self.id


class SNP(Locus):
    def __init__(self,chrom,pos,alt=None,ref=None,id=None):
        super(SNP,self).__init__(chrom,pos,pos,id)
        self.chrom = chrom
        self.pos = pos
        self.id = id
        self.alt = alt
        self.ref = ref

    def nearby_genes(self,window=100000):
        downstream = np.floor(self.pos+(window/2))
        upstream   = np.floor(self.pos-(window/2))
        genes = self.query('''SELECT * FROM mzn_gene_loci WHERE chromosome = {}
            AND chromo_start > {}
            AND chromo_end < {}
            AND gene_build = '{}'
            ORDER BY chromo_start''', self.chrom, upstream, downstream, self.gene_build)
        return genes

    def nearest_gene(self,window=100000):
        genes = self.query('''SELECT * FROM mzn_gene_loci WHERE chromosome = {}
            AND chromo_start < {}
            AND chromo_end > {}
            AND gene_build = '{}'
            ORDER BY chromo_start''', self.chrom, self.pos, self.pos, self.gene_build)
        if genes.empty:
            downstream = np.floor(self.pos+(window/2))
            upstream   = np.floor(self.pos-(window/2))
            genes = self.query('''SELECT * FROM mzn_gene_loci WHERE chromosome = {}
                AND chromo_start > {}
                AND chromo_end < {}
                AND gene_build = '{}'
                ORDER BY chromo_start''', self.chrom, upstream, downstream, self.gene_build)
            if len(genes) > 2:
                if sum(genes.chromo_end < self.pos) == len(genes.chromo_end):
                    # All the genes are before the SNP, return the last
                    return genes.tail(1)
                elif sum(genes.chromo_start > self.pos) == len(genes.chromo_start):
                    # All genes are before SNP, return the first
                    return genes.head(1)
                else:
                    downstream_index = ((genes.chromo_start > self.pos) * 1).idxmax()
                    upstream_index = downstream_index -1
                    return genes.ix[[upstream_index,downstream_index]]


class QTL(Locus):
    def __init__(self,chrom,start,end,id=None):
        if id == None:
            self.id = "QTL-chr{}:{}-{}".format(chrom,start,end)
        super(QTL,self).__init__(chrom,start,end,self.id)
    def __str__(self):
        return "{} - {} - {} - {}".format(self.id,self.chrom,self.start,self.end)

class COB(COBDatabase):
    def __init__(self,network="Developmental"):
        super(COB,self).__init__()
        # HouseKeeping
        self.dataset_info = self.query("SELECT * FROM datasets WHERE name = '{}'",network)
        self.network   = network
        self.sig_edges = self.dataset_info.sig_edges[0]
        self.raw       = self.dataset_info.raw_edges[0]
        self.id        = self.dataset_info.id[0]
        self.log_file  = sys.stderr

    def id2gene(self,id):
        ''' returns a gene name from its id '''
        if isinstance(id,(str)):
            return id
        ser = self.query("SELECT gene_name FROM mzn_gene_name WHERE gene_id = {}",id)
        return(ser['gene_name'][0])

    def ids2genes(self,id_list):
        ''' returns a list of gene names from a list of ids '''
        return [self.id2gene(id) for id in id_list]

    def gene2id(self,gene,add_if_missing=False): 
        ''' returns an id from a gene name, option to add gene name to database '''
        if isinstance(gene,(int,long)):
                return gene
        ser = self.query("SELECT gene_id as id FROM mzn_gene_name WHERE gene_name = '{}'",gene)
        if ser.empty:
            ser = self.query("SELECT common_id as id FROM mzn_gene_common WHERE common_name = '{}'",gene)
        if ser.empty: 
            if add_if_missing == True:
                pass
                #self.query("INSERT INTO mzn_gene_name (gene_name) VALUES ({})",name) 
                #cur = self.query("SELECT gene_id FROM mnz_gene_name WHERE gene_name = '{}'",gene)
            else:
                return None
        return int(ser['id'])
    
    def genes2ids(self,gene_list):
        return [self.gene2id(gene) for gene in gene_list]

    def in_network(self,gene_list):
        gene_ids = self.genes2ids(gene_list)
        raw_exp = self.query('''SELECT * FROM mzn_rawexp_values 
            WHERE gene_id in ({})
            AND dataset_id = {}
        ''',','.join(map(str,gene_ids)),self.id)
        return raw_exp.gene_id.unique()

    def neighbors(self,gene_list):
        ''' Returns the neighbors for a list of genes '''
        gene_list = self.genes2ids(gene_list)
        return pd.concat([
            self.query('''SELECT gene_b as gene, gene_name as neighbor_name, gene_a as neighbor, score
            FROM {} net JOIN mzn_gene_name name ON net.gene_a = name.gene_id
            WHERE gene_b IN ({})''', self.sig_edges, ",".join(map(str,gene_list))),
            self.query('''SELECT gene_a as gene, gene_name as neighbor_name, gene_b as neighbor, score
            FROM {} net JOIN mzn_gene_name name ON net.gene_b = name.gene_id
            WHERE gene_a IN ({})''', self.sig_edges, ",".join(map(str,gene_list)))
        ])

    def num_neighbors(self,gene_list):
        ''' return the number of significant global and local neighbors '''
        gene_list = self.genes2ids(gene_list)
        neighbors = self.neighbors(gene_list)
        return pd.DataFrame({
            'gene'  : [gene for gene,group in neighbors.groupby('gene')],
            'global': [len(group) for gene,group in neighbors.groupby('gene')],
            'local' : [len(set(gene_list).intersection(group['neighbor'])) for gene,group in neighbors.groupby('gene')]
        })

    def density(self,gene_list):
        ''' calculate the density of the network between a set of genes '''
        if len(gene_list) == 0:
            return 0
        gene_list = self.genes2ids(gene_list)    
        scores = self.query("SELECT score FROM {} WHERE gene_a IN ({}) AND gene_b IN ({})",
            self.raw,
            ",".join(map(str,gene_list)),
            ",".join(map(str,gene_list))
        ).score
        if len(scores) == 0:
            return 0
        return (scores.mean()/(1/np.sqrt(len(scores))))
   
    def neighbors_score(self,gene_list):
        ''' returns a series containing the strongest connected neighbors '''
        neighbors = self.neighbors(gene_list)
        if len(neighbors) > 0:
            scores = neighbors.groupby('neighbor').score.sum()
            scores.sort(ascending=False)
        return scores

    def seed(self, gene_list, max_show = 65): 
        ''' given a set of nodes, add on the next X strongest connected nodes ''' 
        if len(gene_list) == 0:
            return []
        gene_list = self.genes2ids(gene_list)
        neighbors = self.neighbors_score(gene_list)
        seed_set =  set(gene_list).union(neighbors.index[0:min(len(neighbors),max_show)])
        return seed_set

    def subnetwork(self,gene_list):
        ''' Calculates a subnetwork based on connected nodes within the input gene list '''
        if len(gene_list) == 0:
            return []
        self.log("Analyzing {} Genes for {} Network",len(gene_list),self.network)
        self.log("Found {} genes in COB",len(gene_list))
        neighbors = self.num_neighbors(gene_list)
        local = neighbors[neighbors.local >= 1]
        self.log("Found {} genes in subnetwork",len(local))
        return(self.ids2genes(local.gene))
   
 
    def graph(self,gene_list):
        ''' Finds the largest connected component from the gene list '''
        if len(gene_list) == 0:
            return []
        neighbors = self.num_neighbors(gene_list)
        neighbors = neighbors[neighbors.local >= 1]
        local = self.query('''SELECT gene_a, gene_b, score FROM {}
            WHERE gene_a in ({}) and gene_b in ({})''',
            self.sig_edges,
            ",".join(map(str,neighbors.gene)),
            ",".join(map(str,neighbors.gene))
        )
        nodes = list(set(local.gene_a).union(local.gene_b))
        nodes.sort()
        graph = ig.Graph()
        for node in nodes:
            graph.add_vertex(name=node,label=self.id2gene(node))
        for edge in local[['gene_a','gene_b','score']].itertuples():
            graph.add_edge(source=graph.vs.find(name=edge[1]),target=graph.vs.find(name=edge[2]),weight=edge[3])
        return graph 
        
    def gene_expr_vals(self,gene_list):
        ids = self.genes2ids(gene_list)
        expr_vals = self.query(
         '''SELECT  val.gene_id, name, tissue_desc, organ_name, growth_stage , log_value, avg_log_value
            FROM mzn_rawexp_values val
            JOIN mzn_rawexp_avg avg ON val.gene_id = avg.gene_id AND val.dataset_id = avg.dataset_id 
            JOIN mzn_rawexp_accessions acc ON val.accession_id = acc.accession_id 
            WHERE val.gene_id IN ({})
            AND val.dataset_id = {}''',
            ",".join(map(str,ids)), 
            self.id
        )
        def zscore(group):
            group.index = group.name
            norm = group.log_value-group.avg_log_value
            return norm/np.std(norm)
        return pd.DataFrame({gid:zscore(group) for gid,group in expr_vals.groupby('gene_id')})

    def heatmap(self,dm):
       
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
        axmatrix.set_xticklabels(self.ids2genes(dm.columns[idx2]))
        axmatrix.xaxis.tick_bottom()
        axmatrix.tick_params(axis='x',labelsize='xx-small')
        axmatrix.yaxis.tick_right()
        axmatrix.set_yticks(np.arange(D.shape[0]))
        axmatrix.set_yticklabels(dm.index[idx1])
        axmatrix.tick_params(axis='y',labelsize='x-small')

        axColorBar = f.add_axes([0.09,0.75,0.2,0.05])
        f.colorbar(im,orientation='horizontal',cax=axColorBar,ticks=np.arange(np.ceil(vmin),np.ceil(vmax),2))

        return {'ordered' : D, 'rorder' : Z1['leaves'], 'corder' : Z2['leaves']} 

    def url(self,gene_list,size=65):
        gene_list = self.ids2genes(gene_list)
        return "http://lovelace.cs.umn.edu/cob?results_source={}&action=seed&q_list={}&neighborhood_size={}".format(
            self.network,
            ",".join(map(str,gene_list)),
            size
        )   
           
 

    ############################################################################
    # Helper Functions
    ############################################################################
    def log(self, string, *args, **kwargs):
        print('[LOG]',time.ctime(), "-", string.format(*args),file=self.log_file)


    @property 
    def cmap(self):
        heatmapdict = {'red': ((0.0, 1.0, 1.0),
                    (0.5, 1.0, 1.0),
                    (1.0, 0.0, 0.0)),
            'green': ((0.0, 1.0, 1.0),
                    (0.5, 1.0, 1.0),
                    (1.0, 0.0, 0.0)),
            'blue': ((0.0, 0.0, 0.0),
                    (0.5, 1.0, 1.0),
                    (1.0, 1.0, 1.0))}
        heatmap_cmap = plt.colors.LinearSegmentedColormap('my_colormap',heatmapdict,256)
        return heatmap_cmap


