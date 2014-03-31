#!/usr/bin/python

from __future__ import print_function

import matplotlib as plt
plt.use('Agg')
import pandas as pd
import pandas.io.sql as psql
import numpy as np
import sys
import os.path as path
import MySQLdb
import time
import igraph as ig
import time

from scipy.stats import hypergeom
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform

class Cobject(object):
    log_file = sys.stderr
    def init(self,log_file=None):
        if log_file:
            self.log_file = open(log_file,'w')
    def log(self,*args):
        ''' shared object logging '''
        print(time.ctime(),'-',*args,file=self.log_file)

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
    def execute(self,query,*variables):
        ''' perform a database execution '''
        cur = self.db.cursor()
        cur.execute(query.format(*variables))
    def add_gene(self,gene_name):
        '''adds the gene name to the database '''
        cur = self.db.cursor()
        cur.execute("INSERT IGNORE INTO mzn_gene_name (gene_name) VALUES ('{}')".format(gene_name))
    def gene_id(self,gene_name):
        pass
    

class COBGene(COBDatabase):
    def __init__(self,id):
        self.id = id   

    @property
    def gramene(self):
        pass
 

class COBDataset(COBDatabase):
    def __init__(self,name,description):
        super(COBDataset,self).__init__()
        self.name = name
        self.description = description
        self.exp_vals = pd.DataFrame()
        self.id = None

    def save(self):
        self.id = self.query("SELECT MAX(id) as MID FROM datasets;").iloc[0]['MID'] + 1
        # add the entry into the database
        self.execute("INSERT INTO datasets (id, name, description) VALUES ({}, '{}', '{}')",self.id,self.name,self.description)

    
    

def from_csv(filename,FPKM=True,sep="\t"):
    expr_vals = pd.read_table(filename,sep=sep)
    try:
        expr_vals[expr_vals.columns] = expr_vals[expr_vals.columns].convert_objects(convert_numeric = True)
    except e:
        exit("csv expression values must be numbers")
    

    # Register new dataset
    # Need to import rawexp values
        # register new accessions
        # register new genes/probes
    # Import new accessions
    # 
    # Need to import 


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

    def genes(self,gene_build=None):
        ''' returns genes within the locus '''
        if not gene_build:
           gene_build = self.gene_build 
        genes = self.query('''SELECT * FROM mzn_gene_loci WHERE chromosome = {}
            AND chromo_start > {}
            AND chromo_end < {}
            AND gene_build = '{}'
            ORDER BY chromo_start''', self.chrom, self.start, self.end, gene_build)
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
        return "Locus Type Object: {}".format(str(self.id))
    def __repr__(self):
        return str(self.id)


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
        else:
            self.id = id
        super(QTL,self).__init__(chrom,start,end,self.id)
    def __str__(self):
        return "{}".format(self.id)


class Gene(Locus):
    def __init__(self,id,gene_build='4a.53'):
        info = self.query('''SELECT chromosome, chromo_start, chromo_end FROM mzn_gene_loci 
            WHERE gene_id = {} and gene_build = '{}' ''',id,gene_build)
        super(Gene,self).__init__(info.iloc[0]['chromosome'],info.iloc[0]['chromo_start'],info.iloc[0]['chromo_end'])
        self.id = int(id)

    @property
    def gramene_id(self):
        return self.query("SELECT gene_name FROM mzn_gene_name WHERE gene_id = {}",self.id).iloc[0]['gene_name']

    @property
    def common_name(self):
        info = self.query("SELECT common_name FROM mzn_gene_common WHERE common_id = {}",self.id)
        return None if info.shape[0] == 0 else info.iloc[0]['common_name']

    @property
    def arab_ortho(self):
        info = self.query('''SELECT * FROM mzn_arab_orthologs orth
            LEFT JOIN mzn_arab_gene info ON info.arab_id = orth.arab_id 
            LEFT JOIN mzn_arab_gene_types type ON info.type_id = type.type_id
            LEFT JOIN mzn_arab_short_desc short ON info.short_id = short.short_id
            LEFT JOIN mzn_arab_curator_desc cur ON info.curated_id = cur.curator_id
            LEFT JOIN mzn_arab_comp_desc comp ON info.comp_id = comp.comp_id
            WHERE gene_id = {} ''',self.id)
        return info

    @property
    def go_terms(self):
        info = self.query('''SELECT term_name, term_short_desc, term_long_desc, space_desc
            FROM mzn_gene_go_terms
            LEFT JOIN mzn_go_terms ON go_id = term_id
            LEFT JOIN mzn_go_space ON term_space = space_id
            WHERE gene_id = {}
        ''',self.id)
        return info


            

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
        self.verbose   = True
        self.num_go_genes_total = 25288

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
        try:
            id = int(gene)
            return id
        except:
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
        axmatrix.set_xticklabels(self.ids2genes(dm.columns[idx2]))
        axmatrix.xaxis.tick_bottom()
        axmatrix.tick_params(axis='x',labelsize='xx-small')
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


