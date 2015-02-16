#!/usr/bin/python3


from camoco.Camoco import Camoco
from camoco.RefGen import RefGen
from camoco.Locus import Locus,Gene
from camoco.Expr import Expr
from camoco.Tools import memoize

from numpy import matrix,arcsinh
from collections import defaultdict
from itertools import chain
from subprocess import Popen, PIPE

import pandas as pd
import igraph as ig
import numpy as np
import time
import math as math
import multiprocessing as multiprocessing
import itertools
import matplotlib.pylab as plt
import matplotlib
import io
import base64
import os

from scipy.stats import hypergeom,pearsonr

import collections

class COB(Expr):
    def __init__(self,name=None,description=None,basedir="~/.camoco"):
        if name is None:
            self.log('You must provide a name')
        else:
            super().__init__(name=name,description=description,basedir=basedir)
            self._create_tables()

    #@memoize
    def num_edges(self,sig_only=True):
        clause = ' WHERE significant =1;' if sig_only else ''
        try:
            return self.db.cursor().execute("SELECT COUNT(*) FROM coex"+clause+';').fetchone()[0]
        except ValueError as e:
            self.log('No edges in database!')

 
    def edges(self,min_distance=None,sig_only=False):
        ''' returns edges '''
        try:
            query = '''SELECT gene_a,gene_b,score,distance
                    FROM coex '''
            if min_distance:
                query += ' WHERE distance < {}'.format(min_distance)
            if sig_only:
                query += ' WHERE significant = 1' if min_distance is None else 'AND significant = 1'
            return pd.DataFrame(self.db.cursor().execute(query).fetchall(),
                columns=['gene_a','gene_b','score','distance'])
        except ValueError as e:
            self.log('Oops! there are no signficant edges in dataset...')

    #@memoize
    def sig_genes(self):
        sig_edges = self.edges(sig_only=True)
        return list(self.refgen.from_ids(set(sig_edges.gene_a).union(sig_edges.gene_b)))

    def coexpression(self,gene_a,gene_b,trans_only=True):
        ''' returns coexpression z-score between two genes '''
        return self.db.cursor().execute(
            ''' SELECT score FROM coex WHERE 
                (gene_a = ? AND gene_b = ?) 
                OR (gene_b = ? AND gene_a = ?)'''
            ,(gene_a.id,gene_b.id,gene_a.id,gene_b.id)
        ).fetchone()[0]
    
    def global_degree(self,gene):
        try:
            return self.db.cursor().execute(''' 
                SELECT degree FROM degree WHERE gene = ?
            ''',(gene.id,)).fetchone()[0]
        except TypeError as e:
            return -1

            return len(self.neighbors([gene]))
    def neighbors(self,gene_list,min_distance=None,sig_only=True):
        ''' Input : a list of COBGene Objects
            Output: Returns the neighbors for a list of genes as a DataFrame
            Columns returned: query_name, queryID, neighbor_name, neighbor, score
        '''
        # Ids need to be extracted from gene list
        gene_list = "','".join([x.id for x in gene_list])
        cur = self.db.cursor()
        query = ''' 
            SELECT {} AS source, {} AS target, score, distance FROM coex 
            WHERE 
            {} IN ('{}')  
            ''' 
        # Tag on qualifiers for filtering 
        if min_distance:
            query += ' AND distance > {} '.format(min_distance)
        if sig_only:
            query += ' AND significant = 1'
        query += ';'
        # The way we do the queries makes sure that the genes in the gene_list are always
        # in the first column. Neighbors are always in the second column. 
        data=(cur.execute(query.format('gene_a','gene_b','gene_a',gene_list)).fetchall() +
              cur.execute(query.format('gene_b','gene_a','gene_b',gene_list)).fetchall())
        if len(data) == 0:
            return pd.DataFrame(columns=['source','target','score','distance'])
        else:
            return pd.DataFrame(data=data,
                columns=['source','target','score','distance'] 
            ) 

    def degree(self,gene_list,min_distance=None):
        ''' 
            Input: a list of Gene Objects
            Output: Dataframe containing the number of significant global
                    and local neighbors for each gene in list
        '''
        neighbors = self.neighbors(gene_list,min_distance=min_distance,sig_only=True)
        if len(neighbors) == 0:
            return pd.DataFrame(columns=["global",'local'])
        def get_nums(df):
            global_length = len(df)
            local_length = len(set(neighbors['source']).intersection(set(df['target'])))
            return pd.Series({
                'global' : global_length,
                'local'  : local_length
            })
        return neighbors.groupby('source').apply(get_nums)

    def next_neighbors(self,gene_list):
        ''' returns a list containing the strongest connected neighbors '''
        neighbors = self.neighbors(gene_list)
        if len(neighbors) > 0:
            scores = neighbors.groupby('target').score.sum()
            scores.sort(ascending=False)
            return scores
        else:
            return pd.Series()

    def neighborhood(self,gene_list):
        ''' Input: A gene List
            Output: a Dataframe containing gene ids which have at least one edge
                    with another gene in the input list. Also returns global degree '''
        if len(gene_list) == 0:
            return []
        self.log("Analyzing {} Genes for {} Network",len(gene_list),self.name)
        self.log("Found {} genes in COB",len(gene_list))
        gene_degree = self.degree(gene_list)
        local = gene_degree[gene_degree.local >= 1]
        self.log("Found {} genes in subnetwork",len(local))
        return(local)

    def subnetwork(self,gene_list,sig_only=True,min_distance=100000):
        ''' Input: a gene list
            Output: a dataframe containing all edges EXCLUSIVELY between genes within list
        '''
        if len(gene_list) == 0:
            return pd.DataFrame()
        # filter for only genes within network
        #gene_list = list(filter(lambda x: x in self, gene_list))
        ids = "','".join([x.id for x in gene_list])
        query = '''
            SELECT gene_a,gene_b,score,distance FROM coex 
            WHERE gene_a IN ('{}') AND gene_b IN ('{}')
            '''
        if min_distance:
            query += ' AND distance > {} '.format(min_distance)
        if sig_only:
            query += ' AND significant = 1;'
        else:
            query += ';'
        data = self.db.cursor().execute(query.format(ids,ids)).fetchall()
        if len(data) == 0:
            return pd.DataFrame(columns=['source','target','score','distance'])
        else:
            return pd.DataFrame(
                data = data,
                columns = ['source','target','score','distance']
            )   
        return subnet

    def lcc(self,gene_list,min_distance=None):
        ''' returns an igraph of the largest connected component in graph '''
        try:
            return self.graph(gene_list).clusters().giant()
        except ValueError as e:
            return ig.Graph() 

    def locality(self,gene_list,min_distance=None):
        ''' Returns the log ratio of median local degree to meadian global degree. 
            This gives you an idea of a group of genes affinity towards one another,
            or if you just happened to stumble into a high degree set.'''
        gene_list = list(filter(lambda x: x in self, gene_list))
        if len(gene_list) == 0:
            return np.nan
        gene_degrees = self.degree(gene_list)
        if sum(gene_degrees['global']) == 0:
            return np.nan
        return(sum(gene_degrees['local'])/sum(gene_degrees['global']))

    def density(self,gene_list,return_mean=True,min_distance=50000):
        ''' calculates the denisty of the non-thresholded network amongst genes
            not within a certain distance of each other. This corrects for
            cis regulatory elements increasing noise in coexpression network '''
        # filter for only genes within network
        edges = self.subnetwork(gene_list,min_distance=min_distance,sig_only=False)
        if len(edges) == 0:
            return np.nan
        if len(edges) == 1:
            return edges.score[0]
        if return_mean:
            return (edges.score.mean()/((edges.score.std())/np.sqrt(len(edges))))
        else:
            return edges

    def seed(self, gene_list, limit=65): 
        ''' Input: given a set of nodes, add on the next X strongest connected nodes ''' 
        if len(gene_list) == 0:
            return []
        neighbors = self.next_neighbors(gene_list)
        seed_set =  list(self.refgen.from_ids(
            set([x.id for x in gene_list]).union(neighbors.index[0:min(len(neighbors),limit)])
        ))
        return seed_set

    def graph(self,gene_list,min_distance=None):
        ''' Input: a gene list
            Output: a iGraph object '''
        if len(gene_list) == 0:
            return ig.Graph()
        # retrieve all first neighbors
        edges = self.subnetwork(gene_list,min_distance=min_distance)
        nodes = list(set(edges.source).union(edges.target))
        idmap = dict()
        for i,node in enumerate(nodes):
            idmap[node] = i
        graph = ig.Graph(
            edges = [(idmap[source],idmap[target]) for source,target 
                in edges[['source','target']].itertuples(index=False)],
            directed = False,
            vertex_attrs = {
                "label" : nodes,
            },
            edge_attrs = {
                "weight" : edges.score
            }
        )
        return graph 


    def squaredist(self,gene_list=None):
        ''' Returns a gene x gene square distance dataframe ''' 
        if gene_list is None:
            gene_list = self.sig_genes()
        subnet = self.subnetwork(gene_list,sig_only=False,min_distance=None).sort(['source','target'])
        # There are at least one gene in each column who doesnt have a corresponding entry in the 
        # 
        lone_source = set(subnet.source).difference(subnet.target).pop()
        lone_target = set(subnet.target).difference(subnet.source).pop()
        lone_mask = (subnet.source == lone_source)&(subnet.target == lone_target)
        subnet.ix[lone_mask,'source'] = lone_target
        subnet.ix[lone_mask,'target'] = lone_source
        dat = pd.pivot(subnet.source,subnet.target,subnet.score) 
        dat[dat.isnull()] = 0
        return dat + dat.T

    def to_treeview(self, gene_list=None):
        ''' outputs treeview files for dataset '''
        if gene_list is None:
            gene_list = self.sig_genes()
        gene_expr_vals = self.expr(gene_list,zscore=True) 
        
    def pcc_heatmap(self,gene_list=None):
        pass
        
    def coordinates(self,gene_list,layout=None):
        ''' returns the static layout, you can change the stored layout by passing 
            in a new layout object. If no layout has been stored or a gene does not have
            coordinates, returns (0,0) for each mystery gene'''
        if not layout:
            # return the current layout for gene list
            return ig.Layout([
                self.db.cursor().execute(
                    '''SELECT x,y FROM coor WHERE gene = ?''',
                    (gene.id,)).fetchone() or (0,0) for gene in gene_list
            ])
        else:
            self.db.cursor().executemany('''
                INSERT OR REPLACE INTO coor (gene,x,y) VALUES (?,?,?)''',
                [(gene.id,layout[i][0],layout[i][1]) for i,gene in enumerate(gene_list)]
            )

    def mcl(self,gene_list,I=2.0,scheme=7,min_distance=100000,min_size=0,max_size=10e10):
        ''' A *very* thin wrapper to the MCL program. The MCL program must
            be accessible by a subprocess (i.e. by the shell).
            Returns clusters (as list) as designated by MCL. 
            Input: a gene list
            Output: a list of lists of genes within each cluster
        '''
        MCLstr = "\n".join(["{}\t{}\t{}".format(a,b,c) for a,b,c in \
            self.subnetwork(gene_list,min_distance=min_distance)[['source','target','score']].itertuples(index=False)])
        cmd = "mcl - --abc -scheme {} -I {} -o -".format(scheme,I)
        self.log("running MCL: {}",cmd)
        try:
            p = Popen(cmd, stdout=PIPE, stdin=PIPE, stderr=self.log_file, shell=True)
            sout,serr = p.communicate(MCLstr.encode('utf-8'))
            p.wait()
            if p.returncode==0 and len(sout)>0:
                # Filter out cluters who are smaller than the min size
                return list(filter(lambda x: len(x) > min_size and len(x) < max_size,
                    # Generate ids from the refgen
                    [ self.refgen.from_ids([gene.decode('utf-8') for gene in line.split()]) for line in sout.splitlines() ]
                ))
            else:
                self.log( "MCL return code: {}".format(p.returncode))
        except FileNotFoundError as e:
            self.log('Could not find MCL in PATH. Make sure its installed and shell accessible as "mcl".')

    def cluster(self,raw=False):
        # Convert to tsv
        rawtype = 'raw' if raw else 'norm'
        filename  = '{}_{}.tsv'.format(self.name, rawtype)
        self.expr(raw=raw,zscore=True).to_csv(filename,sep='\t')
        try: 
            cmd = ['cluster', '-f', filename, '-g', '2', '-e', '2']
            self.log("Executing {}",' '.join(cmd))
            p = Popen(cmd, stderr=self.log_file, shell=True)
            self.log('Waiting for {} cluster...'.format(filename))
            p.wait()
        except FileNotFoundError as e:
            self.log('Could not find cluster command in PATH. Make sure its installed and shell accessible as "cluster".')


    def to_dat(self,gene_list,filename):
        with open(filename, 'w') as OUT:
            for a,b,c in self.subnetwork(gene_list).itertuples(index=False):
                print("{}\t{}\t{}".format(a,b,c),file=OUT)

    def plot_pcc_hist(self,filename=None):
        ''' Plot the histogram of PCCs. Right now this is HUGELY inefficient in memory.'''
        if filename is None:
            filename = self.name+'.png'
        plt.clf()
        plt.hist(np.tanh(self.edges(sig_only=False,min_distance=None).score),bins=50)
        plt.xlabel('PCC')
        plt.ylabel('Freq')
        plt.savefig(filename) 
        
    def plot(self,gene_list,filename=None,width=3000,height=3000,layout=None,**kwargs):
        if not filename:
            filename = "{}.png".format(self.name)
        if not layout:
            layout = self.coordinates(gene_list)
        ig.plot(self.graph(gene_list),layout=layout,target=filename,bbox=(width,height),**kwargs)

    def compare_to_COB(self,list_of_other_COBs,filename=None,gridsize=100,extent=[-10,10,-10,10]):
        ''' Compare the edge weights in this COB to another COB. Prints out edge weights to file'''
        for other_COB in list_of_other_COBs:
            self.log("Comparing {} to {}",self.name,other_COB.name)
            filename = "{}_to_{}".format(self.name,other_COB.name)
            # Print out the edge comparisons for each common gene
            self.log("Printing out common gene edges")
            if not os.path.exists(filename+'.tsv'):
                with open(filename+'.tsv','w') as OUT:
                    # Only compare the common genes
                    for i,common_gene in enumerate(set(self.genes()).intersection(other_COB.genes())):
                        x = self.neighbors([common_gene],sig_only=False)
                        y = other_COB.neighbors([common_gene],sig_only=False)
                        # Merge on the target column and print
                        x[['target','score']].merge(
                            y[['target','score']], on='target'
                        ).to_csv(OUT, sep="\t", header=False, index=None)
                        if i % 1000 == 0:
                            self.log("Done Processing {} genes",i)
            if not os.path.exists(filename+'.png'):
                # Read in the table
                self.log("Reading in {}",filename+'.tsv')
                tbl = pd.read_table(filename+".tsv",names=['target','x','y'])
                from matplotlib import cm
                self.log('Generating hexbin')
                # http://matplotlib.org/examples/pylab_examples/hexbin_demo.html
                plt.clf()
                plt.figure(figsize=(8,8),dpi=200)
                plt.hexbin(tbl['x'],tbl['y'],gridsize=gridsize,cmap=cm.afmhot,bins='log',extent=extent)
                plt.xlabel('{} Edge Z Score'.format(self.name))
                plt.ylabel('{} Edge Z Score'.format(other_COB.name))
                plt.savefig(filename+'.png')
                cb = plt.colorbar()
                cb.set_label('log edge Z-score counts')
            self.log("Done")

    def compare_to_dat(self,filename,sep="\t",score_cutoff=3):
        ''' Compare the number of genes with significant edges as well as degree with a DAT file '''
        self.log("Reading in {}",filename)
        ref_dat = pd.read_table(filename,sep=sep,names=['gene_a','gene_b','score'])
        self.log("Retrieving network edges...")
        dat = pd.DataFrame(self.db.cursor().execute(
            '''SELECT gene_a, gene_b, score 
            FROM  coex 
            WHERE score >= ?;
            ''',(score_cutoff,)).fetchall(),
            columns=['gene_a','gene_b','score']
        )
        ref_set = set(list(ref_dat.gene_a)+list(ref_dat.gene_b))
        # Get a list of genes in the dataset as well as reference
        gene_freqs = defaultdict(lambda:0)
        for gene_a,gene_b in (dat[["gene_a",'gene_b']].itertuples(index=False)):
            if gene_a in ref_set and gene_b in ref_set:
                gene_freqs[gene_a] += 1
                gene_freqs[gene_b] += 1
        ref_gene_freqs = defaultdict(lambda:0)
        for gene_a,gene_b in (ref_dat[["gene_a",'gene_b']].itertuples(index=False)):
            if gene_a in ref_set and gene_b in ref_set:
                ref_gene_freqs[gene_a] += 1
                ref_gene_freqs[gene_b] += 1
        # report list sizes
        self.log("Found {} genes in {}",len(gene_freqs),self.name)
        self.log("Found {} genes in {}",len(ref_gene_freqs),filename)
        # calculate the sets of significant genes
        self.log("Merging degree sets...")
        degrees = pd.DataFrame(pd.Series(gene_freqs),columns=["Dat"]).merge(
            pd.DataFrame(pd.Series(ref_gene_freqs),columns=["RefDat"]),
            how='right',
            left_index=True,
            right_index=True
        )
        self.log("...Done")
        return degrees

    def compare_degree(self,obj,score_cutoff=3):
        self.log("Retrieving network edges...")
        ref_dat = pd.DataFrame(obj.db.cursor().execute(''' 
            SELECT gene_a,gene_b,score FROM coex
            WHERE score > ?
        ''',(score_cutoff,)).fetchall(),columns=['gene_a','gene_b','score'])
        dat = pd.DataFrame(self.db.cursor().execute(
            '''SELECT gene_a, gene_b, score 
            FROM  coex 
            WHERE score >= ?;
            ''',(score_cutoff,)).fetchall(),
            columns=['gene_a','gene_b','score']
        )
        ref_set = set(list(ref_dat.gene_a)+list(ref_dat.gene_b))
        # Get a list of genes in the dataset as well as reference
        gene_freqs = defaultdict(lambda:0)
        for gene_a,gene_b in (dat[["gene_a",'gene_b']].itertuples(index=False)):
            if gene_a in ref_set and gene_b in ref_set:
                gene_freqs[gene_a] += 1
                gene_freqs[gene_b] += 1
        ref_gene_freqs = defaultdict(lambda:0)
        for gene_a,gene_b in (ref_dat[["gene_a",'gene_b']].itertuples(index=False)):
            if gene_a in ref_set and gene_b in ref_set:
                ref_gene_freqs[gene_a] += 1
                ref_gene_freqs[gene_b] += 1
        # report list sizes
        self.log("Found {} genes in {}",len(gene_freqs),self.name)
        self.log("Found {} genes in {}",len(ref_gene_freqs),obj.name)
        # calculate the sets of significant genes
        self.log("Merging degree sets...")
        degrees = pd.DataFrame(pd.Series(gene_freqs),columns=["Dat"]).merge(
            pd.DataFrame(pd.Series(ref_gene_freqs),columns=["RefDat"]),
            how='right',
            left_index=True,
            right_index=True
        )
        self.log("...Done")
        return degrees

    def summary(self):
        print( '''
            COB Dataset: {} - {} - {}
                Desc: {}
                FPKM: {}
                Num Genes: {}
                Num Accessions: {}
                Sig./Total Interactions: {}
        '''.format(self.name, self.organism, self.build,
                self.description,
                self.FPKM == 'RNASEQ',
                self.num_genes(),
                self.num_accessions,
                self.num_edges(sig_only=True)/self.num_edges(sig_only=False)
        ))

    def __repr__(self):
        return '''
            COB Dataset: {} - {} - {}
                Desc: {}
                FPKM: {}
        '''.format(self.name, self.organism, self.build,
                self.description,
                self.FPKM == 'FPKM',
        )

    def __str__(self):
        return self.__repr__()

    ''' ------------------------------------------------------------------------------------------
            Class Methods
    '''
    def _calculate_coexpression(self,significance_thresh=3):
        cur = self.db.cursor() 
        try:
            self._drop_coex()
            cur.execute("BEGIN TRANSACTION")
            self._drop_indices()
            tbl = pd.DataFrame(
                # The sorted call here is KEY! It means the columns in the database will be guaranteed
                # to be sorted so we can access them faster
                list(itertools.combinations([x.id for x in sorted(self.genes())],2)),
                columns=['gene_a','gene_b']
            )
            # Now add coexpression data
            self.log("Calculating Coexpression")
            # NOTE we have to call this outside _coex function because it parallizes things
            tbl['score'] = self._coex(self.expr().as_matrix())
            # correlations of 1 dont transform well
            tbl['score'][tbl['score'] == 1] = 0.99999999
            # Perform fisher transform on PCCs
            tbl['score'] = np.arctanh(tbl['score'])
            # Sometimes, with certain datasets, the NaN mask overlap completely for the
            # two genes expression data making its PCC a nan. This affects the mean and std fro the gene.
            valid_scores = np.ma.masked_array(tbl['score'],np.isnan(tbl['score']))
            # Calculate Z Scores
            tbl['score'] = (valid_scores-valid_scores.mean())/valid_scores.std()
            tbl['significant'] = pd.Series(list(tbl['score'] >= significance_thresh),dtype='int_')
            self.log("Calculating Gene Distance")
            tbl['distance'] = self._gene_distances(self.genes())
            self.log("Building Database")
            cur.executemany(
                ''' INSERT INTO coex (gene_a,gene_b,score,significant,distance) VALUES (?,?,?,?,?)''',
                map(lambda x: (x[0],x[1],float(x[2]),int(x[3]),x[4]), tbl.itertuples(index=False))
            )
            cur.execute("END TRANSACTION")
            self.log("Building indices")
            self._build_indices()
            self.log("Done")
        except Exception as e:
            self.log("Something bad happened:{}",e)
            cur.execute("ROLLBACK")
            raise e
        # update the reference genome
        self._filter_refgen()  
        return self
 
    @classmethod
    def from_Expr(cls,expr):
        ''' Create a coexpression object from an Expression object '''
        # The Expr object already exists, just get a handle on it
        self = cls(name=expr.name,description=expr.description,basedir=expr.basedir)
        # Grab a coffee
        self._calculate_coexpression()
        self._build_degree()
        return self

    @classmethod
    def from_DataFrame(cls,df,name,description,refgen,rawtype=None,basedir='~/.camoco',**kwargs):
        # Create a
        expr = super().from_DataFrame(
            df,name,description,refgen,rawtype,basedir=basedir,**kwargs
        )
        return cls.from_Expr(expr)

    @classmethod
    def from_csv(cls,filename,name,description,refgen,rawtype=None,basedir='~/.camoco',sep='\t',**kwargs):
        ''' Build a COB Object from an FPKM or Micrarray CSV. '''
        return cls.from_DataFrame(pd.read_table(filename),name,description,refgen,rawtype=rawtype,basedir=basedir,**kwargs)
        
    def _coex(self,matrix,method='multiprocessing'):
        ''' 
            Coexpression abstraction method. If you want to try to make calculation faster,
                optimize this method. 
            Input: a numpy matrix
            Output: a vector of PCC correlations among all rows of matrix 
                    (matches itertools.combinations for all row labels)
        '''
        progress = progress_bar(multiprocessing.Manager())
        pool = multiprocessing.Pool()
        scores = np.array(list(itertools.chain(
            *pool.imap(_PCCUp,[(i,matrix,progress) for i in range(len(matrix))])
        )))
        pool.close()
        pool.join()
        return scores

    def _gene_distances(self,genes):
        progress = progress_bar(multiprocessing.Manager())
        pool = multiprocessing.Pool()
        distances = np.array(list(itertools.chain(
            *pool.imap(_DISUp,[(i,genes,progress) for i in range(len(genes))])
        )))
        pool.close()
        pool.join()
        return distances


    def _build_degree(self):
        degree = defaultdict(lambda : 0)
        for gene_a,gene_b in self.db.cursor().execute('''
                SELECT gene_a,gene_b FROM coex WHERE significant = 1;
            ''').fetchall():
            degree[gene_a] += 1
            degree[gene_b] += 1
        self.db.cursor().executemany(''' 
            INSERT INTO degree (gene,degree) VALUES (?,?);
        ''',degree.items())

    def _drop_coex(self):
        cur = self.db.cursor()
        cur.execute('BEGIN TRANSACTION')
        cur.execute(''' 
            DROP INDEX IF EXISTS coex_abs;
            DROP INDEX IF EXISTS coex_gene_b;
            DROP INDEX IF EXISTS coex_gen_ab;
            DROP INDEX IF EXISTS coex_score;
            DROP INDEX IF EXISTS coex_significant;
            DROP TABLE IF EXISTS coex;
            DROP TABLE IF EXISTS degree;
        ''')
        cur.execute('END TRANSACTION')
    def _build_indices(self):
        cur = self.db.cursor()
        cur.execute('BEGIN TRANSACTION')
        cur.execute(''' 
            CREATE INDEX IF NOT EXISTS coex_abs ON coex (gene_a); 
            CREATE INDEX IF NOT EXISTS coex_gene_b ON coex (gene_b); 
            CREATE INDEX IF NOT EXISTS coex_gen_ab ON coex (gene_a,gene_b); 
            CREATE INDEX IF NOT EXISTS coex_score ON coex (score); 
            CREATE INDEX IF NOT EXISTS coex_significant ON coex (significant); 
            CREATE INDEX IF NOT EXISTS coor_gene ON coor(gene);
            CREATE INDEX IF NOT EXISTS coor_x ON coor(x);
            CREATE INDEX IF NOT EXISTS coor_y ON coor(y);
            ANALYZE;
        ''')
        cur.execute('END TRANSACTION')

    def _drop_indices(self):
        self.db.cursor().execute(''' 
            DROP INDEX IF EXISTS coex_gene_a;
            DROP INDEX IF EXISTS coex_gene_b;
            DROP INDEX IF EXISTS coex_score;
            DROP INDEX IF EXISTS coex_significant;
        ''')

    def _create_tables(self):
        # build the expr tables too
        super()._create_tables()
        cur = self.db.cursor()
        cur.execute(''' 
            CREATE TABLE IF NOT EXISTS degree (
                gene TEXT,
                degree INTEGER
            );
            CREATE TABLE IF NOT EXISTS coex ( 
                gene_a TEXT,
                gene_b TEXT,
                score REAL,
                significant INTEGER,
                distance INTEGER
            );
            CREATE TABLE IF NOT EXISTS coor (
                gene TEXT,
                x REAL,
                y REAL
            );
        ''') 


class progress_bar(object):
    def __init__(self,manager,interval=10): 
        self.interval = interval
        self.lock = manager.Lock()
        self.old_per = manager.Value(float,0.0)
    def update(self,cur_per):
        cur_per = math.floor(cur_per)
        if cur_per > self.old_per.value and cur_per % self.interval == 0:
            self.lock.acquire()
            print("\t\t{} {}% Complete".format(time.ctime(),cur_per))
            self.old_per.set(cur_per)
            self.lock.release()


def _PCCUp(tpl):
    ''' Returns a tuple containing indices i,j and their Pearson Correlation Coef '''
    i,m,pb = (tpl) # values are index, exprMatrix, and lock
    total = ((m.shape[0]**2)-m.shape[0])/2 # total number of calculations we need to do
    left = ((m.shape[0]-i)**2-(m.shape[0]-i))/2 # How many we have left
    percent_done = (1-left/total)*100 # percent complete based on i and m
    pb.update(percent_done)
    vals = list()
    for j in range(i+1,len(m)): 
        mask = np.logical_and(np.isfinite(m[i,:]),np.isfinite(m[j,:]))
        if all(mask == False):
            vals.append(np.nan)
        else:
            vals.append(pearsonr(m[i,mask],m[j,mask])[0] )
    return vals

def _DISUp(tpl):
    ''' Returns a tuple containing the indices i,j and their distance '''
    i,genes,pb = (tpl) # values are index, gene_list, and lock
    total = ((len(genes)**2)-len(genes))/2 # total number of calculations we need to do
    left = ((len(genes)-i)**2-(len(genes)-i))/2 # How many we have left
    percent_done = (1-left/total)*100 # percent complete based on i and m
    pb.update(percent_done)
    vals = list()
    for j in range(i+1,len(genes)): 
        vals.append(abs(genes[i]-genes[j]))
    return vals

   
