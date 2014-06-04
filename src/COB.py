#!/usr/bin/python3

import pandas as pd

from Camoco import Camoco
from numpy import matrix,arcsinh
from Locus import Locus,Gene
from scipy.stats import pearsonr

import igraph as ig
import numpy as np
import pandas as pd
import time
import math as math
import multiprocessing as multiprocessing
import itertools
import matplotlib as plt

from scipy.stats import hypergeom
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform
 

class COB(Camoco):
    def __init__(self,name=None,basedir="~/.camoco"):
        if name is None:
            self.log('You must provide a name')
        else:
            super().__init__(name=name,type="COB",basedir=basedir)

    @property
    def num_genes(self):
        return self.db.cursor().execute("SELECT COUNT(DISTINCT(gene)) FROM expression").fetchone()[0]

    @property
    def num_accessions(self):
        return self.db.cursor().execute("SELECT COUNT(DISTINCT(accession)) FROM expression").fetchone()[0]
   
    @property
    def num_sig_interactions(self):
        return self.db.cursor().execute("SELECT COUNT(*) FROM coex WHERE significant = 1;").fetchone()[0]


    def neighbors(self,gene_list):
        ''' 
            Input : a list of COBGene Objects
            Output: Returns the neighbors for a list of genes as a DataFrame
            Columns returned: query_name, queryID, neighbor_name, neighbor, score
        '''
        gene_list = "','".join([x.id for x in gene_list])
        cur = self.db.cursor()
        return pd.DataFrame(
            data=(cur.execute(''' 
                SELECT gene_a AS source, gene_b AS target, score FROM coex 
                WHERE 
                gene_a IN ('{}')  
                AND significant = 1;
                '''.format(gene_list,gene_list)).fetchall() + 
                cur.execute(''' 
                SELECT gene_b AS source, gene_a AS target, score FROM coex 
                WHERE 
                gene_b IN ('{}')  
                AND significant = 1;
                '''.format(gene_list,gene_list)).fetchall()
            ),
            columns=['source','target','score'] 
        ) 

    def num_neighbors(self,gene_list):
        ''' 
            Input: a list of Gene Objects
            Output: Dataframe containing the number of significant global
                    and local neighbors for each gene in list
        '''
        neighbors = self.neighbors(gene_list)
        def get_nums(df):
            global_length = len(df)
            local_length = len(set(neighbors['source']).intersection(set(df['target'])))
            return pd.Series({
                'global' : global_length,
                'local'  : local_length
            })
        return neighbors.groupby('source').apply(get_nums)

    def degree(self,gene_list):
        return self.num_neighbors(gene_list)

    def neighbors_score(self,gene_list):
        ''' returns a series containing the strongest connected neighbors '''
        neighbors = self.neighbors(gene_list)
        if len(neighbors) > 0:
            scores = neighbors.groupby('target').score.sum()
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
            Output: a dataframe containing all edges EXCLUSIVELY between genes within list
        '''
        if len(gene_list) == 0:
            return pd.DataFrame()
        gene_list = "','".join([x.id for x in gene_list])
        subnet = pd.DataFrame(
            self.db.cursor().execute(''' 
            SELECT gene_a, gene_b, score FROM coex
            WHERE significant = 1 
            AND gene_a IN ('{}')
            AND gene_b IN ('{}') '''.format(gene_list,gene_list)).fetchall(),
            columns = ['source','target','score']
        )   
        return subnet

    def density(self,gene_list):
        ''' calculate the density of the ENTIRE network between a set of genes 
            Input: A list of COBGene Objects
            Output: a scalar value representing the density of genes
        '''
        if len(gene_list) == 0:
            return 0
        ids = "','".join([x.id for x in gene_list])
        scores = np.array(list(itertools.chain(*self.db.cursor().execute('''
            SELECT score FROM coex WHERE gene_a IN ('{}') AND gene_b IN ('{}')
            '''.format(ids,ids)
        ).fetchall())))
        if len(scores) == 0:
            return 0
        return (scores.mean()/(1/np.sqrt(len(scores))))

    def seed(self, gene_list, max_show = 65): 
        ''' Input: given a set of nodes, add on the next X strongest connected nodes ''' 
        if len(gene_list) == 0:
            return []
        neighbors = self.neighbors_score(gene_list)
        seed_set =  set([x.id for x in gene_list]).union(neighbors.index[0:min(len(neighbors),max_show)])
        return seed_set

    def gene_expression_vals(self,gene_list,zscore=True):
        ''' Input: A list of COBGenes
            Output: A dataframe containing normalized expression values for gene_list
        '''
        expr = pd.DataFrame(
            data = self.db.cursor().execute(
            '''SELECT gene,accession,value FROM expression WHERE gene in ('{}')
            '''.format("','".join([x.id for x in gene_list])), 
            ).fetchall(),columns=['gene','accession','value']
        ).pivot(index="accession",columns="gene",values='value')
        if zscore:
            expr = expr.apply(lambda x: (x-x.mean())/x.std(),axis=0)
        return expr
        
    def graph(self,gene_list):
        ''' Input: a gene list
            Output: a iGraph object '''
        if len(gene_list) == 0:
            return []
        # retrieve all first neighbors
        edges = self.subnetwork(gene_list)
        nodes = list(set(edges.source).union(edges.target))
        idmap = dict()
        for i,node in enumerate(nodes):
            idmap[node] = i
        graph = ig.Graph([ (idmap[source],idmap[target]) for source,target,score in edges[['source','target','score']].itertuples(index=False)])
        #graph.add_vertex(name=i,label=node)
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
        f.show()

        return {'ordered' : D, 'rorder' : Z1['leaves'], 'corder' : Z2['leaves']} 

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

    def compare_to_dat(self,filename,sep="\t",score_cutoff=3):
        ''' Compare the number of genes with significant edges as well as degree with a DAT file '''
        self.log("Reading in {}",filename)
        ref_dat = pd.read_table(filename,sep=sep,names=['gene_a','gene_b','score'])
        self.log("Retrieving network edges...")
        dat = pd.DataFrame(self.db.cursor().execute('''SELECT gene_a, gene_b, score 
            FROM  coex 
            WHERE score >= ?;
            ''',(score_cutoff,)).fetchall(),
            columns=['gene_a','gene_b','score']
        )
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
            return (gene,len(neighbors),len(ref_neighbors))
        return pd.DataFrame(
                    [compare_edges(gene) for gene in genes.intersection(ref_genes)],
                    columns = ['gene','degree','ref_degree']
                )

    def __repr__(self):
        return '''
            COB Dataset: {} - {} - {}
                Desc: {}
                FPKM: {}
                Num Genes: {}
                Num Accessions: {}
                Num Sig. Interactions: {}
        '''.format(self.name, self.organism, self.build,
                self.description,
                self.FPKM == 'FPKM',
                self.num_genes,
                self.num_accessions,
                self.num_sig_interactions
        )

class COBBuilder(Camoco):
    def __init__(self,name,description,FPKM,gene_build,organism,basedir="~/.camoco"):
        super().__init__(name,description,type="COB",basedir=basedir)
        self._global('FPKM',FPKM)
        self._global('build',gene_build)
        self._global('organism',organism)
        self.expr = matrix(0)
        self.genes = []
        self._create_tables()

    def add_accession(self,name,type="",description=""):
        ''' adds an accession to the database, useful for updating accessions with 
            details after importing them from a DataFrame '''
        self.db.cursor().execute(''' 
            INSERT OR UPDATE INTO accessions (name,type,description) VALUES (?,?,?)''',(name,type,description)
        )

    def from_DataFrame(self,df,transform=np.arctanh,significance_thresh=3):
        ''' assumes genes as rows and accessions as columns '''
        genes = df.index 
        accessions = df.columns
        expr = df.as_matrix()
        self.log("Imported {} genes and {} accession",len(genes),len(accessions))
        if self.FPKM == 'FPKM':
            self.log("Performing FPKM normalization using arcsinh")
            expr = np.arcsinh(expr) 
        else:
            self.log("WARNING! Make sure you normalized according to best practices for {}",self.FPKM)
        cur = self.db.cursor()
        try:
            cur.execute("BEGIN TRANSACTION")
            self._drop_indices()
            self.log("Adding Accession Values")
            cur.executemany("INSERT OR IGNORE INTO accessions (name) VALUES (?)",
                [(a,) for a in accessions]) 
            self.log("Adding Expression Values")
            cur.executemany("INSERT OR IGNORE INTO expression (gene,accession,value) VALUES (?,?,?)",
                # This comprehension iterates over matrix cells
                [(gene,acc,val) for acc,genevals in df.iteritems() for gene,val in genevals.iteritems()] 
            )
            self.log("Calculating Coexpression")
            tbl = pd.DataFrame(
                list(itertools.combinations(genes,2)),
                columns=['gene_a','gene_b']
            )
            # Now add coexpression data
            tbl['score'] = self._coex(expr)
            tbl['score'][tbl['score'] == 1] = 0.99999999 
            tbl['score'] = transform(tbl['score'])
            tbl['score'] = (tbl['score']-tbl.score.mean())/tbl.score.std()
            tbl['significant'] = pd.Series(list(tbl['score'] >= significance_thresh),dtype='int_')
            cur.executemany(
                ''' INSERT INTO coex (gene_a,gene_b,score,significant) VALUES (?,?,?,?)''',
                map(lambda x: (x[0],x[1],x[2],int(x[3])), tbl.itertuples(index=False))
            )
            cur.execute("END TRANSACTION")
            self.log("Building indices")
            self._build_indices()
            self.log("Done")
        except Exception as e:
            self.log("Something bad happened:{}",e)
            cur.execute("ROLLBACK")
        
    def _coex(self,matrix,parallelize=True):
        if parallelize: 
            progress = progress_bar(multiprocessing.Manager())
            pool = multiprocessing.Pool()
            scores = np.array(list(itertools.chain(
                *pool.imap(_PCCUp,[(i,matrix,progress) for i in range(len(matrix))])
            )))
            pool.close()
            pool.join()
            return scores
        else:
            scores = itertools.imap( _PCCUp,
                ( (i,self.genes,self.expr,UseGramene)
                 for i in range(self.num_genes))
            )

    def from_csv(self,filename,sep="\t"):
        ''' assumes genes as rows and accessions as columns '''
        df = pd.read_table(filename,sep=sep)
        self.from_DataFrame(df)
    
    def _build_indices(self):
        self.db.cursor().execute(''' 
            CREATE INDEX IF NOT EXISTS coex_abs ON coex (gene_a); 
            CREATE INDEX IF NOT EXISTS coex_gene_b ON coex (gene_b); 
            CREATE INDEX IF NOT EXISTS coex_gen_ab ON coex (gene_a,gene_b); 
            CREATE INDEX IF NOT EXISTS coex_score ON coex (score); 
            CREATE INDEX IF NOT EXISTS coex_significant ON coex (significant); 
            ANALYZE;
        ''')

    def _drop_indices(self):
        self.db.cursor().execute(''' 
            DROP INDEX IF EXISTS coex_gene_a;
            DROP INDEX IF EXISTS coex_gene_b;
            DROP INDEX IF EXISTS coex_score;
            DROP INDEX IF EXISTS coex_significant;
        ''')

    def _create_tables(self):
        cur = self.db.cursor()
        cur.execute('PRAGMA page_size = 1024;') 
        cur.execute('PRAGMA cache_size = 100000;') 
        cur.execute(''' 
            CREATE TABLE IF NOT EXISTS accessions (
                name TEXT,
                type TEXT,
                description TEXT
            );
            CREATE TABLE IF NOT EXISTS expression (
                gene TEXT,
                accession TEXT,
                value REAL
            );
            CREATE TABLE IF NOT EXISTS coex (
                gene_a TEXT,
                gene_b TEXT,
                score REAL,
                significant INTEGER 
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
    return [pearsonr(m[i,:],m[j,:])[0] for j in range(i+1,len(m))] 
