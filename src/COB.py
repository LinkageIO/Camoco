#!/usr/bin/python3

import pandas as pd

from Camoco import Camoco
from RefGen import RefGen,RefGenBuilder
from numpy import matrix,arcsinh
from Locus import Locus,Gene
from Tools import memoize
from scipy.stats import pearsonr
from collections import defaultdict

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
from scipy.spatial.distance import pdist, squareform, euclidean

import collections

class COB(Camoco):
    def __init__(self,name=None,basedir="~/.camoco"):
        if name is None:
            self.log('You must provide a name')
        else:
            super().__init__(name=name,type="COB",basedir=basedir)
            self.refgen = RefGen(self.refgen) 

    @property
    @memoize
    def num_genes(self):
        return self.db.cursor().execute("SELECT COUNT(DISTINCT(gene)) FROM expression").fetchone()[0]

    @property
    @memoize
    def num_sig_edges(self):
        return self.db.cursor().execute("SELECT COUNT(*) FROM coex WHERE significant = 1;").fetchone()[0]

    @property
    @memoize
    def num_edges(self):
            return self.db.cursor().execute("SELECT COUNT(*) FROM coex").fetchone()[0]

    @property
    @memoize
    def num_accessions(self):
        return self.db.cursor().execute("SELECT COUNT(DISTINCT(accession)) FROM expression").fetchone()[0]
   
    @property
    @memoize
    def num_sig_interactions(self):
        return self.db.cursor().execute("SELECT COUNT(*) FROM coex WHERE significant = 1;").fetchone()[0]

    @memoize
    def accessions(self):
        ''' returns a list of accessions used to build COB '''
        return self.db.cursor().execute("SELECT * FROM accessions").fetchall()

    @memoize
    def genes(self):
        ''' returns a list of genes used to build COB '''
        return self.refgen.from_ids([x[0] for x in 
            self.db.cursor().execute('SELECT DISTINCT(gene) FROM expression').fetchall()
        ])

    @memoize
    def sig_edges(self,limit=None):
        ''' returns a dataframe containing significant edges '''
        if not limit:
            return pd.DataFrame(
                self.db.cursor().execute(
                '''SELECT gene_a,gene_b,score
                FROM coex WHERE significant = 1''').fetchall(),
                columns = ['gene_a','gene_b','score']
            )
        else:
            return pd.DataFrame(
                self.db.cursor().execute(
                '''SELECT gene_a,gene_b,score
                FROM coex ORDER BY score DESC LIMIT ?''',(limit,)).fetchall(),
                columns = ['gene_a','gene_b','score']
            )

    @memoize
    def sig_genes(self,limit=None):
        sig_edges = self.sig_edges(limit=limit)
        return self.refgen.from_ids(set(sig_edges.gene_a).union(sig_edges.gene_b))

    def coexpression(self,gene_a,gene_b):
        ''' returns coexpression z-score between two genes '''
        return self.db.cursor().execute(
            ''' SELECT score FROM coex WHERE 
                (gene_a = ? AND gene_b = ?) 
                OR (gene_b = ? AND gene_a = ?)'''
            ,(gene_a.id,gene_b.id,gene_a.id,gene_b.id)
        ).fetchone()[0]
        

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

    def density(self,gene_list,pos_limit=50000,return_mean=True,trans_only=True):
        ''' calculates the denisty of the ENTIRE network amongst genes
            not within a certain distance of each other. This corrects for
            cis regulatory elements increasing noise in coexpression network '''
        # filter for only genes within network
        gene_list = list(filter(lambda x: x in self, gene_list))
        scores = []
        if trans_only:
            for gene_a,gene_b in itertools.combinations(gene_list,2):
                if abs(gene_a-gene_b) > pos_limit:
                    scores.append((gene_a,gene_b,self.coexpression(gene_a,gene_b)) )
        else:
            ids = "','".join([x.id for x in gene_list])
            scores = self.db.cursor().execute('''
                SELECT gene_a,gene_b,score FROM coex WHERE gene_a IN ('{}') AND gene_b IN ('{}')
                '''.format(ids,ids)
            ).fetchall()
        if len(scores) == 0:
            return np.nan
        if len(scores) == 1:
            return scores[0][2]
        scores = pd.DataFrame(scores,columns=['gene_a','gene_b','score'])
        if return_mean:
            return (scores.score.mean()/((scores.score.std())/np.sqrt(len(scores))))
        else:
            return scores

    def bootstrap_density(self,gene_list,miniter=20,maxiter=200,trans_only=True,pos_limit=50000):
        ''' Returns the number of bootstrapped densities higher that the 
            empirical value '''
        gene_list = list(filter(lambda x: x in self, gene_list))
        if not gene_list:
            return (np.nan,np.nan)
        emp_density = self.density(gene_list,trans_only=trans_only,pos_limit=pos_limit)
        bootstraps = []
        while True:
            bootstraps.append(emp_density < self.density(
                itertools.chain.from_iterable(
                    [self.refgen.bootstrap_flanking_genes(gene,gene_filter=self) for gene in gene_list]
                )
            ))
            if len(bootstraps) % miniter == 0 and np.sum(bootstraps) > 0.2*maxiter:
                break
            if len(bootstraps) >= maxiter:
                break
        return (np.sum(bootstraps)/len(bootstraps),len(bootstraps))
                

    def seed(self, gene_list, limit=65): 
        ''' Input: given a set of nodes, add on the next X strongest connected nodes ''' 
        if len(gene_list) == 0:
            return []
        neighbors = self.neighbors_score(gene_list)
        seed_set =  self.refgen.from_ids(
            set([x.id for x in gene_list]).union(neighbors.index[0:min(len(neighbors),limit)])
        )
        return seed_set

    def gene_expression_vals(self,gene_list,accession_list=None,zscore=True,transform=None):
        ''' Input: A list of COBGenes
            Output: A dataframe containing normalized expression values for gene_list
        '''
        expr = pd.DataFrame(
            data = self.db.cursor().execute(
            '''SELECT gene,accession,value FROM expression WHERE gene in ('{}')
            '''.format("','".join([x.id for x in gene_list])), 
            ).fetchall(),columns=['gene','accession','value']
        ).pivot(index="accession",columns="gene",values='value')
        if accession_list:
            expr = expr.ix[accession_list]
        if transform:
            expr = transform(expr)
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


    def heatmap(self,dm,filename=None,figsize=(16,16), maskNaNs=True, cluster_x=True, cluster_y=True,
        cluster_method="euclidian", title=None, heatmap_unit_label='Expression Z Score'):
        ''' Draw clustered heatmaps of an expression matrix'''
        D = np.array(dm) 
        row_labels = dm.index
        col_labels = dm.columns
        f = plt.pylab.figure(figsize=figsize,facecolor='white')
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
            plt.pylab.title(title)
        vmax = max(np.nanmin(abs(D)),np.nanmax(abs(D)))
        vmin = vmax*-1
        self.log("Extremem Values: {}",vmax)
        cmap = self.cmap
        im = axmatrix.matshow(D, aspect='auto', cmap=cmap,vmax=vmax,vmin=vmin)
        # Handle NaNs
        if maskNaNs:
            nan_mask = np.ma.array(D,mask=np.isnan(D))
            cmap.set_bad('grey',1.0)
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
        # Add color bar
        axColorBar = f.add_axes([0.09,0.75,0.2,0.05])
        f.colorbar(im,orientation='horizontal',cax=axColorBar,
            ticks=np.arange(np.ceil(vmin),np.ceil(vmax),int((vmax-vmin)/2))
        )
        plt.pylab.title(heatmap_unit_label)
        if filename:
            plt.pyplot.savefig(filename)
            plt.pyplot.close()
        return pd.DataFrame(
            data=D,
            index=row_labels,
            columns=col_labels
        )

  
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
            self.db.cursor().executemany('''INSERT OR REPLACE INTO coor (gene,x,y) VALUES (?,?,?)''',
                [(gene.id,layout[i][0],layout[i][1]) for i,gene in enumerate(gene_list)]
            )

    def mcl(self,gene_list,I=2.0):
        from subprocess import Popen, PIPE
        MCLstr = "\n".join(["{}\t{}\t{}".format(a,b,c) for a,b,c in self.subnetwork(gene_list).itertuples(index=False)])
        cmd = "mcl - --abc -scheme 6 -I {} -o -".format(I)
        p = Popen(cmd, stdout=PIPE, stdin=PIPE, stderr=self.log_file, shell=True)
        sout,serr = p.communicate(MCLstr.encode('utf-8'))
        p.wait()
        if p.returncode==0 and len(sout)>0:
            return [ self.refgen.from_ids([gene.decode('utf-8') for gene in line.split()]) for line in sout.splitlines() ]
        else:
            self.log( "MCL return code: {}".format(p.returncode))

    def to_dat(self,gene_list,filename):
        with open(filename, 'w') as OUT:
            for a,b,c in self.subnetwork(gene_list).itertuples(index=False):
                print("{}\t{}\t{}".format(a,b,c),file=OUT)
 
    def plot(self,gene_list,filename=None,width=3000,height=3000,**kwargs):
        if not filename:
            filename = "{}.png".format(self.name)
        ig.plot(self.graph(gene_list),layout=self.coordinates(gene_list),target=filename,bbox=(width,height),**kwargs)

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
                self.FPKM == 'FPKM',
                self.num_genes,
                self.num_accessions,
                self.num_sig_edges/self.num_edges
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

    def __contains__(self,obj):
        try:
            if self.db.cursor().execute('''
                SELECT COUNT(*) FROM expression WHERE gene = ?''',(obj.id,)).fetchone()[0] > 0:
                return True
            else:
                return False
        except Exception as e:
            pass
        try:
            # Can be a string object
            if self.db.cursor().execute('''
                SELECT COUNT(*) FROM expression WHERE gene = ?''',(str(obj),)).fetchone()[0] > 0:
                return True
            else:
                return False
        except Exception as e:
            pass
        self.log("object '{}' is not correct type to test for membership in {}",obj,self.name)
        

class COBBuilder(Camoco):
    def __init__(self, name, description, FPKM, refgen, organism, basedir="~/.camoco"):
        super().__init__(name,description,type="COB",basedir=basedir)
        self._global('FPKM',FPKM)
        self._global('refgen',refgen.name)
        self.expr = matrix(0)
        self.genes = []
        self._create_tables()

    def add_accession(self,name,type="",description=""):
        ''' adds an accession to the database, useful for updating accessions with 
            details after importing them from a DataFrame '''
        self.db.cursor().execute(''' 
            INSERT OR UPDATE INTO accessions (name,type,description) VALUES (?,?,?)''',(name,type,description)
        )

    def add_filter_refgen(self,refgen):
        ''' Filter the refgen to only contain genes available in COB. Only do this after the expression table
            has been populated!!'''
        filtered_refgen = RefGenBuilder.import_from_RefGen("{}RefGen".format(self.name),
            "RefGen for {} filtered from {}".format(self.name,refgen.name),
            gene_filter = self,
            basedir = self.basedir
        )
        # Remember for next time
        self._global('refgen',filtered_refgen.name)

    def from_DataFrame(self, df, transform=np.arctanh, significance_thresh=3, min_expr=1, 
        min_missing_data=0.2, min_single_sample_expr=5, membership=None ):
        ''' Import a COB dataset from a pandas dataframe. Assumes genes/traits as rows and accessions/genotypes as columns.
            Options:
                transform                   - a vectorizable (numpy) function used to normalize expression data. 
                                            Default is inverse hyperbolic sin which transforms larger values more than smaller ones. 
                                            see DOI:10.1371/journal.pone.0061005
                significance_threshold      - a Z-score threshold at which each interaction will be called significant.
                min_expr                    - expr (usually FPKM) lower than this will be set as NaN and not used in correlation.
                min_missing_data            - The threshold of percent missing data for a gene/trait to be thrown out.
                min_single_sample_expr      - The minimum threshold a gene/trait needs to meet in a single sample
                                              to be considered expressed and worth correlating.
                membership                  - a class instance or object implementing the contains method
                                              in which to filter out genes which the object does not
                                              contain. Useful with RefGen objects to get rid of unmapped genes.
            '''
        # ----- Perform QC steps ----
        # Remember build options for later
        self._global('build_transform',transform.__name__)
        self._global('build_significance_thresh',significance_thresh)
        self._global('build_min_expr',min_expr)
        self._global('build_min_missing_data',min_missing_data)
        self._global('build_min_single_sample_expr',min_single_sample_expr)
        self._global('build_membership',str(membership))
        # Filter out genes which are not in the membership object
        if membership:
            self.log("Filtering out genes not in {}",membership)
            df = df[[x in membership for x in df.index]]
        # Set minimum FPKM threshold
        self.log("Filtering out expression values lower than {}",min_expr)
        df.loc[df < min_expr] = np.nan
        # filter out genes with too much missing data
        self.log("Filtering out genes with > {} missing data",min_missing_data)
        df = df[df.apply(lambda x : ((sum(np.isnan(x))) < len(x)*min_missing_data),axis=1)] 
        # filter out genes which do not meet a minimum expr threshold in at least one sample
        self.log("Filtering out genes which do not have one sample above {}",min_single_sample_expr)
        df = df[df.apply(lambda x: any(x >= min_single_sample_expr),axis=1)]
        # ----- import into Database -----
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
            self.log("Building Database")
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
            CREATE INDEX IF NOT EXISTS expression_gene_name ON expression (gene);
            CREATE INDEX IF NOT EXISTS coor_gene ON coor(gene);
            CREATE INDEX IF NOT EXISTS coor_x ON coor(x);
            CREATE INDEX IF NOT EXISTS coor_y ON coor(y);
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
        vals.append(pearsonr(m[i,mask],m[j,mask])[0] )
    return vals
