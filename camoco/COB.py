#!/usr/bin/python3

import pandas as pd

from camoco.Camoco import Camoco
from camoco.RefGen import RefGen,RefGenBuilder
from camoco.Locus import Locus,Gene
from camoco.Tools import memoize
from numpy import matrix,arcsinh
from scipy.stats import pearsonr
from collections import defaultdict
from itertools import chain

import igraph as ig
import numpy as np
import pandas as pd
import time
import math as math
import multiprocessing as multiprocessing
import itertools
import matplotlib.pylab as plt
import matplotlib

from scipy.stats import hypergeom
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform, euclidean

import collections

class COB(Camoco):
    def __init__(self,name=None,basedir="~/.camoco"):
        if name is None:
            self.log('You must provide a name')
        else:
            try:
                super().__init__(name=name,type="COB",basedir=basedir)
                self.refgen = RefGen(self.refgen) 
            except Exception as e:
                return None

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
    def sig_edges(self,limit=None,inc_dis=False):
        ''' returns a dataframe containing significant edges '''
        if not limit:
            return pd.DataFrame(
                self.db.cursor().execute(
                '''SELECT gene_a,gene_b,score,distance
                FROM coex WHERE significant = 1''').fetchall(),
                columns = ['gene_a','gene_b','score','distance']
            )
        else:
            return pd.DataFrame(
                self.db.cursor().execute(
                '''SELECT gene_a,gene_b,score,distance
                FROM coex ORDER BY score DESC LIMIT ?''',(limit,)).fetchall(),
                columns = ['gene_a','gene_b','score','distance']
            )

    @memoize
    def sig_genes(self,limit=None):
        sig_edges = self.sig_edges(limit=limit)
        return list(self.refgen.from_ids(set(sig_edges.gene_a).union(sig_edges.gene_b)))

    def coexpression(self,gene_a,gene_b,trans_only=True):
        ''' returns coexpression z-score between two genes '''
        return self.db.cursor().execute(
            ''' SELECT score FROM coex WHERE 
                (gene_a = ? AND gene_b = ?) 
                OR (gene_b = ? AND gene_a = ?)'''
            ,(gene_a.id,gene_b.id,gene_a.id,gene_b.id)
        ).fetchone()[0]
        

    def neighbors(self,gene_list,min_distance=None,sig_only=True):
        ''' Input : a list of COBGene Objects
            Output: Returns the neighbors for a list of genes as a DataFrame
            Columns returned: query_name, queryID, neighbor_name, neighbor, score
        '''
        gene_list = "','".join([x.id for x in gene_list])
        cur = self.db.cursor()
        query = ''' 
            SELECT {} AS source, {} AS target, score, distance FROM coex 
            WHERE 
            {} IN ('{}')  
            ''' 
        if min_distance:
            query += ' AND distance > {} '.format(min_distance)
        if sig_only:
            query += ' AND significant = 1;'
        else:
            query += ';'
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
        gene_list = list(filter(lambda x: x in self, gene_list))
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

    def density(self,gene_list,return_mean=True,min_distance=None):
        ''' calculates the denisty of the ENTIRE network amongst genes
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

    def gene_expression_vals(self,gene_list,accession_list=None,zscore=True,transform=None):
        ''' Input: A list of COBGenes
            Output: A dataframe containing normalized expression values for gene_list
        '''
        try:
            expr = pd.DataFrame(
                data = self.db.cursor().execute(
                '''SELECT gene,accession,value FROM expression WHERE gene in ('{}')
                '''.format("','".join([x.id for x in gene_list])), 
                ).fetchall(),columns=['gene','accession','value']
            ).pivot(index="accession",columns="gene",values='value')
        except ValueError as e:
            self.log("No expression values present")
            return
        if accession_list:
            expr = expr.ix[accession_list]
        if transform:
            expr = transform(expr)
        if zscore:
            expr = expr.apply(lambda x: (x-x.mean())/x.std(),axis=0)
        return expr
        
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


    def heatmap(self,dm,filename=None,figsize=(16,16), maskNaNs=True, cluster_x=True, cluster_y=True,
        cluster_method="euclidian", title=None, heatmap_unit_label='Expression Z Score'):
        ''' Draw clustered heatmaps of an expression matrix'''
        D = np.array(dm) 
        row_labels = dm.index
        col_labels = dm.columns
        f = plt.figure(figsize=figsize,facecolor='white')
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
            plt.title(title)
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
        plt.title(heatmap_unit_label)
        if filename:
            plt.savefig(filename)
            plt.close()
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

    def mcl(self,gene_list,I=2.0,scheme=7,min_distance=100000):
        from subprocess import Popen, PIPE
        MCLstr = "\n".join(["{}\t{}\t{}".format(a,b,c) for a,b,c in \
            self.subnetwork(gene_list,min_distance=min_distance)[['source','target','score']].itertuples(index=False)])
        cmd = "mcl - --abc -scheme {} -I {} -o -".format(scheme,I)
        self.log("running MCL: {}",cmd)
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
 
    def plot(self,gene_list,filename=None,width=3000,height=3000,layout=None,**kwargs):
        if not filename:
            filename = "{}.png".format(self.name)
        if not layout:
            layout = self.coordinates(gene_list)
        ig.plot(self.graph(gene_list),layout=layout,target=filename,bbox=(width,height),**kwargs)

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
        heatmap_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',heatmapdict,256)
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

    def filter_refgen(self):
        ''' Filter the refgen to only contain genes available in COB. Only do this after the expression table
            has been populated!!'''
        filtered_refgen = RefGenBuilder.import_from_RefGen(
            "{}RefGen".format(self.name),
            "RefGen for {} filtered from {}".format(self.name,self.refgen.name),
            self.refgen,
            gene_filter = self,
            basedir = self.basedir
        )
        # Remember for next time
        self._global('refgen',filtered_refgen.name)



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
        

    def save(self, name, description, FPKM, refgen, organism, basedir="~/.camoco"):
        super().__init__(name,description,type="COB",basedir=basedir)
        self._global('FPKM',FPKM)
        self._global('refgen',refgen.name)
        self.expr = matrix(0)
        self.refgen = RefGen(self.refgen)
        self._create_tables()

    def add_accession(self,name,type="",description=""):
        ''' adds an accession to the database, useful for updating accessions with 
            details after importing them from a DataFrame '''
        self.db.cursor().execute(''' 
            INSERT OR UPDATE INTO accessions (name,type,description) VALUES (?,?,?)''',(name,type,description)
        )


    def from_DataFrame(self, df, transform=np.arctanh, significance_thresh=3, min_expr=1, 
        max_gene_missing_data=0.2, min_single_sample_expr=5, 
        max_accession_missing_data=0.5,membership=None,dry_run=False):
        ''' Import a COB dataset from a pandas dataframe. Assumes genes/traits as rows and accessions/genotypes as columns.
            Options:
                transform                   - a vectorizable (numpy) function used to normalize expression data. 
                                              Default is inverse hyperbolic sin which transforms larger values more than smaller ones. 
                                              see DOI:10.1371/journal.pone.0061005
                significance_threshold      - a Z-score threshold at which each interaction will be called significant.
                min_expr                    - expr (usually FPKM) lower than this will be set as NaN and not used in correlation.
                max_gene_missing_data       - The threshold of percent missing data for a gene/trait to be thrown out.
                                              Genes with missing more missing data than this threshold will be removed.
                min_single_sample_expr      - The minimum threshold a gene/trait needs to meet in a single sample
                                              to be considered expressed and worth correlating.
                max_accession_missing_data  - Accessions with more than this percent missing data will be removed.
                membership                  - a class instance or object implementing the contains method
                                              in which to filter out genes which the object does not
                                              contain. Useful with RefGen objects to get rid of unmapped genes.
            '''
        # ----- Perform QC steps ----
        # Remember build options for later
        self._global('build_transform',transform.__name__)
        self._global('build_significance_thresh',significance_thresh)
        self._global('build_min_expr',min_expr)
        self._global('build_max_gene_missing_data',max_gene_missing_data)
        self._global('build_min_single_sample_expr',min_single_sample_expr)
        self._global('build_membership',str(membership))
        self._global('build_max_accession_missing_data',max_accession_missing_data)
        # Filter out genes which are not in the membership object
        self.log('Starting set: {} genes {} accessions'.format(len(df.index),len(df.columns)))
        if membership:
            self.log("Filtering out genes not in {}",membership)
            df = df[[x in membership for x in df.index]]
        self.log('Kept: {} genes {} accessions'.format(len(df.index),len(df.columns)))
        # Set minimum FPKM threshold
        self.log("Filtering out expression values lower than {}",min_expr)
        df_flt = df.copy()
        df_flt[df < min_expr] = np.nan
        df = df_flt
        self.log('Kept: {} genes {} accessions'.format(len(df.index),len(df.columns)))
        # filter out genes with too much missing data
        self.log("Filtering out genes with > {} missing data",max_gene_missing_data)
        df = df.loc[df.apply(lambda x : ((sum(np.isnan(x))) < len(x)*max_gene_missing_data),axis=1),:] 
        self.log('Kept: {} genes {} accessions'.format(len(df.index),len(df.columns)))
        self.log("Filtering out accessions with > {} missing data",max_accession_missing_data)
        df = df.loc[:,df.apply(lambda x : ((sum(np.isnan(x))) < len(x)*max_accession_missing_data),axis=0)] 
        self.log('Kept: {} genes {} accessions'.format(len(df.index),len(df.columns)))
        # filter out genes which do not meet a minimum expr threshold in at least one sample
        self.log("Filtering out genes which do not have one sample above {}",min_single_sample_expr)
        df = df[df.apply(lambda x: any(x >= min_single_sample_expr),axis=1)]
        self.log('Kept: {} genes {} accessions'.format(len(df.index),len(df.columns)))
        #

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
        if dry_run:
            return
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
            tbl = pd.DataFrame(
                list(itertools.combinations(genes,2)),
                columns=['gene_a','gene_b']
            )
            # Now add coexpression data
            self.log("Calculating Coexpression")
            tbl['score'] = self._coex(expr)
            tbl['score'][tbl['score'] == 1] = 0.99999999 
            tbl['score'] = transform(tbl['score'])
            # Sometimes, with certain datasets, the NaN mask overlap completely for the
            # two genes expression data making its PCC a nan. This affects the mean and std fro the gene.
            valid_scores = np.ma.masked_array(tbl['score'],np.isnan(tbl['score']))
            tbl['score'] = (valid_scores-valid_scores.mean())/valid_scores.std()
            tbl['significant'] = pd.Series(list(tbl['score'] >= significance_thresh),dtype='int_')
            self.log("Calculating Gene Distance")
            tbl['distance'] = self._gene_distances()
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
        # update the reference geneom
        self.add_filter_refgen(self.refgen)
        
    def _coex(self,matrix):
        progress = progress_bar(multiprocessing.Manager())
        pool = multiprocessing.Pool()
        scores = np.array(list(itertools.chain(
            *pool.imap(_PCCUp,[(i,matrix,progress) for i in range(len(matrix))])
        )))
        pool.close()
        pool.join()
        return scores

    def _gene_distances(self):
        progress = progress_bar(multiprocessing.Manager())
        pool = multiprocessing.Pool()
        genes = list(self.genes())
        distances = np.array(list(itertools.chain(
            *pool.imap(_DISUp,[(i,genes,progress) for i in range(len(genes))])
        )))
        pool.close()
        pool.join()
        return distances


    @classmethod
    def from_csv(cls, name, description, FPKM, refgen, organism, filename, basedir="~/.camoco", sep="\t",**kwargs):
        ''' assumes genes as rows and accessions as columns '''
        df = pd.read_table(filename,sep=sep)
        cob = cls(name,description,FPKM,refgen,organism,basedir)
        cob.from_DataFrame(df,**kwargs)
   
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

   
