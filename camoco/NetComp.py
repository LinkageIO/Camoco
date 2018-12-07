import types
import six
import random
import pandas as pd
import logging

from .COB import COB
from .Camoco import Camoco

from functools import wraps 
from collections.abc import Iterable
from itertools import permutations
from minus80 import Freezable
from contextlib import contextmanager


def accepts_iterable(fn):
    '''
    This decorator detects when an iterable is passed into the method call
    instead of a single element (first arg only). Then, instead of calling the
    function on the arg (normally) the method is applied to each element of the
    iterable. This essentially turns fn(arg) into [f(x) for x in arg]. 
    '''
    @wraps(fn)
    def wrapped(self,arg,*args,**kwargs):
        # 
        if (isinstance(arg,Iterable) and not \
            isinstance(arg, six.string_types)):
            return [fn(self,x,*args,**kwargs) for x in arg]      
        else:
            fn(self,arg,*args,**kwargs)
    return wrapped


class NetComp(Freezable):

    def __init__(self,name,networks=None):
        # init core objects
        super().__init__(name=name)
        self._initialize_tables()
        self.log = logging.getLogger(f'NetComp.{name}')
        self.networks = set()
        # Retrieve stored networks
        net_names = [x[0] for x in self._db.cursor().execute('''
            SELECT name FROM networks
        ''')]
        for name in net_names:
            self.networks.add(COB(name))
        # Handle args
        if networks is None:
            networks = []
        for n in networks:
            self.add_network(n)

    def _initialize_tables(self):
        cur = self._db.cursor()
        cur.execute('''
            CREATE TABLE IF NOT EXISTS networks (
                name TEXT
            )
        ''')

        cur.execute('''
            CREATE TABLE IF NOT EXISTS net_comp (
                source TEXT,
                target TEXT,
                source_cluster TEXT,
                n_genes INT,
                comp_method TEXT,
                source_coex FLOAT,
                target_coex FLOAT,
                target_coex_pval FLOAT
            ) 
        ''')

    @accepts_iterable
    def add_network(self,net):
        '''
            Add a network (COB) to the 
            NetComp object.
        '''
        if isinstance(net,str):
            net = COB(net)
        if not isinstance(net,COB):
            raise ValueError(f'a valid network must be provided')
        # check if name is already in networks
        if net.name not in [x.name for x in self.networks]:
            self.networks.add(net)
        self._db.cursor().execute('''
            INSERT INTO networks(name) VALUES (?)
        ''',(net.name,))


    @contextmanager 
    def _net_comp_buffer(self):
        '''
            Add a comparison record to the internal database
        '''
        cur = self._db.cursor()
        record_buffer = []
        yield record_buffer
        cur.executemany('''
            INSERT INTO net_comp VALUES (?,?,?,?,?,?,?,?)
        ''',record_buffer)


    def compare_cluster_coex(self,min_cluster_size=10,max_cluster_size=300,num_bootstrap=100,method='density'):
        '''
            Compare the co-expression of genes within clusters between networks.
            This will compute the strength of coexpression of genes from clusters
            in network A in network B. This comparison tests that clusters in network
            A are *also* clustered in network B.
        '''
        with self._net_comp_buffer() as results:
            for source,target in permutations(self.networks,2):
                print(f'comparing {source.name} to {target.name}')
                if source.name == target.name:
                    continue 
                common_genes = set(source.genes()).intersection(target.genes())
                for cid in source.clusters.cluster.unique():
                    cluster_genes = source.cluster_genes(cid)
                    cluster_genes = common_genes.intersection(cluster_genes)
                    if (len(cluster_genes) >= max_cluster_size \
                    or len(cluster_genes) <= min_cluster_size):
                        continue
                    # Calculate Co-expression
                    if method == 'density':
                        source_coex = source.density(cluster_genes)
                        target_coex = target.density(cluster_genes)
                        target_coex_pval = sum(
                            [target.density(random.sample(common_genes,len(cluster_genes))) >= target_coex \
                            for _ in range(num_bootstrap)]
                        ) / num_bootstrap
                    elif method == 'locality':
                        source_coex = source.locality(cluster_genes)
                        target_coex = target.locality(cluster_genes,include_regression=True).resid.mean()
                        target_coex_pval = sum(
                            [target.locality(random.sample(common_genes,len(cluster_genes))).resid.mean() >= target_coex \
                            for _ in range(num_bootstrap)]
                        ) / num_bootstrap
                    results.append((
                        source.name,target.name,str(cid),len(cluster_genes),
                        method, source_coex, target_coex, target_coex_pval
                    ))

    def _percent_shared_mcl(self, pval_cutoff=0,comp_method='density'):
        results = []
        sig_clusters = self._query(f'SELECT * FROM net_comp WHERE comp_method = "{comp_method}" and target_coex_pval <= {pval_cutoff}')
        sig_clusters['cluster'] = sig_clusters.source + '_' + sig_clusters.source_cluster
        
        for a,b in permutations(sig_clusters.target.unique(),2):
            # Get the number that are in both

            # Check each 
            a_sig = set(sig_clusters.query(f'target == "{a}" and source != "{a}" and source != "{b}"').cluster)
            b_sig = set(sig_clusters.query(f'target == "{b}" and source != "{b}" and source != "{a}"').cluster)
            sig_both = set(
                sig_clusters.query(f'target == "{b}" and source == "{a}"'))\
                .union(sig_clusters.query(f'target == "{a}" and source == "{b}"')
            )
            total = a_sig.union(b_sig).union(sig_both)
            numerator = len(a_sig.intersection(b_sig)) + len(sig_both)
            denominator = len(total)
            results.append((a,b,numerator/denominator))
        results = pd.DataFrame(results,columns=['source','target','percent_sig'])
        return pd.pivot_table(results,index='source',columns='target')
