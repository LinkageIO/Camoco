import sys
import os
import math

import pandas as pd
import numpy as np
import camoco as co

from collections import OrderedDict
from camoco.Tools import log


class simulateGWAS(object):
    def __init__(self,cob=None,go=None):
        # GWAS Simulations needs a COB and an Gene Ontology
        self.cob = cob
        self.go = go
        self.method = 'density'
        # Give a place to store the results
        self.results = pd.DataFrame()

    def generate_output_name(self):
        # Handle the different output schemes
        if self.args.out is None:
            self.args.out = '{}_{}_{}_{}_{}'.format(
                self.cob.name,
                self.go.name,
                self.args.candidate_window_size,
                self.args.candidate_flank_limit,
                ':'.join(self.args.terms)
            )
        self.args.out = os.path.splitext(self.args.out)[0] + '.{}.GWASsimulation.tsv'.format(self.args.method)
        if os.path.dirname(self.args.out) != '':
            os.makedirs(os.path.dirname(self.args.out),exist_ok=True)
        if os.path.exists(self.args.out) and self.args.force != True:
            print(
                "Output for {} exists! Skipping!".format(
                    self.args.out
                ),file=sys.stderr
            )
            return

    def overlap(self,loci,bootstrap=False,iter_name=None):
        '''
            Calculate Network-Term Overlap based on the method in CLI args

            This method must return a DataFrame 
        '''
        # generate emirical and bootstrap values 
        if self.args.method == 'density':
            return self.cob.trans_locus_density(
                loci,
                flank_limit=self.args.candidate_flank_limit,
                by_gene = True,
                bootstrap = bootstrap,
                iter_name = iter_name
            ) 
        elif self.args.method == 'locality':
            return self.cob.trans_locus_locality(
                loci,
                flank_limit=self.args.candidate_flank_limit,
                by_gene=True,
                bootstrap=bootstrap,
                iter_name=iter_name,
                include_regression=True,
            ).rename(columns={'resid':'score'})


    def generate_bootstraps(self, loci, overlap):
        '''
            Bootstrapping procedure. Our target here is to provide enough bootstraps
            to identify loci that are significant at n==1000 bootstraps. The auto 
            procedure will continue untill we meet n==1000 OR we find 50 bootstraps
            that are have higher score such that we will never be significant at 1000
            boostraps (0.05 * 1000 = 50).
        '''
        target_score = overlap.score.mean()
        max_bs = 1000
        num_bs = 0
        bs = []
        if self.args.num_bootstraps == 'auto':
            # Create a bullshit generator... err bootstraps
            bs_generator = (self.overlap(loci,bootstrap=True,iter_name=x) \
                for x in range(max_bs)
            )
            while num_bs <= 50 and len(bs) < 1000: 
                # Add 50 bootstraps to current bootstraps
                bs = [next(bs_generator) for x in range(50)] + bs
                # Find the number of bs more extreme than the empirical score
                num_bs = sum([df.score.mean() >= target_score for df in bs])
                log(
                    "Iteration: {} -- current pval: {} {}% complete",
                    len(bs), num_bs/len(bs), max(len(bs)/10,100*(num_bs/50))
                )
        else:
            # Be a lil broke back noodle and explicitly bootstrap
            bs = [self.overlap(loci,bootstrap=True,iter_name=x) \
                for x in range(int(self.args.num_bootstraps))
            ]
        return pd.concat(bs)

    def simulate_missing_candidates(self, eloci, MCR=0):
        '''
            Simulate the effects of Missing Candidate Rates (MCR)

            Parameters
            ----------
            eloci : (iterable of loci objects)
                Set of true starting loci 
            MCR : (float - default = 0)
                Missing candidate rate. Can provide either a 
                floating point percentage or a whole number.
                (i.e. method will convert 30 <-> 0.30)

            Returns
            -------
            an iterable of loci
        '''
        # Convert between percentage and float
        if MCR < 1 and MCR > 0:
            MCR = MCR * 100
        # MCR: Remove a percentage of SNPs to simulate false negatives
        if MCR > 0:
            # Calulate the index needed to hit percent missing
            missing_index = math.ceil(len(eloci) * (1-(MCR/100)))
            if missing_index < 2:
                missing_index = 2 
            new_eloci = np.random.permutation(eloci)[0:missing_index]
            log('Simulating {}% of SNPs missed by GWAS ({} SNPs -> {})',MCR,len(eloci),len(new_eloci))
            eloci = new_eloci
        return eloci

    def simulate_false_candidates(self, eloci, FCR=0):
        '''
            Simulate the effects of False Candidate Rates (MCR)

            Parameters
            ----------
            eloci : (iterable of loci objects)
                Set of true starting loci 
            FCR : (float - default = 0)
                False candidate rate. Can provide either a 
                floating point percentage or a whole number.
                (i.e. method will convert 30 <-> 0.30)

            Returns
            -------
            an iterable of loci
        '''
        # Convert between percentage and float
        if FCR < 1 and FCR > 0:
            FCR = FCR*100
        # FCR: Replace a percentage of SNPs with false positives
        if FCR > 0:                                            
            # replace some loci with random genes if FDR specified              
            num_fcr = math.ceil(len(eloci) * (FCR/101))               
            fcr_loci = self.cob.refgen.random_genes(num_fcr,window=self.args.candidate_window_size) 
            log('Simulating {}% of SNPs as false positive -> adding {} SNPs',FCR,len(fcr_loci))
            # permute and truncate the loci then add fcr loci                   
            eloci = np.concatenate([                                            
                eloci,
                np.array(list(fcr_loci))
            ])   
        return eloci

    @classmethod
    def from_CLI(cls,args):
        '''
            Implements an interface to the CLI to perform GWAS simulation
        '''
        self = cls()
        # Build the base objects
        self.args = args
        # Load camoco objects
        self.go = co.GOnt(self.args.GOnt)
        self.cob = co.COB(self.args.cob)
        self.generate_output_name()

        # Generate an iterable of GO Terms
        if 'all' in self.args.terms:
            # Create a list of all terms within the size specification
            terms = list(self.go.iter_terms(
                min_term_size=self.args.min_term_size,
                max_term_size=self.args.max_term_size
            ))
        elif os.path.exists(self.args.terms[0]):
            # If parameter is a filename, read term name from a filenamie
            terms = list([self.go[x.strip()] for x in open(args.terms[0]).readlines()])
        else:
            # Generate terms from a parameter list
            terms = list([ self.go[x] for x in self.args.terms ])
        # Iterate and calculate
        log("Simulating GWAS for {} GO Terms",len(terms))
        min_term_size = np.min([len(x) for x in terms])
        max_term_size = np.max([len(x) for x in terms])
        log("All terms are between {} and {} 'SNPs'", min_term_size, max_term_size)

        results = []
        for i,term in enumerate(terms):
            log('-'*75)
            window_size = self.args.candidate_window_size
            flank_limit =  self.args.candidate_flank_limit
            # Generate a series of densities for parameters
            num_genes = len([x for x in term.loci if x in self.cob])
            eloci = [ x for x  in term.effective_loci(
                window_size=window_size
            ) if x in self.cob ]
            eloci = self.simulate_missing_candidates(eloci,self.args.percent_mcr)
            eloci = self.simulate_false_candidates(eloci,self.args.percent_fcr)
            log(
                'GWAS Simulation {}: {} ({}/{} genes in {})',
                i,term.id,len(eloci),num_genes,self.cob.name
            )   
            # Make sure that the number of genes is adequate
            if num_genes > self.args.max_term_size:
                log("Too many genes... skipping")
                continue
            elif num_genes < self.args.min_term_size:
                log("Too few genes... skipping")
                continue
            elif num_genes == 0:
                continue
            # Generate candidate genes from the effecive loci 
            candidates = self.cob.refgen.candidate_genes(
                eloci,
                flank_limit=flank_limit
            )
            log(
                "SNP to gene mapping finds {} genes at window:{} bp, "
                "flanking:{} genes", len(candidates),
                self.args.candidate_window_size,
                self.args.candidate_flank_limit
            )
            overlap = self.overlap(eloci)
            # Dont bother bootstrapping on terms with overlap score below 0
            if overlap.score.mean() < 0:
                continue
            bootstraps = self.generate_bootstraps(eloci,overlap)
            bs_mean = bootstraps.groupby('iter').score.apply(np.mean).mean()
            bs_std  = bootstraps.groupby('iter').score.apply(np.std).mean()
            # Calculate z scores for density
            overlap['zscore'] = (overlap.score-bs_mean)/bs_std
            bootstraps['zscore'] = (bootstraps.score-bs_mean)/bs_std
            overlap_pval = (
                (sum(bootstraps.groupby('iter').apply(lambda x: x.score.mean()) >= overlap.score.mean()))\
                / len(bootstraps.iter.unique())
            )
            # Create a results object
            overlap['COB'] = self.cob.name
            overlap['Ontology'] = self.go.name
            overlap['Term'] = term.id
            overlap['WindowSize'] = self.args.candidate_window_size
            overlap['FlankLimit'] = self.args.candidate_flank_limit
            overlap['FCR'] = args.percent_fcr
            overlap['MCR'] = args.percent_mcr
            overlap['NumRealGenes'] = num_genes
            overlap['NumEffective'] = len(eloci)
            overlap['NumCandidates'] = len(candidates)
            overlap['TermSize'] = len(term)
            overlap['TermCollapsedLoci'] = len(eloci)
            overlap['TermPValue'] = overlap_pval
            overlap['NumBootstraps'] = len(bootstraps.iter.unique())
            overlap['Method'] = self.args.method
            results.append(overlap.reset_index())

        self.results = pd.concat(results)
        self.results.to_csv(args.out,sep='\t',index=False)
