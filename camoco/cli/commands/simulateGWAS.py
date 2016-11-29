import sys
import os
import math

import pandas as pd
import numpy as np
import camoco as co
from collections import OrderedDict


class simulateGWAS(object):
    def __init__(self):
        self.results = pd.DataFrame()

    def generate_output_name(self):
        # Handle the different output schemes
        if self.args.out is None:
            self.args.out = '{}_{}_{}_{}_{}'.format(
                self.cob.name,
                self.ont.name,
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

    def generate_bootstraps(self,loci,overlap):
        '''
            bootstrapping procedure. Our target here is to provide enough bootstraps
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
                self.cob.log(
                    "Iteration: {} -- current pval: {} {}% complete",
                    len(bs),num_bs/len(bs),num_bs/50
                )
        else:
            # Be a lil broke back noodle and explicitly bootstrap
            bs = [self.overlap(loci,bootstrap=True,iter_name=x) \
                for x in range(int(self.args.num_bootstraps))
            ]
        return pd.concat(bs)

    def simulate_missing_candidates(self, eloci, MCR=None):
        '''
            Simulate the effects of Missing Candidate Rates (MCR)

            Parameters
            ----------
            eloci : (iterable of loci objects)
                Set of true starting loci 
            MCR : (float)
                Missing candidate rate. Can provide either a 
                floating point percentage or a whole number.
                (i.e. method will convert 30 <-> 0.30)

            Returns
            -------
            an iterable of loci
        '''
        # Convert between percentage and float
        if MCR < 1:
            MCR = MCR * 100
        # MCR: Remove a percentage of SNPs to simulate false negatives
        if MCR != None and MCR > 0:
            # Calulate the index needed to hit percent missing
            missing_index = math.ceil(len(eloci) * (1-(MCR/100)))
            if missing_index < 2:
                missing_index = 2 
            new_eloci = np.random.permutation(eloci)[0:missing_index]
            #cob.log('Simulating {}% of SNPs missed by GWAS ({} SNPs -> {})',MCR,len(eloci),len(new_eloci))
            eloci = new_eloci
        return eloci

    def simulate_false_candidates(self, eloci, FCR=None):
        '''
            Simulate the effects of False Candidate Rates (MCR)

            Parameters
            ----------
            eloci : (iterable of loci objects)
                Set of true starting loci 
            FCR : (float)
                False candidate rate. Can provide either a 
                floating point percentage or a whole number.
                (i.e. method will convert 30 <-> 0.30)

            Returns
            -------
            an iterable of loci
        '''
        if FCR < 1:
            FCR = FCR*100
        # FCR: Replace a percentage of SNPs with false positives
        if FCR != None and FCR > 0:                                            
            # replace some loci with random genes if FDR specified              
            num_fcr = math.ceil(len(eloci) * (FCR/101))               
            fcr_loci = self.cob.refgen.random_genes(num_fcr,window=window_size) 
            #cob.log('Simulating {}% of SNPs as false positive -> adding {} SNPs',FCR,len(fcr_loci))
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
        go = co.GOnt(self.args.GOnt)
        cob = co.COB(self.args.cob)
        self.generate_output_name()

        # Generate an iterable of GO Terms
        if 'all' in self.args.terms:
            # Create a list of all terms within the size specification
            terms = list(go.iter_terms(
                min_term_size=self.args.min_term_size,
                max_term_size=self.args.max_term_size
            ))
        elif os.path.exists(self.args.terms[0]):
            # If parameter is a filename, read term name from a filenamie
            terms = list([go[x.strip()] for x in open(args.terms[0]).readlines()])
        else:
            # Generate terms from a parameter list
            terms = list([ go[x] for x in self.args.terms ])
        # Iterate and calculate
        cob.log("Simulating GWAS for {} GO Terms",len(terms))
        min_term_size = np.min([len(x) for x in terms])
        max_term_size = np.max([len(x) for x in terms])
        cob.log("All terms are between {} and {} 'SNPs'", min_term_size, max_term_size)

        results = []
        for i,term in enumerate(terms):
            cob.log('-'*75)
            window_size = self.args.candidate_window_size
            flank_limit =  self.args.candidate_flank_limit
            # Generate a series of densities for parameters
            num_genes = len([x for x in term.loci if x in cob])
            eloci = [ x for x  in term.effective_loci(
                window_size=window_size
            ) if x in cob ]
            cob.log(
                'Simulation {}: {} ({}/{} genes in {})',
                i,term.id,len(eloci),num_genes,cob.name
            )   
            if num_genes > self.args.max_term_size:
                cob.log("Too many genes... skipping")
                continue
            elif num_genes < self.args.min_term_size:
                cob.log("Too few genes... skipping")
                continue
            elif num_genes == 0:
                continue
            
            candidates = cob.refgen.candidate_genes(
                eloci,
                flank_limit=flank_limit
            )
            cob.log(
                "SNP to gene mapping finds {} genes at window:{} bp, "
                "flanking:{} genes", len(candidates),
                self.args.candidate_window_size,
                self.args.candidate_flank_limit
            )
            if len(candidates) == 0:
                density = np.NaN
                trans_density = np.NaN
                density_pval = np.NaN
                locality = np.NaN
                trans_locality = np.NaN
                locality_pval = np.NaN
            else:
                # Density
                density = cob.density(
                    candidates
                )
                trans_density = cob.trans_locus_density(
                    eloci,
                    flank_limit=flank_limit
                )
                density_bs = np.array([
                    cob.trans_locus_density(
                        eloci,
                        flank_limit=flank_limit,
                        iter_name = x,
                        bootstrap=True
                    ) for x in range(args.num_bootstraps)
                ])
                density_pval = sum([x>=trans_density for x in density_bs])/args.num_bootstraps
                # Locality
                locality = cob.locality(
                    candidates, include_regression=True
                ).resid.mean()
                trans_locality = cob.trans_locus_locality(
                   eloci,
                   flank_limit=flank_limit,
                   include_regression=True
                ).resid.mean()
                locality_bs = np.array([
                    cob.trans_locus_locality(
                        eloci,
                        flank_limit=flank_limit,
                        iter_name=x,
                        bootstrap=True,
                        include_regression=True
                    ).resid.mean() for x in range(args.num_bootstraps)    
                ])
                locality_pval = sum([x>=trans_locality for x in locality_bs])/args.num_bootstraps
            results.append(OrderedDict([
               ('COB',cob.name),
               ('GO',term.id),
               ('TermSize',len(term)),
               ('NumRealGenes',num_genes),
               ('WindowSize',window_size),
               ('FlankLimit',flank_limit),
               ('FCR',args.percent_fcr),
               ('MCR',args.percent_SNPs_missed),
               ('NumEffective',len(eloci)),
               ('NumCandidates',len(candidates)),
               ('Density',density),
               ('TransDensity',trans_density),
               ('Density_pval',density_pval),
               ('Locality',locality),
               ('TransLocality',trans_locality),
               ('Locality_pval',locality_pval),
               ('NumBootstraps',args.num_bootstraps)
            ]))
        results = pd.DataFrame.from_dict(results)
        results.to_csv(args.out,sep='\t',index=False)
