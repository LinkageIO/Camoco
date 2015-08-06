#!/usr/bin/python3

from camoco.Ontology import Term,Ontology
from camoco.Locus import Locus,Gene
from camoco.RefGen import RefGen
from camoco.Tools import log
from camoco.COB import COB

from pandas import DataFrame
import numpy as np

class CompCOB(object):
    def __init__(self,cob,ont,min_members=3,max_members=300):
        self.cob = cob
        self.ont = ont
        self.min = min_members
        self.max = max_members
        self.ont_loci = self._find_term_loci()

    def _find_term_loci(self):
        print('Finding the ontology terms.')
        ont_list = self.ont.terms()
        ont_loci = dict()
        print(len(ont_list))
        print('Finding the genes associated with the terms in the ontology.')
        count = 0
        for term in ont_list:
            # Filtering out the terms with too few or too many genes
            if len(term) >= self.min and len(term) <= self.max:
                ont_loci[term.id] = term.locus_list
                count += 1
        print('Found '+str(count)+' terms elegible for analysis.')
        return ont_loci

    def run_comparison(self,random_trials=100,sig_threshold=0.05,debug=False):
        print('Running '+str(random_trials)+' random sets for each term and comparing them.')
        dens = dict()
        dud_count = 0
        significant_terms = 0
        n = 1

        # Use only 25 terms for testing purposes
        if debug:
            ont_loci = list(self.ont_loci.items())[:25]
        else:
            ont_loci = list(self.ont_loci.items())

        # Iterate through all terms
        for term,loci in ont_loci:
            # Log how many teerms are done, to ensure people it hasn't crashed
            if n%100 == 0:
                print('Compared '+str(n)+' terms.')

            # Add the term to the dict and find the density
            loci = [x for x in loci if x in self.cob]
            if loci == []:
                dud_count += 1
                continue
            dens[term] = [self.cob.density(loci)]

            # Run the random samples
            loci_count = len(loci)
            scores = []
            for x in range(random_trials):
                # Get the random genes
                loci_list = []
                for i in range(loci_count):
                    loci_list.append(self.cob.refgen.random_gene())
                # Find the density and save it
                scores.append(self.cob.density(loci_list))

            # Add on the states from the scores
            dens[term].append(np.mean(scores))
            dens[term].append(np.std(scores))

            # Find how many random trials are more dense than the term
            aboves = 0
            for x in scores:
                if x >= dens[term][0]:
                    aboves += 1
            dens[term].append(aboves)

            # Figure out if that makes it significant
            if dens[term][-1] <= (random_trials*sig_threshold):
                dens[term].append(1)
                significant_terms += 1
            else:
                dens[term].append(0)
            n += 1
        print('Compared all '+str(n)+' terms.')

        # Convert the dict to a DataFrame
        ans = DataFrame.from_items(dens.items(),columns=[self.cob.name+' Density','Random Density Mean','Random STD','Items >= Test Mean', 'Significant'],orient='index')

        #self.densities = ans
        print('Number of Significant Terms: '+ str(significant_terms))
        print('Number Random Significants Expected: '+str(len(dens)*0.05))
        return ans
