#!/usr/bin/python3

from camoco.Ontology import Ontology
from camoco.Camoco import Camoco
from camoco.COB import COB

from collections import defaultdict
from pandas import DataFrame
import networkx as nx
import numpy as np
import re
import os


class CompCOB(Camoco):
    def __init__(self, name):
        super().__init__(name, type='CompCOB')
        if self.cob:
            self.cob = COB(self.cob)
        if self.ont:
            self.ont = Ontology(self.ont)
        try:
            # Open the HDF5 store
            self.hdf5 = self._hdf5(name)
            self.densities = self.hdf5['densities']
        except KeyError:
            self.densities = DataFrame()

    @classmethod
    def compare(cls, cob, ont, name, description, min_members=3,
                max_members=300, random_trials=100, sig_threshold=0.05,
                debug=False):
        self = super().create(name, description, type='CompCOB')

        # Assign the COB and Ontology classes
        self.cob = cob
        self.ont = ont

        # Save the arguments for later retreival
        self._global('cob', cob.name)
        self._global('ont', ont.name)
        self._global('min_members', min_members)
        self._global('max_members', max_members)
        self._global('sig_threshold', sig_threshold)
        self._global('random_trials', random_trials)

        # Finding the terms and their genes
        self.ont_loci = self._find_term_loci(
            min_members=min_members, max_members=max_members
        )

        # Runnign the actual comparison
        self.densities = self._run_comparison(
            random_trials=random_trials,
            sig_threshold=sig_threshold, debug=debug
        )

        # Save the results
        self.log('Writing the results to the database.')
        # self.hdf5['ont_loci'] = DataFrame(self.ont_loci.items())
        self.hdf5['densities'] = self.densities
        return self

    def _find_term_loci(self, min_members=3, max_members=300):
        self.log('Finding the ontology terms.')
        ont_list = self.ont.terms()
        ont_loci = dict()
        self.log('Found {} total terms in the onotlogy.', len(ont_list))
        self.log('Filtering the terms and finding their genes.')
        count = 0
        for term in ont_list:
            # Filtering out the terms with too few or too many genes
            filt_loci = [x for x in term.locus_list if x in self.cob]
            if len(filt_loci) >= min_members and len(filt_loci) <= max_members:
                ont_loci[term.id] = filt_loci
                count += 1
        self.log('Found {} terms elegible for analysis.', count)
        return ont_loci

    def _run_comparison(self, random_trials=100,
                        sig_threshold=0.05, debug=False):
        self.log(
            'Running ' + str(random_trials) +
            ' random sets for each term and comparing them.')
        dens = dict()
        significant_terms = 0
        n = 0

        # Use only 25 terms for testing purposes
        if debug:
            ont_loci = list(self.ont_loci.items())[:25]
        else:
            ont_loci = list(self.ont_loci.items())

        # Iterate through all terms
        for term, loci in ont_loci:
            # Log how many teerms are done, to ensure people it hasn't crashed
            if n % 100 == 0:
                self.log('Compared {}/{} terms so far.', n, len(ont_loci))

            real_density = self.cob.density(loci)
            dens[term] = [real_density]

            # Run the random samples
            loci_count = len(loci)
            loci_tot_list = self.cob.refgen.random_genes(
                loci_count*random_trials
            )
            scores = []
            aboves = 0
            for x in range(random_trials):
                # Get the random genes
                loci_list = [loci_tot_list.pop() for x in range(loci_count)]

                # Find the density and save it
                score = self.cob.density(loci_list)
                if score >= real_density:
                    aboves += 1
                scores.append(score)

            # Add on the states from the scores
            dens[term].append(np.mean(scores))
            dens[term].append(np.std(scores))
            dens[term].append(aboves)

            # Figure out if that makes it significant
            if dens[term][-1] <= (random_trials*sig_threshold):
                dens[term].append(1)
                significant_terms += 1
            else:
                dens[term].append(0)
            n += 1
        self.log('Compared all {} terms.', n)

        # Convert the dict to a DataFrame
        ans = DataFrame.from_items(
            dens.items(),
            columns=[
                self.cob.name+' Density', 'Random Density Mean',
                'Random STD', 'Items >= '+self.cob.name, 'Significant'
            ],
            orient='index'
        )

        self.log('Number of Significant Terms: ' + str(significant_terms))
        self.log('Number Random Significants Expected: '+str(len(dens)*0.05))
        return ans

    def to_cytoscape(self,obo_file,out_file):
        self.log('Importing OBO: {}',obo_file)
        terms = defaultdict(dict)
        cur_term = ''
        alt_terms = []
        alts = set()
        isa_re = re.compile('is_a: (.*) !.*')
        with open(obo_file,'r') as INOBO:
            for line in INOBO:
                line = line.strip()
                if line.startswith('id: '):
                    for alt in alt_terms:
                        terms[alt] = terms[cur_term].copy()
                    alt_terms = []
                    cur_term = line.replace('id: ','')
                elif line.startswith('name: '):
                    terms[cur_term]['name'] = line.replace('name: ','')
                    terms[cur_term]['desc'] = ''
                    terms[cur_term]['is_a'] = []
                    terms[cur_term]['genes'] = set()
                elif line.startswith('namespace: '):
                    terms[cur_term]['type'] = line.replace('namespace: ','')
                elif line.startswith('alt_id: '):
                    alt_terms.append(line.replace('alt_id: ',''))
                    alts.add(line.replace('alt_id: ',''))
                elif line.startswith('def: '):
                    terms[cur_term]['desc'] += line.replace('def: ','')
                elif line.startswith('comment: '):
                    terms[cur_term]['desc'] += line.replace('comment: ','')
                elif line.startswith('is_a: '):
                    terms[cur_term]['is_a'].append(isa_re.match(line).group(1))

        def getISA(terms,term,alts):
            rels = []
            if terms[term]['is_a'] == []:
                return None
            else:
                for isa in terms[term]['is_a']:
                    if isa not in alts:
                        rels.append((term,isa))
                for isa in terms[term]['is_a']:
                    if isa not in alts:
                        x = getISA(terms,isa,alts)
                        if x:
                            rels += x
                return rels

        self.log('Adding the edges.')
        graph = nx.DiGraph()
        for term in self.densities.index:
            if term not in alts:
                x = getISA(terms,term,alts)
                if x:
                    graph.add_edges_from(x)

        self.log('Adding the node information.')
        all_nodes = set(list(graph.nodes()) + list(self.densities.index))
        for id in all_nodes:
            if id not in alts:
                if id in self.densities.index:
                    sig = int(self.densities['Significant'][id])
                else:
                    sig = 0
                graph.add_node(id, significant=sig, name=terms[id]['name'], desc=terms[id]['desc'])
        self.log('Done.')

        nx.write_graphml(graph,out_file)
        return graph

    def report(self, folder=''):
        # Figuring out the names for everything
        zipname = os.path.join(folder, self.cob.name+'EvalReport.zip')
        rawname = os.path.join(folder, self.cob.name+'RawEval.csv')
        termsname = os.path.join(folder, 'sigTerms.txt')
        readmename = os.path.join(folder, 'EvalREADME.txt')

        self.log('Writing the raw data.')
        self.densities.to_csv(path_or_buf=rawname)

        self.log('Writing the significant terms.')
        terms = self.densities[self.densities['Significant'] > 0]
        termsfile = open(termsname, 'w')
        for term in list(terms.index):
            termsfile.write(str(self.ont[term])+'\n')
        termsfile.close()

        self.log('Writing the README.')
        readme = open(readmename, 'w')
        readme.write((
            'This contains two files that were made from doing a comparison '
            'of density in the COB network based on GO terms vs. random '
            'trials. The method was just finding all of the terms that are in the'
            ' GO ontology that are anotated for maize, giving the list of genes to'
            ' the density function in the COB class and getting that back.'
            ' Then we picked the same ammount of random genes that are in the term'
            ' and did the density measurement. After repeating this 100 times,'
            ' we count how many densities are bigger than the real network, if'
            ' it is less than 5, it is considered significant.'

            'Contained are two files, the raw data and the list of terms '
            'that were considered significant. Here are some other relevant '
            'stats: Total Testable Terms: '+str(len(self.densities))+'\n'
            'Total Significant Terms: '+str(len(terms))+'\n'
            'The Expected Number of False Significants: '
            + str(len(self.densities)*0.05) + '\n'

            'Let me know if you want any other stats or information, or '
            'have suggestions as to what is or is not useful in a '
            'report about a dataset.'

            'Joe Jeffers'
            'jeffe174@umn.edu'
        ))
        readme.close()

        self.log('Zipping it all up.')
        os.system(
            'zip -j -q '+zipname+' '+rawname+' '+termsname+' '+readmename
        )
        os.system('rm '+rawname+' '+termsname+' '+readmename)
        self.log('Done, please check '+folder+' for your results.')
        return
