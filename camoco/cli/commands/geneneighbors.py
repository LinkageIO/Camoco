#!/usr/bin/env python3
import camoco as co
import argparse
import sys
import os
import copy

import pandas as pd
import numpy as np
import scipy as sp

from itertools import chain

def geneneighbors(args):
    args.out = os.path.splitext(args.out)[0] + '_neighbors.txt'
    if os.path.dirname(args.out) != '':
        os.makedirs(os.path.dirname(args.out),exist_ok=True)
    if os.path.exists(args.out):
        print(
            "Output for {} exists! Skipping!".format(
                args.out
            ),file=sys.stderr
        )
        return

    cob = co.COB(args.cob)
    cob.set_sig_edge_zscore(int(args.zscore))
    genes = cob.refgen.iter_genes()
    print("Generating neighbors for {} ".format(
            cob.name,
        ),file=sys.stderr)
    #make empty list to store our results
    GENE = []
    NumNeighbors = []
    results = []
    #iterate through each gene in the network
    for i in genes:
        #get the list of neighbors
        NB = cob.neighbors(i)
        #pandas is weird so we need to get gene names like this
        NB = NB.reset_index()
        x = set(NB.gene_a)
        x = x.union(NB.gene_b)
        #sometimes it lists itself as a neighbor so remove
        if i.id in x:
            x.remove(i.id)
        else:
            continue
        #store how many neighbors a gene has
        NumNeighbors.append(str(len(x)))
        GeneNeighbors = []
        GENE.append(i.id)
        SCORE = []
        DIST = []
        JGENE = []
        #for each of the gene neighbors
        for j in x:
            #get the ID and see the co-expression results 
            #between the neighbor and original gene
            gene2 = cob.refgen.from_ids(j)
            y = cob.coexpression(i, gene2)
            score = y[0]
            significant = y[1]
            distance = y[2]
            #store all of the information
            JGENE.append(j)
            SCORE.append(score)
            DIST.append(distance)
        #zip those results so we can sort it together
        ZIPPED = zip(JGENE ,SCORE, DIST) 
        #sort by the Z-score 
        SZIPPED = sorted(ZIPPED, key=lambda x: x[1], reverse = True)
        #grab the top 10 genes
        TOP10 = SZIPPED[0:int(args.numneighbors)]
        #unzip them so we can combine them for writing
        x,y,z= zip(*TOP10)

        #combin the lists and add to results
        for l,m,n in zip(
        x,y,z):
            Temp = (str(l),str(m),str(n))
            NeighborInfo = ','.join(Temp)
            GeneNeighbors.append(NeighborInfo)
        GeneNeighbors = '\t'.join(GeneNeighbors)
        #print(GeneNeighbors)
        results.append(GeneNeighbors)

    #write to a file
    output = open(args.out, 'w')
    output.write('Gene'+'\t'+'Number of Neighbors'+'\t'+'Gene,ZScore,Significant,distance'+'\n')

    for a, b, c in zip(
        GENE, NumNeighbors, results):
        final = (a,b,c)
        output.write('\t'.join(final))
        output.write('\n')
