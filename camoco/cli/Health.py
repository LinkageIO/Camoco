import camoco as co
import numpy as np
from os import path

import matplotlib.pylab as plt

def cob_health(args):
    cob = co.COB(args.cob)
    if args.out is None:
        args.out = '{}_Health'.format(cob.name)

    cob.log('Plotting Scores') #----------------------------------------------
    if not path.exists('{}_CoexPCC_raw.png'.format(args.out)):
        cob.plot_scores(
            '{}_CoexPCC_raw.png'.format(args.out),
            pcc=True
        )
    if not path.exists('{}_CoexScore_zscore.png'.format(args.out)):
        cob.plot_scores(
            '{}_CoexScore_zscore.png'.format(args.out),
            pcc=False
        )
    cob.log('Plotting Expression') #------------------------------------------
    if not path.exists('{}_Expr_raw.png'.format(args.out)):
        cob.plot(
            '{}_Expr_raw.png'.format(args.out),
            raw=True,
            cluster_method=None
        )
    if not path.exists('{}_Expr_norm.png'.format(args.out)):
        cob.plot(
            '{}_Expr_norm.png'.format(args.out),
            raw=False
        )

    if args.RefGen is not None:
        if not path.exists('{}_qc_gene.txt'.format(args.out)):
            with open('{}.summary.txt'.format(args.out)) as OUT:
                # Print out the network summary
                print(cob.summary,file=OUT)
                # Print out the breakdown of QC Values
                gene_qc = cob.hdf5['qc_gene']
                gene_qc = gene_qc[gene_qc.pass_membership]
                gene_qc['chrom'] = ['chr'+str(Zm5bFGS[x].chrom) for x in gene_qc.index]
                x.groupby('chrom').agg(sum,axis=0).to_csv('{}_qc_gene.txt'.format(args.out))
        

    cob.log('Plotting GO') #------------------------------------------
    if args.GO is not None: #-------------------------------------------------
        if not path.exists('{}_GO.png'.format(args.out)):
            go = co.GOnt(args.GO)
            emp_z = []
            pvals  = []
            terms_tested = 0
            for term in go.iter_terms():
                term.loci = list(filter(lambda x: x in cob, term.loci))
                if len(term) < 3 or len(term) > 300:
                    continue
                terms_tested += 1
                emp = cob.density(term.loci)
                emp_z.append(emp)
                # Calculate PValue
                bs = np.array([
                    cob.density(cob.refgen.random_genes(n=len(term.loci))) \
                    for x in range(args.num_bootstraps)
                ])
                if emp > 0:
                    pval = sum(bs>=emp)/args.num_bootstraps
                else:
                    pval = sum(bs<=emp)/args.num_bootstraps
                pvals.append(pval)
            plt.scatter(emp_z,-1*np.log10(pvals))
            plt.xlabel('Empirical Z Score')
            plt.ylabel('Bootstraped -log10(p-value)')
            fold = sum(np.array(pvals)<=0.05)/(0.05 * terms_tested)
            plt.title('{} x {}')
            plt.text(
                'right', 'bottom',
                '{} Fold Enrichement'.format(cob.name,go.name,fold),
            )
            plt.savefig('{}_GO.png'.format(args.out))
            
