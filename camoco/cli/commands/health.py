import pandas as pd
import camoco as co
import numpy as np
import powerlaw

from os import path

from camoco.Tools import log as coblog

import matplotlib.pylab as plt

def cob_health(args):
    log = coblog()
    log('\n'
        '-----------------------\n'
        '   Network Health      \n'
        '-----------------------\n'
    )
    cob = co.COB(args.cob)
    if args.out is None:
        args.out = '{}_Health'.format(cob.name)

    log('Plotting Scores ----------------------------------------------------') 
    if not path.exists('{}_CoexPCC_raw.png'.format(args.out)):
        cob.plot_scores(
            '{}_CoexPCC_raw.png'.format(args.out),
            pcc=True
        )
    else:
        log('Skipped Raw.')

    if not path.exists('{}_CoexScore_zscore.png'.format(args.out)):
        cob.plot_scores(
            '{}_CoexScore_zscore.png'.format(args.out),
            pcc=False
        )
    else:
        log('Skipped Norm.')

    log('Plotting Expression ------------------------------------------------')
    if not path.exists('{}_Expr_raw.png'.format(args.out)):
        cob.plot(
            '{}_Expr_raw.png'.format(args.out),
            raw=True,
            cluster_method=None
        )
    else:
        log('Skipped raw.')

    if not path.exists('{}_Expr_norm.png'.format(args.out)):
        cob.plot(
            '{}_Expr_norm.png'.format(args.out),
            raw=False
        )
    else:
        log('Skipped norm.')

    log('Printing Summary ---------------------------------------------------')
    if not path.exists('{}.summary.txt'.format(args.out)):
        with open('{}.summary.txt'.format(args.out),'w') as OUT:
            # Print out the network summary
            cob.summary(file=OUT)
    else:
        log('Skipped summary.')

    log('Printing QC Statistics ---------------------------------------------')
    if args.refgen is not None:
        if not path.exists('{}_qc_gene.txt'.format(args.out)):
                # Print out the breakdown of QC Values
                refgen = co.RefGen(args.refgen)
                gene_qc = cob.hdf5['qc_gene']
                gene_qc = gene_qc[gene_qc.pass_membership]
                gene_qc['chrom'] = ['chr'+str(refgen[x].chrom) for x in gene_qc.index]
                gene_qc = gene_qc.groupby('chrom').agg(sum,axis=0)
                # Add totals at the bottom
                totals = gene_qc.ix[:,slice(1,None)].apply(sum)
                totals.name = 'TOTAL'
                gene_qc = gene_qc.append(totals)
                gene_qc.to_csv('{}_qc_gene.txt'.format(args.out),sep='\t')
        else:
            log('Skipped QC summary.')

    #if not path.exists('{}_CisTrans.png'.format(args.out)):
        # Get trans edges

    log('Plotting Degree Distribution ---------------------------------------')
    if not path.exists('{}_DegreeDist.png'.format(args.out)):
        degree = cob.degree['Degree'].values
        fit = powerlaw.Fit(degree,discrete=True,xmin=1)
        # get an axis
        ax = plt.subplot()
        # Calculate log ratios
        t2p = fit.distribution_compare('truncated_power_law', 'power_law')
        t2e = fit.distribution_compare('truncated_power_law', 'exponential')
        p2e = fit.distribution_compare('power_law','exponential')
        # Plot!
        emp = fit.plot_ccdf(ax=ax,color='r',linewidth=3, label='Empirical Data')
        pwr = fit.power_law.plot_ccdf(ax=ax, color='b', linestyle='--', label='Power law')
        tpw = fit.truncated_power_law.plot_ccdf(ax=ax, color='k', linestyle='--', label='Truncated Power')
        exp = fit.exponential.plot_ccdf(ax=ax, color='g', linestyle='--', label='Exponential')
        ####
        ax.set_ylabel("p(Degree≥x)")
        ax.set_xlabel("Degree Frequency")
        ax.legend(
            loc='best'
        )
        plt.title('{} Degree Distribution'.format(cob.name))
        # Save Fig
        plt.savefig('{}_DegreeDist.png'.format(args.out))
    else:
        log('Skipping Degree Dist.')

    log('Plotting GO --------------------------------------------------------')
    if args.go is not None:
        if not path.exists('{}_GO.csv'.format(args.out)):
            go = co.GOnt(args.go)
            term_ids = []
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
                term_ids.append(term.id)
                pvals.append(pval)
            go_enrichment = pd.DataFrame({
                'id':term_ids,
                'density' : emp_z,
                'pval':pvals
            })
            go_enrichment.sort('pval',ascending=True)\
                .to_csv('{}_GO.csv'.format(args.out),index=False)
        else:
            go_enrichment = pd.read_table('{}_GO.csv'.format(args.out))

        if not path.exists('{}_GO.png'.format(args.out)):
            plt.clf()
            plt.scatter(go_enrichment['density'],-1*np.log10(go_enrichment['pval']))
            plt.xlabel('Empirical Z Score')
            plt.ylabel('Bootstraped -log10(p-value)')
            fold = sum(np.array(pvals)<=0.05)/(0.05 * (go_enrichment.id.unique())
            plt.title('{} x {}'.format(cob.name,go.name))
            plt.axhline(y=-1*np.log10(0.05),color='red')
            plt.text(
                1, 0.1,
                '{:.3g} Fold Enrichement'.format(fold),
            )
            plt.savefig('{}_GO.png'.format(args.out))
        else:
            log('Skipping GO Volcano.')
            
