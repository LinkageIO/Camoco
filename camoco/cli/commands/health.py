import pandas as pd
import camoco as co
import numpy as np
import powerlaw

from os import path

from camoco.Tools import log as coblog

import matplotlib
matplotlib.style.use('ggplot')
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

    if args.edge_zscore_cutoff is not None:
        log("Changing Z-Score cutoff to {}",args.edge_zscore_cutoff)
        cob.set_sig_edge_zscore(args.edge_zscore_cutoff)

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
            include_accession_labels=True,
            raw=True,
            cluster_method=None
        )
    else:
        log('Skipped raw.')
    if not path.exists('{}_Expr_norm.png'.format(args.out)):
        cob.plot(
            '{}_Expr_norm.png'.format(args.out),
            include_accession_labels=True,
            raw=False,
            cluster_method='leaf',
            cluster_accessions=True
        )
    else:
        log('Skipped norm.')
    log('Plotting Cluster Expression-----------------------------------------')
    if not path.exists('{}_Expr_cluster.png'.format(args.out)):
        cob.plot(
            '{}_Expr_cluster.png'.format(args.out),
            include_accession_labels=True,
            raw=False,
            cluster_accessions=True,
            avg_by_cluster=True
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
            gene_qc = cob._bcolz('qc_gene')
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
        #Using powerlaw makes run-time warning the first time you use it.
        #This is still an open issue on the creators github.
        #The creator recommends removing this warning as long as there is a fit.
        np.seterr(divide='ignore', invalid='ignore')
        fit = powerlaw.Fit(degree,discrete=True,xmin=1)
        # get an axis
        ax = plt.subplot()
        # Calculate log ratios
        t2p = fit.distribution_compare('truncated_power_law', 'power_law')
        t2e = fit.distribution_compare('truncated_power_law', 'exponential')
        p2e = fit.distribution_compare('power_law','exponential')
        # Plot!
        emp = fit.plot_ccdf(ax=ax,color='r',linewidth=3, label='Empirical Data')
        pwr = fit.power_law.plot_ccdf(ax=ax,linewidth=2, color='b', linestyle=':', label='Power law')
        tpw = fit.truncated_power_law.plot_ccdf(ax=ax,linewidth=2, color='k', linestyle='-.', label='Truncated Power')
        exp = fit.exponential.plot_ccdf(ax=ax,linewidth=2, color='g', linestyle='--', label='Exponential')
        ####
        ax.set_ylabel("p(Degree≥x)")
        ax.set_xlabel("Degree Frequency")
        ax.legend(
            loc='best'
        )
        plt.title('{} Degree Distribution'.format(cob.name))
        # Save Fig
        try:
            plt.savefig('{}_DegreeDist.png'.format(args.out))
        except FutureWarning as e:
            # This is a matplotlib bug
            pass
    else:
        log('Skipping Degree Dist.')

    log('Plotting GO --------------------------------------------------------')
    if args.go is not None:
        if not path.exists('{}_GO.csv'.format(args.out)):
            go = co.GOnt(args.go)
            term_ids = []
            density_emp = []
            density_pvals  = []
            locality_emp = []
            locality_pvals = []
            term_sizes = []
            term_desc = []
            terms_tested = 0
            if args.max_terms is not None:
                log('Limiting to {} GO Terms',args.max_terms)
                terms = go.rand(
                    n=args.max_terms,
                    min_term_size=args.min_term_size,
                    max_term_size=args.max_term_size
                )
            else:
                terms = go.iter_terms(
                    min_term_size=args.min_term_size,
                    max_term_size=args.max_term_size
                )
            for term in terms:
                term.loci = list(filter(lambda x: x in cob, term.loci))
                if len(term) < args.min_term_size or len(term) > args.max_term_size:
                    continue
                #set density value for two tailed go so we only test it once
                density = cob.density(term.loci)
                #one tailed vs two tailed test
                if args.two_tailed_GO is False:
                    #run one tail for only positive values
                    if density > 0:
                        density_emp.append(density)
                    #skip negative density values
                    else:
                        continue
                #if two_tailed_go is not none
                else:
                    density_emp.append(density)
                term_ids.append(term.id)
                term_sizes.append(len(term))
                term_desc.append(str(term.desc))
                # ------ Density 
                # Calculate PVals
                density_bs = np.array([
                    cob.density(cob.refgen.random_genes(n=len(term.loci))) \
                    for x in range(args.num_bootstraps)
                ])
                if density > 0:
                    pval = sum(density_bs >= density)/args.num_bootstraps
                else:
                    pval = sum(density_bs <= density)/args.num_bootstraps
                density_pvals.append(pval)

                # ------- Locality
                locality = cob.locality(
                    term.loci,include_regression=True
                ).resid.mean()
                locality_emp.append(locality)
                # Calculate PVals
                locality_bs = np.array([
                    cob.locality(
                        cob.refgen.random_genes(n=len(term.loci)),
                        include_regression=True
                    ).resid.mean() \
                    for x in range(args.num_bootstraps)
                ])
                if locality > 0:
                    pval = sum(locality_bs >= locality)/args.num_bootstraps
                else:
                    pval = sum(locality_bs <= locality)/args.num_bootstraps
                locality_pvals.append(pval)
                # -------------
                terms_tested += 1
                if terms_tested % 100 == 0 and terms_tested > 0:
                    log('Processed {} terms'.format(terms_tested)) 
            go_enrichment = pd.DataFrame({
                'GOTerm' : term_ids,
                'desc' : term_desc,
                'size' : term_sizes,
                'density' : density_emp,
                'density_pval' : density_pvals,
                'locality' : locality_emp,
                'locality_pval' : locality_pvals
            })
            go_enrichment\
                .sort_values(by='density_pval',ascending=True)\
                .to_csv('{}_GO.csv'.format(args.out),index=False)
            if terms_tested == 0:
                log.warn('No GO terms met your min/max gene criteria!')
        else:
            go_enrichment = pd.read_table('{}_GO.csv'.format(args.out),sep=',')

        if not path.exists('{}_GO.png'.format(args.out)):
            # Convert pvals to log10
            with np.errstate(divide='ignore'):
                # When no bootstraps are more extreme than the term, the minus log pval yields an infinite
                go_enrichment['density_pval'] = -1*np.log10(go_enrichment['density_pval'])
                go_enrichment['locality_pval'] = -1*np.log10(go_enrichment['locality_pval'])
                # Fix the infinites so they are plotted
                max_density = np.max(go_enrichment['density_pval'][np.isfinite(go_enrichment['density_pval'])])
                max_locality = np.max(go_enrichment['locality_pval'][np.isfinite(go_enrichment['locality_pval'])])
                go_enrichment.loc[np.logical_not(np.isfinite(go_enrichment['density_pval'])),'density_pval'] = max_density + 1
                go_enrichment.loc[np.logical_not(np.isfinite(go_enrichment['locality_pval'])),'locality_pval'] = max_locality + 1
            plt.clf()
            figure,axes = plt.subplots(3,2,figsize=(12,12))
            # Calculate alpha based on one- or two-tails
            if args.two_tailed_GO == True:
                alpha = 0.05 / 2
            else:
                alpha = 0.05
            # -----------
            # Density
            # ----------
            axes[0,0].scatter(
                go_enrichment['density'],
                go_enrichment['density_pval'],
                alpha=0.05,
                color='blue'
            )
            axes[0,0].set_xlabel('Empirical Density (Z-Score)')
            axes[0,0].set_ylabel('Bootstraped -log10(p-value)')
            fold = sum(np.array(go_enrichment['density_pval'])>np.log10(alpha))/(alpha * len(go_enrichment))
            axes[0,0].axhline(y=-1*np.log10(0.05),color='red')
            axes[0,0].text(
                min(axes[0,0].get_xlim()),-1*np.log10(alpha),
                '{:.3g} Fold Enrichement'.format(fold),
                color='red'
            )
            axes[1,0].scatter(
                go_enrichment['size'],
                go_enrichment['density_pval'],
                alpha=0.05,
                color='blue'
            )
            axes[1,0].set_ylabel('Bootstrapped -log10(p-value)')
            axes[1,0].set_xlabel('Term Size')
            axes[1,0].axhline(y=-1*np.log10(alpha),color='red')
            axes[2,0].scatter(
                go_enrichment['size'],
                go_enrichment['density'],
                alpha=0.05,
                color='blue'
            )
            axes[2,0].scatter(
                go_enrichment.query('density_pval>1.3')['size'],
                go_enrichment.query('density_pval>1.3')['density'],
                alpha=0.05,
                color='r'
            )
            axes[2,0].set_ylabel('Density')
            axes[2,0].set_xlabel('Term Size')
            # ------------
            # Do Locality
            # ------------
            axes[0,1].scatter(
                go_enrichment['locality'],
                go_enrichment['locality_pval'],
                alpha=0.05,
                color='blue'
            )
            axes[0,1].set_xlabel('Empirical Locality (Residual)')
            axes[0,1].set_ylabel('Bootstraped -log10(p-value)')
            fold = sum(np.array(go_enrichment['locality_pval'])>1.3)/(0.05 * len(go_enrichment))
            axes[0,1].axhline(y=-1*np.log10(0.05),color='red')
            axes[0,1].text(
                min(axes[0,1].get_xlim()),-1*np.log10(0.05),
                '{:.3g} Fold Enrichement'.format(fold),
                color='red'
            )
            axes[1,1].scatter(
                go_enrichment['size'],
                go_enrichment['locality_pval'],
                alpha=0.05,
                color='blue'
            )
            axes[1,1].set_xlabel('Term Size')
            axes[1,1].set_ylabel('Bootstrapped -log10(p-value)')
            axes[1,1].axhline(y=-1*np.log10(0.05),color='red')
            axes[2,1].scatter(
                go_enrichment['size'],
                go_enrichment['locality'],
                alpha=0.05,
                color='blue' 
            )
            axes[2,1].scatter(
                go_enrichment.query('locality_pval>1.3')['size'],
                go_enrichment.query('locality_pval>1.3')['locality'],
                alpha=0.05,
                color='r'
            )
            axes[2,1].set_ylabel('Locality')
            axes[2,1].set_xlabel('Term Size')
            # Save Figure
            plt.tight_layout()
            try:
                plt.savefig('{}_GO.png'.format(args.out))
            except FutureWarning as e:
                pass
        else:
            log('Skipping GO Volcano.')
            
