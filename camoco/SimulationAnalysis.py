import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib.pylab as plt
matplotlib.style.use('ggplot')
#matplotlib.use('Agg')

import glob as glob
import camoco as co
import pandas as pd
import numpy as np
pd.set_option('display.width',300)

#from ggplot import *

class SimulationAnalysis(object):

    def __init__(self,dir='./',sep='\t'):
        self.sep = sep
        self.build_data(dir=dir)

    # Read in Raw data
    def build_data(self,dir='./'):
        self.df = pd.concat(
            [pd.read_table(x,sep=self.sep) \
                for x in glob.glob(dir+'*GWASsimulation.tsv')
            ]
        ).fillna(0)
        for col in ['MCR','FCR','WindowSize','FlankLimit','COB','Term']:
            if col not in self.df.columns:
                self.df[col] = 0
        self.all_results = self.df
        self.df = pd.pivot_table(
            self.df,
            index=['COB','Term','WindowSize','FlankLimit','MCR','FCR'],
            dropna=True,
            aggfunc=np.mean
        )
        self.df = self.df.reset_index()
        # Add an id column
        self.df['id'] = self.df.COB+'/'+self.df.Term
        self.df['window_id'] = ['{}/{}'.format(x,y) for x,y in zip(self.df.FlankLimit,self.df.WindowSize)]

    def terms_with_signal(self, max_pval=0.05, min_pval=0,
                          min_term_size=0, max_term_size=10e10
                          ):
        '''
            Return only terms with a starting pvalue (no noise) 
            between the specified arguments.

            Parameters
            ----------
            min_pval : float (default: 0.0)
                The minimum pval that a term can have and still be considered
                to have 'signal'. This is useful to stratify terms into groups,
                for instance moderate signal: min_pval=0.01 (pval_max: 0.05)
            max_pval : float (default: 0.05)
                The maximum pval that a term can have and still be considered
                to have 'signal'
            min_term_size : int (default: 0)
                The minimum starting term size to consider
            max_term_size : int (default: 0)
                The maximum term size to consider
            max_WindowSize : int (default: 1)
                Stratify results terms by window_size as specified in the term
                candidate genes function (snp2genes mapping). Default will 
                constrain genes to only those which are truly in the simulated
                term.
            max_FlankLimit : int (default: 0)
                Stratify terms by flanklimit as specified in the term candidate
                genes function (snp2gene mapping). Default will constrain genes
                to only include those which are truly in the simulated term.
        '''
        fdr = self.df.set_index(['COB','Term','MCR','FCR'],drop=False).sort_index()
        # Split out groups
        terms_with_signal = fdr.query(
            'MCR == 0 and FCR==0 and FlankLimit==0 and WindowSize<=1 '
            'and TermPValue<{max_pval} '
            'and TermPValue>={min_pval} '
            'and NumCandidates>={min_term_size} '
            'and NumCandidates<={max_term_size}'.format(**locals())
        )
        # Only return terms that were significant at FCR of 0 
        return pd.concat(
            [fdr.loc[x,y,:,:] for x,y in terms_with_signal[['COB','Term']].values]
        )
    
    def plot_signal_vs_noise(self,filename=None, figsize=(16,6), alpha=0.05,
                             figaxes=None, label=None, noise_source='MCR',
                             max_pval=0.05, min_pval=0, 
                             normalize_num_sig=False,
                             min_term_size=0,max_term_size=10e10):
        '''
            Plot the number of terms significant at varying FCR/MCR levels.

            THIS IS FOR WHEN YOU HAVE THE SAME TERMS AND YOU VARY THEIR
            FCR/MCR OR NOISE LEVEL!!!
        '''
        df = self.terms_with_signal(
            max_pval=max_pval, min_pval=min_pval, 
            min_term_size=min_term_size, max_term_size=max_term_size
        )
        # Aggregate by COB and by noise source and calculate the 
        # number of GO terms that have a sig pval below alpha
        breakdown = df.reset_index(drop=True)\
            .groupby(['COB',noise_source])\
            .apply(lambda x: sum(x.TermPValue <= alpha))\
            .reset_index()\
            .rename(columns={0:'num_sig'}).copy()
        breakdown['norm'] = breakdown.groupby(['COB']).apply(
            lambda x: x.num_sig/x.num_sig.max()
        ).values
        # Get the number of unique cobs in the dataset
        cobs = breakdown.COB.unique()
        if figaxes == None:
            plt.clf()
            fig, axes = plt.subplots(1,len(cobs),figsize=figsize,sharey=True)
        else:
            fig,axes = figaxes
            if len(cobs) != len(axes):
                raise ValueError('Cannot overlay plots with different number of COBs')
        for i,cob in enumerate(cobs):
            data = breakdown[breakdown.COB==cob]
            if normalize_num_sig:
                signal = data['norm']
            else:
                signal = data['num_sig']
            # Plot Signal vs Noise
            axes[i].plot(
                np.append(data[noise_source],1),
                np.append(signal,0),
                label=label,
                marker='o'
            )
            axes[i].set_title("{} signal vs {}".format(cob,noise_source))
            if i == 0:
                if normalize_num_sig:
                    axes[i].set_ylim(0,1.05)
                    axes[i].set_ylabel('Percent Significant Terms')
                else:
                    axes[i].set_ylabel('Number Significant Terms')
            axes[i].set_xlabel(noise_source)
        lgd = axes[i].legend(bbox_to_anchor=(2.0,1))
        if filename is not None:
            plt.savefig(filename,bbox_extra_artists=(lgd,),bbox_inches='tight',dpi=300)
        print(breakdown)
        return (fig,axes)

    def plot_pval_heatmap(self,filename=None,pval_cutoff=0.05):
        '''
            Generate a heatmap based on Density PVal
        '''
        cobs = self.df.COB.unique()
        fig,axes = plt.subplots(len(cobs),1,sharex=True)
        fig.set_size_inches(8,11)
        signal = self.terms_with_signal()
        for i,cob in enumerate(cobs):
            data = pd.pivot_table(
                signal.ix[signal.COB==cob],
                columns=['WindowSize','FlankLimit'],
                index=['GO'],
                values='Density_pval'
            )
            #data[data > pval_cutoff] = np.nan
            #data[data < pval_cutoff] = 0
            axes[i].set_frame_on(False)
            cmap = plt.cm.Greens_r
            cmap.set_bad('white',1.0)
            im = axes[i].pcolormesh(
                np.ma.masked_invalid(data),
                cmap=cmap,
                #linewidth=0.1,
                edgecolor='lightgrey',
                vmin=0,
                vmax=0.05
            )
            # Make the layout more natural
            axes[i].set_ylabel(cob,fontsize=20)
            axes[i].yaxis.set_label_position("right")
            axes[i].invert_yaxis()
            axes[i].set(
                yticks=np.arange(len(data.index))+0.5
            )
            axes[i].yaxis.set_ticks_position('left')
            axes[i].xaxis.set_ticks_position('none')
            if i == 0:
                axes[i].xaxis.tick_top()
                axes[i].set_xticks(np.arange(len(data.columns))+0.5)
                axes[i].set_xticklabels(data.columns.values, rotation=90)
            axes[i].set_yticklabels(data.index.values)
        # Create a new axis to append the color bar to
        divider = make_axes_locatable(axes[len(cobs)-1])
        cax = divider.append_axes("bottom", size="5%", pad=0.05)
        cbar = fig.colorbar(
            im, cax=cax, orientation='horizontal', ticks=[0,0.025,0.05]
        )
        cbar.set_ticklabels(['0','0.025','≥0.05'])
        plt.tight_layout(pad=0.4,w_pad=0.5, h_pad=1.0)
        if filename is not None:
            plt.savefig(filename,dpi=300)
            plt.close()
        return (fig,axes)   
    
    # FCR Analysis
    # ----------------------------------------------
    @staticmethod
    def plot_all_terms(df,filename='SimulatedFCR_Signal_vs_Noise.png',figsize=(16,8)):
        # Get the number of COBs there are
        raise NotImplementedError()
        return ggplot(df.reset_index(drop=True), aes(x='FCR',y='-logpval',color='id')) +\
            geom_line() + geom_point(alpha=0.1) + xlab('FCR') +\
            ylab('-log10(PVal)') + facet_wrap('COB') + ggtitle('Signal/Noise')

    def plot_windowed_terms(self,filename=None,figsize=(16,8),alpha=0.05,
                            figaxes=None,label=None,bins=10,
                            max_pval=0.05, min_pval=0, 
                            normalize_num_sig=False,
                            min_term_size=0,max_term_size=10e10):
        '''
            Plot the number of terms that retain 
        '''
        df = self.terms_with_signal(
            max_pval=max_pval, min_pval=min_pval, 
            min_term_size=min_term_size, max_term_size=max_term_size
        )
        # Generate the effective FCR rates
        real = df.set_index('id').query('FlankLimit==0 and FCR==0.00').NumCandidates
        df['EffectiveFCR'] = 1 - (real[df.id].values / df.NumCandidates)
        df['FCRBin'] = df['EffectiveFCR'].apply(lambda x: int(x*bins))
        df = df.sort_values('FCRBin')
        # Get the number of columns
        cobs = df.COB.unique()
        # Get the number of 
        if figaxes == None:
            fig, axes = plt.subplots(2,len(cobs),figsize=figsize,sharex=True)
        else:
            fig,axes = figaxes
            if len(cobs) != axes.shape[1]:
                raise ValueError('Cannot overlay plots with different number of COBs')
        for i,cob in enumerate(cobs):
            data = df[df.COB==cob]
            #starting_signal = sum(data.TermPValue[data.FCRBin >= 0] <= alpha)
            starting_signal = len(data.Term.unique())
            fcr = list()
            signal = list()
            num = list()
            for bin in data.FCRBin.unique():
                fcr.append(bin/bins)
                # Find out how many terms still have PVals < 0.05
                if normalize_num_sig:
                    signal.append(
                        len(data[(data.FCRBin >= bin) & (data.TermPValue <= alpha)]\
                            .Term.unique())/starting_signal
                    )
                else:
                    signal.append(
                        len(data[(data.FCRBin >= bin) & (data.TermPValue <= alpha)]\
                            .Term.unique())
                    )
            fcr.append(1)
            signal.append(0)
            axes[0,i].plot(fcr,signal,label=label,marker='o')
            if normalize_num_sig:
                axes[0,i].set_ylim((0,1))
            axes[0,i].set_title("{} Signal vs. FCR".format(cob))
            # do a box-plot too
            boxes = []
            labels = []
            for wid,wd in data.query('WindowSize>1')\
                    .sort_values(['WindowSize','FlankLimit'])\
                    .groupby(['WindowSize','FlankLimit']):
                print('Plotting {}'.format(wid))
                boxes.insert(0,wd.EffectiveFCR)
                labels.insert(0,'/'.join(map(str,wid)))
            if filename is not None:
                axes[1,i].boxplot(boxes,labels=labels,vert=False)
                axes[1,i].set_xlabel('FCR')
            # if we are on the left, plot the label
            if i == 0:
                if normalize_num_sig:
                    axes[0,i].set_ylabel('Signal (Fraction Terms Significant)')
                else:
                    axes[0,i].set_ylabel('Signal (Number Terms Significant)')
                axes[1,i].set_ylabel('Windowing Parameters')
        lgd = axes[0,i].legend(bbox_to_anchor=(2.0,1))
        if filename is not None:
            plt.savefig(filename,bbox_extra_artists=(lgd,),bbox_inches='tight',dpi=300)
        return (fig,axes)

       
    
    # Window Analysis
    # ----------------------------------------------
    def plot_windowsize_flanklimit(self):
        x=self.df.query('WindowSize>1').reset_index()\
            .groupby(['COB','WindowSize','FlankLimit','FCR'])\
            .apply(lambda x: sum(x.Pval <= 0.05)/len(x))\
            .reset_index()
        x['id'] = ['{}:{}:{}'.format(x,y,z) \
            for x,y,z in zip(x.COB,x.WindowSize,x.FlankLimit)
        ]
        # Plot the bastards 
        pd.pivot_table(
                df[(df.COB=='ZmPAN')&(df.WindowSize>1)],
                index=['FCR'],
                columns=['COB','WindowSize','FlankLimit'],
                values='Pval',
                aggfunc=lambda x: sum(x<=0.05)
            ).plot(subplots=True,layout=(3,3),sharey=True,figsize=(20,10))
        plt.savefig('ZmPAN_NinePlot.png')
        pd.pivot_table(
                df[(df.COB=='ZmSAM')&(df.WindowSize>1)],
                index=['FCR'],
                columns=['COB','WindowSize','FlankLimit'],
                values='Pval',
                aggfunc=lambda x: sum(x<=0.05)
            ).plot(subplots=True,layout=(3,3),sharey=True,figsize=(20,10))
        plt.savefig('ZmSAM_NinePlot.png')
        pd.pivot_table(
                df[(df.COB=='ZmRoot')&(df.WindowSize>1)],
                index=['FCR'],
                columns=['COB','WindowSize','FlankLimit'],
                values='Pval',
                aggfunc=lambda x: sum(x<=0.05)
            ).plot(subplots=True,layout=(3,3),sharey=True,figsize=(20,10))
        plt.savefig('ZmRoot_NinePlot.png')
                
