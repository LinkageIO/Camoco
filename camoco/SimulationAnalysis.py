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

from ggplot import *

tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)

class SimulationAnalysis(object):

    def __init__(self,dir='./',sep='\t'):
        self.sep = sep
        self.build_data(dir=dir)

    # Read in Raw data
    def build_data(self,dir='./'):
        self.df = pd.concat(
            [pd.read_table(x,sep=self.sep) \
                for x in glob.glob(dir+'*GWASSim.csv')]
        )
        if 'MCR' not in self.df.columns:
            self.df['MCR'] = 0
        if 'FCR' not in self.df.columns:
            self.df['FCR'] = 0
        self.df = pd.pivot_table(
            self.df,
            index=['COB','GO','WindowSize','FlankLimit','MCR','FCR'],
            aggfunc=np.mean
        )
        self.df = self.df.reset_index()
        # Add an id column
        self.df['id'] = self.df.COB+'/'+self.df.GO
        self.df['window_id'] = ['{}/{}'.format(x,y) for x,y in zip(self.df.FlankLimit,self.df.WindowSize)]
    
    def terms_with_signal(self, max_pval=0.05, min_pval=0,
			  pval_col='Density_pval'):
        '''
            Return only terms with a starting pvalue between the 
            specified argument.
        '''
        # Split the data in a few informative ways
        # We only want the True FCR
        fdr = self.df.set_index(['COB','GO','FCR'],drop=False).sort_index()
        # Split out groups
        terms_with_signal = fdr.query(
            'FCR==0 and MCR==0 and FlankLimit==0 and {}<{} and {}>={}'.format(pval_col,max_pval,pval_col,min_pval)
        )
        # Only return terms that were significant at FCR of 0 
        return pd.concat(
            [fdr.loc[x,y,:] \
                for x,y in zip(
                    terms_with_signal.COB,
                    terms_with_signal.GO
                )
            ]
        )

    def plot_signal_vs_noise(self,filename=None, figsize=(16,8), alpha=0.05,
                             figaxes=None, label=None, noise_source='FCR',
                             pval_col='Density_pval', max_pval=0.05, 
                             min_pval=0, normalize_num_sig=False):
        '''
            Plot the number of terms significant at varying FCR/MCR levels.

            THIS IS FOR WHEN YOU HAVE THE SAME TERMS AND YOU VARY THEIR
            FCR/MCR OR NOISE LEVEL!!!
        '''
        df = self.terms_with_signal(max_pval=max_pval,min_pval=min_pval)
        # Aggregate by COB and by noise source and calculate the 
        # number of GO terms that have a sig pval below alpha
        breakdown = df.reset_index(drop=True)\
            .groupby(['COB',noise_source])\
            .apply(lambda x: sum(x[pval_col] <= alpha))\
            .reset_index()\
            .rename(columns={0:'num_sig'})
        # Get the number of unique cobs in the dataset
        cobs = breakdown.COB.unique()
        if figaxes == None:
            fig, axes = plt.subplots(1,len(cobs),figsize=figsize,sharey=True)
        else:
            fig,axes = figaxes
            if len(cobs) != len(axes):
                raise ValueError('Cannot overlay plots with different number of COBs')
        for i,cob in enumerate(cobs):
            data = breakdown[breakdown.COB==cob]
            if normalize_num_sig:
                # Divide by the starting number of significant terms
                starting_signal = data.ix[data[noise_source]==data[noise_source].min()].iloc[0]['num_sig']
                data.loc[:,'num_sig'] = data.loc[:,'num_sig'] / starting_signal
            axes[i].plot(data[noise_source],data['num_sig'],label=label,marker='o')
            axes[i].set_title("{} Terms".format(cob))
            if i == 0:
                if normalize_num_sig:
                    axes[i].set_ylim(0,1.05)
                    axes[i].set_ylabel('Percent Significant Terms')
                else:
                    axes[i].set_ylabel('Number Significant Terms')
            axes[i].set_xlabel(noise_source)
        lgd = axes[i].legend(bbox_to_anchor=(2.0,0.5))
        if filename is not None:
            plt.savefig(filename,bbox_extra_artists=(lgd,),bbox_inches='tight')
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
        return ggplot(df.reset_index(drop=True), aes(x='FCR',y='-logpval',color='id')) +\
            geom_line() + geom_point(alpha=0.1) + xlab('FCR') +\
            ylab('-log10(PVal)') + facet_wrap('COB') + ggtitle('Signal/Noise')

    @staticmethod
    def plot_windowed_terms(df,filename=None,figsize=(16,8),alpha=0.05,
                            figaxes=None,label=None,bins=20):
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
            if len(cobs) != len(axes):
                raise ValueError('Cannot overlay plots with different number of COBs')
        for i,cob in enumerate(cobs):
            data = df[df.COB==cob]
            starting_signal = sum(data.Pval[data.FCRBin >= 0] <= alpha)
            fcr = list()
            signal = list()
            num = list()
            for bin in data.FCRBin.unique():
                fcr.append(bin/bins)
                signal.append(sum(data.Pval[data.FCRBin >= bin] <= alpha)/starting_signal)
                num.append(sum(data.Pval[data.FCRBin >= bin] <= alpha))
            axes[0,i].plot(fcr,signal,label=label,marker='o',color='k')
            axes[0,i].set_title("{} Signal vs. FCR".format(cob))
            # Plot the other axis too
            ax2 = axes[0,i].twinx()
            ax2.plot(fcr,num,label=label,marker='o',color='k')
            # do a scatter plot too
            boxes = []
            labels = []
            for wid,wd in data.query('WindowSize>1').sort_values(['WindowSize','FlankLimit']).groupby(['WindowSize','FlankLimit']):
                print('Plotting {}'.format(wid))
                boxes.insert(0,wd.EffectiveFCR)
                labels.insert(0,'/'.join(map(str,wid)))
            axes[1,i].boxplot(boxes,labels=labels,vert=False)
            # if we are on the left, plot the label
            if i == 0:
                axes[0,i].set_ylabel('Signal (Fraction Terms Significant)')
                axes[1,i].set_ylabel('Windowing Parameters')
        #plt.legend()
        if filename is not None:
            plt.savefig(filename)
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
                
