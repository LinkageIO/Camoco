import matplotlib
matplotlib.use('Agg')

import matplotlib.pylab as plt
matplotlib.style.use('ggplot')

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

class Analysis(object):

    def __init__(self,dir='./'):
        self.build_data(dir=dir)

    # Read in Raw data
    def build_data(self,dir='./'):
        self.df = pd.concat(
            [pd.read_table(x,sep='\t',index_col=0) \
                for x in glob.glob(dir+'*GWASSim.csv')]
        )
        if len(self.df) == 0:
            raise ValueError('No simulation output files were found.')
        self.df = pd.pivot_table(
            self.df.reset_index(),
            index=['COB','GO','WindowSize','FlankLimit','FCR'],
            aggfunc=np.mean
        ).reset_index()
        # Add an id column
        self.df['id'] = self.df.COB+'/'+self.df.GO
        self.df['window_id'] = ['{}/{}'.format(x,y) for x,y in zip(self.df.FlankLimit,self.df.WindowSize)]
        #self.df['-logpval'] = -1*np.log10(self.df['Pval']) 
        #self.df.loc[self.df['-logpval']==np.inf,'-logpval'] = 4

    def terms_with_signal(self, max_pval=0.05, min_pval=0,
                          max_WindowSize=1, max_FlankLimit=0,
                          pval_col='Density_Pval'):
        '''
            Return only terms with a pvalue between the 
            specified arguments.

            Parameters
            ----------

            max_pval : float (default: 0.05)
                The maximum pval that a term can have and still be considered
                to have 'signal'
            min_pval : float (default: 0.0)
                The minimum pval that a term can have and still be considered
                to have 'signal'. This is useful to stratify terms into groups,
                for instance moderate signal: min_pval=0.01 (pval_max: 0.05)
            max_WindowSize : int (default: 1)
                Stratify results terms by window_size as specified in the term
                candidate genes function (snp2genes mapping). Default will 
                constrain genes to only those which are truly in the simulated
                term.
            max_FlankLimit : int (default: 0)
                Stratify terms by flanklimit as specified in the term candidate
                genes function (snp2gene mapping). Default will constrain genes
                to only include those which are truly in the simulated term.
            pval_col : enum Density_Pval,Locality_Pval (default: 'Density_Pval')
                The results pval column to use to define 'signal'.
        '''
        # Split the data in a few informative ways
        # We only want the True FCR
        fdr = self.df[(self.df.WindowSize <= max_WindowSize) \
                & (self.df.FlankLimit <= max_FlankLimit)]
        fdr = fdr.set_index(['COB','GO','FCR'],drop=False).sort_index()
        # Split out groups
        terms_with_signal = fdr.query(
            '''FCR==0 and FlankLimit==0 and {pval_col}<{max_pval} and {pval_col}>={min_pval}'''.format(**locals())
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

       
    @staticmethod
    def plot_aggregate_terms(df,filename=None, figsize=(16,8), alphs=0.05,
                             figaxes=None,label=None):
        '''
            Plot the number of terms significant at varying FCR levels.
            THIS IS FOR WHEN YOU HAVE THE SAME TERMS AND YOU VARY THEIR
            FCR OR NOISE LEVEL!!!
        '''
        breakdown = df.reset_index(drop=True).groupby(['COB','FCR'])\
            .apply(lambda x: sum(x.Pval <= alpha)).reset_index()
        # Get a
        cobs = breakdown.COB.unique()
        if figaxes == None:
            fig, axes = plt.subplots(1,len(cobs),figsize=figsize,sharey=True)
        else:
            fig,axes = figaxes
            if len(cobs) != len(axes):
                raise ValueError('Cannot overlay plots with different number of COBs')
        for i,cob in enumerate(cobs):
            data = breakdown[breakdown.COB==cob]
            axes[i].plot(data['FCR'],data[0],label=label)
            axes[i].set_title("{} Terms".format(cob))
            if i == 0:
                axes[i].set_ylabel('Number Significant Terms')
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
               
    '''
        Static and Helper methods
    '''
                
    @staticmethod
    def log10pval(pvals):
        transformed = -1 * np.log10(pvals)
        inf_mask = np.isinf(transformed)
        transformed[inf_mask] = np.max(transformed[np.logical_not(mask)])


