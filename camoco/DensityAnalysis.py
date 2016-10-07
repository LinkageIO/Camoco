import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib.pylab as plt
matplotlib.style.use('ggplot')
matplotlib.use('Agg')

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

class DensityAnalysis(object):
    def __init__(self,dir='.',sep='\t'):
        self.sep=sep
        self._build_data(dir=dir)

    def _build_data(self,dir='./'):
        self.df = pd.concat(
            [pd.read_table(x,sep=self.sep) \
                for x in glob.glob(dir+"/*.tsv") ]
        )

    def pivot(self,value='PValue'):
        data = pd.pivot_table(
            self.df,
            index=['COB','WindowSize','FlankLimit'],
            columns=['Term'],
            values=value
        )
        return data

   
    def plot_pval_heatmap(self,filename=None,pval_cutoff=0.05):
        '''
            Generate a heatmap based on Density PVal
        '''
        cobs = self.df.COB.unique()
        fig,axes = plt.subplots(len(cobs),1,sharex=True)
        fig.set_size_inches(8,11)
        for i,cob in enumerate(cobs):
            data = pd.pivot_table(
                self.df.ix[self.df.COB==cob],
                index=['WindowSize','FlankLimit'],
                columns=['Term'],
                values='PValue'
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
