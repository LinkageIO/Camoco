import camoco as co
import glob
import pandas as pd
import os
import numpy as np

def read_FDR(glob_path):
    dfs = []
    for x in glob.glob(glob_path):
        df = pd.read_table(x,sep=',')
        net,gwas,win,flank,*junk = os.path.basename(x).replace('.','_').split('_')
        if 'WindowSize' not in df.columns:
            df.insert(0,'WindowSize',win)
        if 'NumFlank' not in df.columns:
            df.insert(0,'NumFlank',flank)
        dfs.append(df)
    df = pd.concat(dfs)
    # I guess we forgot this before
    df.insert(4,'TraitType','Element')
    df.loc[[x.startswith('Log') for x in df.Term],'TraitType'] = 'Log'
    df.loc[[x.startswith('PCA') for x in df.Term],'TraitType'] = 'PCA'
    df.loc[[x.startswith('Trans') for x in df.Term],'TraitType'] = 'Trans'
    return df

def zmax(a):
    if len(a) == 0:
        return 0
    else:
        return np.max(a)

def groupedFDR(df):
    def grouped_agg(x):
        return pd.Series(
            {
                'Tot':  sum(x.numReal),
                'FDR10':zmax(x[x.FDR<=0.1].numReal),
                'FDR35':zmax(x[x.FDR<=0.35].numReal),
                'FDR50':zmax(x[x.FDR<=.5].numReal)
            }
        )
    groups = ['Ontology','COB','WindowSize','NumFlank','TraitType','Term']
    return df.reset_index().groupby(groups).apply(grouped_agg)
