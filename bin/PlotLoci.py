#!/usr/bin/env python3

import argparse
import sys
import os
import copy

import numpy as np
import scipy as sp
import scipy.stats

import camoco as co
import pandas as pd
import matplotlib.pylab as plt

from camoco.Config import cf
cf.logging.log_level = 'quiet'

def main(args):
    
    cob = co.COB(args.cob)
    cob.plot(filename='Health_{}_raw.png'.format(cob.name),
        raw=True,cluster_method=None
    )
    cob.plot(filename='Health_{}_Nomalized.png'.format(cob.name),
        raw=False, cluster_method='mcl'
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()    
    # Data Set Arguments
    parser.add_argument(
        '--cob',
        help='The camoco network to use.'
    )
    import sys
    from IPython.core import ultratb
    sys.excepthook = ultratb.FormattedTB(mode='Verbose',
                 color_scheme='Linux', call_pdb=1)

    args = parser.parse_args()
    sys.exit(main(args))
