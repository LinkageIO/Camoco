#!/usr/env/python3

import os
import configparser


global cf


cf = configparser.ConfigParser()
cf._interpolation = configparser.ExtendedInterpolation()

cf_file = os.path.expanduser('~/.camoco.conf')

default_config = '''
[options]
basedir = ~/.camoco/
testdir = ~/.camoco/

[logging]
log_level = verbose

[test]
refgen   = Zm5bFGS
cob      = NewRoot
ontology = ZmIonome
term     = Fe57
gene     = GRMZM2G000014
'''

# Check to see if
if not os.path.isfile(cf_file):
    with open(cf_file, 'w') as CF:
        print(default_config,file=CF)

cf.read(os.path.expanduser('~/.camoco.conf'))
