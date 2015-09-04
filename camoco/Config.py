#!/usr/env/python3

import os
import yaml
import pprint
import getpass
import configparser

global cf

default_config = '''--- # YAML Camoco Configuration File
options:
    basedir: ~/.camoco/
    testdir: ~/build/{user}/Camoco/tests/

logging:
    log_level: verbose

test:
    force:
        RefGen:   False
        COB:      False
        Ontology: False
    num:      50
    refgen:   Zm5bFGS
    cob:      NewRoot
    ontology: ZmIonome
    term:     Fe57
    gene:     GRMZM2G000014
    
'''.format(**{'user':getpass.getuser()})

class Level(dict):
    '''
        Ha! Take that config parser! I am accessing
        everything like an object.
    '''
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)

    def __getattr__(self,item):
        if isinstance(self[item],dict):
            return Level(self[item])
        else:
            if 'dir' in item and '~' in self[item]:
                return os.path.expanduser(self[item])
            return self[item]

class Config(object):

    def __init__(self,filename):
        filename = os.path.expanduser(filename)
        self.data = Level(yaml.load(open(filename,'r')))

    def __getattr__(self,item):
        return Level(self.data[item])

    def __getitem__(self,item):
        return Level(self.data[item])

    def __repr__(self):
        return pprint.pformat(self.data)

''' -------------------------------------------------------------------------
        Program Logic
'''

cf_file = os.path.expanduser('~/.camoco.conf')

# Check to see if there is a config file available
if not os.path.isfile(cf_file):
    with open(cf_file, 'w') as CF:
        print(default_config, file=CF)

cf = Config(cf_file)
