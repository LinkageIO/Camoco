#!/usr/bin/python3.4
activate_this = '/home/schaefer/Envs/venv2/bin/activate_this.py'
# This AINT PYTHON 2!!!
#execfile(activate_this,dict(__file__=activate_this))
with open(activate_this) as f:
    code = compile(f.read(), activate_this, 'exec')
    exec(code, dict(__file__=activate_this))

import sys
sys.path.insert(0,'/heap/pyCOB/interfaces/')
sys.path.insert(0,'/heap/pyCOB/src/')

from CamocoWeb import app as application
