#!/usr/env/python3 

import os
import configparser

global cf
cf = configparser.ConfigParser()
cf._interpolation = configparser.ExtendedInterpolation()
cf.read(os.path.expanduser('~/.camoco.conf'))

