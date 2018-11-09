import types
import six
import random
import pandas as pd
import logging

from .COB import COB
from .Camoco import Camoco

from functools import wraps 
from collections.abc import Iterable
from itertools import permutations
from minus80 import Freezable
from contextlib import contextmanager


def accepts_iterable(fn):
    '''
    This decorator detects when an iterable is passed into the method call
    instead of a single element (first arg only). Then, instead of calling the
    function on the arg (normally) the method is applied to each element of the
    iterable. This essentially turns fn(arg) into [f(x) for x in arg]. 
    '''
    @wraps(fn)
    def wrapped(self,arg,*args,**kwargs):
        # 
        if (isinstance(arg,Iterable) and not \
            isinstance(arg, six.string_types)):
            return [fn(self,x,*args,**kwargs) for x in arg]      
        else:
            fn(self,arg,*args,**kwargs)
    return wrapped

    def __init__(self,name,networks):
        self.networks = set()

        # Add all the networks
        for n in networks:
            self.add_network(n)

    def add_network(self,net):
        '''
            Add a network (COB) to the 
            NetComp object.
        '''
        if isinstance(net,str):
            net = COB(net)
        if not isinstance(net,COB):
            raise ValueError(f'a valid network must be provided')
        self.networks.add(net)
