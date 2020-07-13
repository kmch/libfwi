"""
Author: Kajetan Chrapkiewicz, 2018. 
All rights reserved. Ask for permision writing to K.Chrapkiewicz17@imperial.ac.uk.

This library provides framework procedures for...
At the moment fully implemented is:

1. ...

"""

## MODULES
import numpy as np
import matplotlib.pyplot as plt

## MY MODULES
#from lib_generic import *
#from lib_generic_CONST import *
##


## CONSTS
this_lib = 'lib_generic_CLASS.py/'


# --------------------------------------------------------------------------


class File:
  '''
  File.
  
  '''
  def __init__(self, name, path, info='Missing info.'):
    self.name = name
    self.path = path
    self.info = info
    
    
# --------------------------------------------------------------------------

