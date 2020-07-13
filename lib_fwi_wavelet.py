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
from lib_generic import *
from lib_generic_CONST import *
from lib_fwi_project_CONST import dt_default, ns_default, f_max_default
##


## CONSTS
this_lib = 'lib_fwi_wavelet.py/'

## FUNCTIONS


# -------------------------------------------------------------------------------


def New_Func(**kwargs):
  """
  
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'New_Func: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')



  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
# CREATE A NEW WAVELET
# -------------------------------------------------------------------------------


def Wavelet_Create(proj_name, **kwargs):
  """
  
  Parameters
  ----------
  proj_name : str
    Project name assumed to be provided 
    as a first argument of the FWI code
    being a prefix of all FWI input files  
  **kwargs : keyword arguments (optional)  
  
  Returns
  -------
  0
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Wavelet_Create: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  #pars_gen = Kwarg('pars_gen', {}, kwargs)
  dt = Kwarg('dt', dt_default, kwargs) # s AS REQUIRED BY su_rickers.sh
  ns = Kwarg('ns', ns_default, kwargs) 
  #pars_src = Kwarg('pars_src', {}, kwargs)
  f_max = Kwarg('f_max', f_max_default, kwargs)
  center = Kwarg('center', 0, kwargs) # IN samples
  shape = Kwarg('shape', 'ricker', kwargs)
  #path = Kwarg('path', './', kwargs)
  
  if shape == 'ricker':
    f_peak = float(f_max) / 2.
    o, e = Bash2('su_ricker.sh ' + proj_name + ' ' + str(f_peak) + ' ' + str(dt) + ' ' + str(ns), **kwargs)
  
  elif shape == 'dirac':
    from lib_io_fullwave import Save_vtr
    s = np.zeros(ns)
    s[center] = 1
    Z = Array1D_Convert_To_3D(s, **kwargs)
    core = proj_name + '-RawSign'
    Save_vtr(core, Z.shape[0], Z.shape[1], Z.shape[2], Z, **kwargs) 
    o, e = Bash2('convert_vtr2sgy.sh ' + core + '.vtr')
    fname = core + '.sgy'
    dt_micros = str(dt * 1e6) # CONVERT seconds -> micro seconds
    o, e = Bash2('segyread tape=' + fname + ' | ' + 
                 'sushw key=dt a=' + dt_micros + ' | ' + 
                 'segyhdrs | ' + 
                 'segywrite tape=' + fname)
    
  
  #print this_func + 'END'
  return 0 #wavelet


# -------------------------------------------------------------------------------
