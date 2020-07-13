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
from lib_generic import Kwarg
from lib_generic_CONST import *
##


## CONSTS
this_lib = 'lib_math_physics/'

## FUNCTIONS & CLASSES


# -------------------------------------------------------------------------------


class New_Class:
  """
  
  
  """
  
  def __init__(self, path, name="Not_defined", info="Missing description,"):
    self.path = path
    self.name = name
    self.info = info


# -------------------------------------------------------------------------------


class Field:
  """
  
  
  """
  
  def __init__(self, potential, magnitude):
    self.potential = potential
    self.magnitude = magnitude


# -------------------------------------------------------------------------------
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


def V_Coulomb(x, y):
  """
  
  """
  r = np.sqrt(x**2 + y**2)
  #print r
  
  #if r < epsi:
    #V = 0.0
  #else:
  V = 1. / (r + epsi)
  
  return V


# -------------------------------------------------------------------------------


def V_Yukawa(x, y): 
  V = lambda x, y : -exp(-np.sqrt(x**2 + y**2)) / np.sqrt(x**2 + y**2)
  V = np.vectorize(V)
  return V(x, y)


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


def Density_Gardner(vp, **kwargs):
  """
  Calculate density based on P-wave velocity 
  using Gardner's relation.
  
  Parameters
  ----------
  vp : float 
    P-wave velocity in m/s.

  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  Returns
  -------
  rho : float 
    Density in g/cm3 !!!
  
  Notes
  -----
  It is used in Fullwave3D.
  
  More on: https://en.wikipedia.org/wiki/Gardner%27s_relation
  
  """
  this_func = this_lib + 'Density_Gardner: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  rho = 0.31 * vp**(0.25)

  #print this_func + 'END'
  return rho


# -------------------------------------------------------------------------------


def Bulk_Modulus(**kwargs):
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
  this_func = this_lib + 'Bulk_Modulus: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')



  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------

