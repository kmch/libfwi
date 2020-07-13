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
this_lib = 'lib_exceptions.py/'

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


class Error(Exception):
  """Base class for exceptions in this module."""
  pass


# -------------------------------------------------------------------------------


class InputError(Error):
  """Exception raised for errors in the input.

  Attributes:
    expression -- input expression in which the error occurred
    message -- explanation of the error
  """

  def __init__(self, expression, message):
    self.expression = expression
    self.message = message


# -------------------------------------------------------------------------------


class ParameterError(Error):
  """Exception raised for errors in the input.

  Attributes:
    expression -- input expression in which the error occurred
    message -- explanation of the error
  """

  def __init__(self, message):
    #self.expression = expression
    self.message = message
    

# -------------------------------------------------------------------------------


class TransitionError(Error):
  """Raised when an operation attempts a state transition that's not
  allowed.

  Attributes:
    previous -- state at beginning of transition
    next -- attempted new state
    message -- explanation of why the specific transition is not allowed
  """

  def __init__(self, previous, next, message):
    self.previous = previous
    self.next = next
    self.message = message


# -------------------------------------------------------------------------------

