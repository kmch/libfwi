"""
Author: Kajetan Chrapkiewicz, 2018. 
All rights reserved. Ask for permision writing to K.Chrapkiewicz17@imperial.ac.uk.

This library serves as an interface between Python scripts and Seismic Unix.

"""

## MODULES
import numpy as np
import matplotlib.pyplot as plt

## MY MODULES
from lib_generic import *
from lib_generic_CONST import verbos_func
##


## CONSTS
this_lib = 'lib_io_su.py/'

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


def SU_Supergather_Split_Lines(fname, **kwargs):
  """
  # FIXME
  Split a file containing gathers for multiple lines
  (a supergather) into files containing a single gather 
  for each of the line.
  
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
  this_func = this_lib + 'SU_Supergather_Split_Lines: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  #key_line = Kwarg('key_line', 'ep', kwargs)
  #line_range = Kwarg('line_range', range(1, 61), kwargs)
  #
  #ext = Ext(fname)
  #
  #for line in [26]:#line_range:
  #  #command = 'sgyread.sh ' + fname + ' | suwind key=' + key_line + ' min=' + str(line) + ' max=' + str(line) + ' | sgyhdrs | sgywrite tape=' + fname[ :-(len(ext) + 1)] + '_line' + str(line) + '.' + ext
  #  line = str(line)
  #  file_in = fname
  #  key = key_line
  #  value = line
  #  file_out = fname[ :-(len(ext) + 1)] + '_line' + line + '.' + ext
  #  command = 'su_supergather_split_lines.sh ' + file_in + ' ' + value + ' ' + file_out
  #  o = Bash(command)
  #  print command
  #  ###print o
  #  break
  
  #o = Bash('surange.sh ' + file_name)
  #ntraces = o[20:23] # FIXME: DOUBLE-CHECK IT'S GENERIC

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def SU_Get_Traces_No(file_name, **kwargs):
  """
  Get no. of traces from the SEGY 
  header.
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  ntraces : int 
    No. of traces in the SEGY file.
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'SU_Get_Traces_No: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  
  from os.path import isfile
  
  if not isfile(file_name):
    raise IOError('File ' + file_name + ' not found.')
 
 
  o = Bash('surange.sh ' + file_name)
  #ntraces = o[20:23] # FIXME: DOUBLE-CHECK IT'S GENERIC
  lines = o.splitlines()
  ntraces = None
  for line in lines:
    if line[1] == 'traces:':
      ntraces = line[0]
      break
      
  if not ntraces:
    raise IOError('Failed to read number of traces')
  
  #for line in o:
    #print line
  #try:
  #  ntraces = int(ntraces)
  #except ValueError:
  #  ntraces = o[20:22]
  #  try:
  #    ntraces = int(ntraces)
  #  except ValueError:
  #    ntraces = o[20:21]
  #    try:
  #      ntraces = int(ntraces)
  #    except ValueError:
  #      eprint(this_func + 'Error. Could not read no. of traces from SEGY header.\n')
  #      quit()
  
  #else:
  #  ntraces = 1 # FIXME: INTERIM WORK-AROUND
    
  #try:
  #  o = Bash('surange.sh ' + file_name)
  #  print o
  #  ntraces = o[20:23] # FIXME: DOUBLE-CHECK IT'S GENERIC
  #except OSError:
  #  eprint(this_func + 'Could not read the number of traces. Assumes 1 trace...\n')
  #  ntraces = 1
  #  #eprint(this_func + 'Error. Could not read the number of traces.\n')
  #  #quit()
  
  #print this_func + 'END'
  return ntraces


# -------------------------------------------------------------------------------

