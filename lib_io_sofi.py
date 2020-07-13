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
#from lib_generic_CONST import *

from lib_generic_PLOTT import *
##


## CONSTS
this_lib = 'lib_io_sofi.py/'


sofi_runfile = '-Runfile.json'
sofi_S = '-Sources.dat'
sofi_R = '-Receivers.dat'

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


def Sofi_SR_Create(proj_name, dx, **kwargs):
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
  Now only 1 fpeak, delay and ampl for all sources.
  
  """
  this_func = this_lib + 'SR_Create: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  s_suffix = '-Sources.dat'
  r_suffix = '-Receivers.dat'
  proj_path = Kwarg('proj_path', './', kwargs)
  Ss = Kwarg('Ss', None, kwargs)
  Rs = Kwarg('Rs', None, kwargs)
  f_peak = Kwarg('f_peak', None, kwargs)
  delay = Kwarg('delay', 0.0, kwargs)
  ampl = Kwarg('ampl', 1.0, kwargs)
  
  if not (Ss and Rs):
    raise ValueError('You must specify both SR lists.')
  
  if not f_peak:
    raise ValueError('You must specify central frequency.')
  
  ##
  for Xs, suffix in zip([Ss, Rs], [s_suffix, r_suffix]):
    fname = proj_path + proj_name + suffix
    
    Xs = np.array(Xs) * dx # NOTE
    
    with open(fname, 'w') as f:
    
      for X in Xs: # NOTE: 2ND COLUMN CONTAINS DEPTH COORDINATE!
        f.write(str(X[0]) + ' ' + str(X[2]) + ' ' + str(X[1]))
        if suffix == s_suffix: # EXTRA INFO FOR SOURCES
          f.write(' ' + str(delay) + ' ' + str(f_peak) + ' ' + str(ampl))
          
        f.write('\n')
  
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Sofi_SR_Read(proj_name, **kwargs):
  """
  FIXME: implement return dims, Ss, Rs
  
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
  this_func = this_lib + 'Sofi_SR_Read: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  fw = Kwarg('fw', None, kwargs)
  if not fw:
    fw = float(Sofi_Runfile_Read(proj_name, 'FW'))

  dx = Kwarg('dx', None, kwargs)
  if not dx:
    dx = float(Sofi_Runfile_Read(proj_name, 'DX'))
    
  Ss_Rs = []
  for suffix in ['-Sources.dat', '-Receivers.dat']:
      c = Read_File(proj_name + suffix)
      points = [[float(i[0])-(fw*dx), float(i[2])-(fw*dx), i[1]] for i in c]
      points = [[float(i)/dx for i in j] for j in points] 
      Ss_Rs.append(points)
  Ss, Rs = Ss_Rs

  #print this_func + 'END'
  return Ss, Rs


# -------------------------------------------------------------------------------


def Sofi_Runfile_Read(proj_name, param, **kwargs):
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
  this_func = this_lib + 'Runfile_Read_Sofi: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  
  file_name = proj_name + sofi_runfile
  
  if verbos > 3:
    print(this_func, 'Runfile to read: ', file_name)
  
  content = Read_File(file_name)
  #params = {}
  for record in content:
    #params[record[0]
    #"FW"
    key = record[0][1:-1] # SKIP ""
    try:
      value = record[2][1:-2] # SKIP : AND ""
    except IndexError:
      continue
    
    if key == param:
      return value
      
    
    #print key, value
    #params[key] = record
    ## GET ONLY RECORDS OF CERTAIN FORMAT:
    #if (len(record) == 3) and (record[1] == ':'):
    #  key = record[0]
    #  value = record[2]
    #  params[key] = value
    #  
    ## SPECIAL FORMAT
    #if (record[0] == 'MAX') and (record[1] == 'TIME'):
    #  key = 'ttime'
    #  value = str(float(record[3]) / 1000.) # CONVERT ms TO s
    #  params[key] = value
  
  raise KeyError('Error. Param ' + param + ' not found')
  
  #print this_func + 'END'  
  return 1


# -------------------------------------------------------------------------------
# PLOTTING
# -------------------------------------------------------------------------------


def Sofi_SR_Plot(proj_name, **kwargs):
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
  this_func = this_lib + 'Sofi_SR_Plot: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_generic_PLOTT import Plot_Points
  
  Ss, Rs = Sofi_SR_Read(proj_name)
  
  c_S, s_S, c_R, s_R = 'grey', 50, 'red', 150
  Plot_Points(Ss, plot_type='proj', c=c_S, s=s_S, label='input Ss', plane='XZ')
  Plot_Points(Rs, plot_type='proj', c=c_R, s=s_R, label='input Rs', plane='XZ')
  plt.gca().set_aspect('equal')
  plt.xlim(0,100)
  plt.ylim(100,0)


  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Sofi_Project_Output_Plot(proj_name, **kwargs):
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
  this_func = this_lib + 'Sofi_Project_Output_Plot: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  pars = {'cmap': 'seismic', 'data_type': 'wavelet', 'verbos': 0}
  
  xlim = Kwarg('xlim', None, kwargs)
  if not xlim:
    xlim = float(Sofi_Runfile_Read(proj_name, 'TIME', **kwargs)) # MAX TIME
  
  dt = Kwarg('dt', None, kwargs)
  if not dt:
    dt = float(Sofi_Runfile_Read(proj_name, 'DT', **kwargs)) 

  fig, ax, i = Subplots(2, 3, figsize=[25,15])
  # 1ST ROW
  i = Subplot(fig, i)
  plt.title('p')
  Plot_Gridded_Data_From_File(proj_name+'-Observed_p.su', **pars)
  i = Subplot(fig, i)
  plt.title('curl p')
  Plot_Gridded_Data_From_File(proj_name+'-Observed_curl.su', **pars)
  i = Subplot(fig, i)
  plt.title('div p')
  Plot_Gridded_Data_From_File(proj_name+'-Observed_div.su', **pars)
  # 2ND ROW
  i = Subplot(fig, i)
  plt.title('vx')
  Plot_Gridded_Data_From_File(proj_name+'-Observed_vx.su', **pars)
  i = Subplot(fig, i)
  plt.title('vy')
  Plot_Gridded_Data_From_File(proj_name+'-Observed_vz.su', **pars)
  i = Subplot(fig, i)
  plt.title('vz (vertical)')
  Plot_Gridded_Data_From_File(proj_name+'-Observed_vy.su', **pars)

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------

