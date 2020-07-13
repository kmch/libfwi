"""
Author: Kajetan Chrapkiewicz, 2018. 
All rights reserved. Ask for permision writing to K.Chrapkiewicz17@imperial.ac.uk.

This library provides framework procedures for...
At the moment fully implemented is:

1. simplest numerical stability, accuracy, resolution, etc. analysis of 
   a FWI forward problem.

"""

## MODULES
import numpy as np
import matplotlib.pyplot as plt

## MY MODULES
from lib_generic import *
from lib_generic_CONST import *
##


## CONSTS
this_lib = 'lib_math_num.py/'

## FUNCTIONS


#--------------------------------------------------------------------------------


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
  return 


# -------------------------------------------------------------------------------
# CALCULUS
# -------------------------------------------------------------------------------


def Derivative(y, **kwargs):
  """
  Calculate derivative.
  
  Parameters
  ----------
  y : array 
    Discrete function to differentiate.
  
  Returns
  -------
  yy : array 
    First derivative of y. Same lenght as y.
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Derivative: '  
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  n = Kwarg('n', 1, kwargs) # 1 IF 1ST DERIVATIVE AND SO ON
  acc = Kwarg('acc', 2, kwargs) # ORDER OF NUMERICAL ACCURACY
  dx = Kwarg('dx', 1, kwargs)
  
  yy = np.zeros(len(y))
  
  if n == 1:
    if acc == 2:
      yy[1:-1] = [(-0.5*y[i-1] + 0.5*y[i+1]) / dx for i in range(1, len(y)-1)]
      
      # 2ND ORDER AS WELL
      yy[0] = ((-3/2.)*y[0] + 2*y[1] - 1/2.*y[2]) / dx # FORWARD FD
      yy[-1] = (3/2.*y[-1] - 2*y[-2] + 1/2.*y[-3]) / dx # NOTE: BACKWARD FD HAS OPPOSITE SIGNS 
      
      
  #print this_func + 'END'
  return yy


# -------------------------------------------------------------------------------


def Integral(y, **kwargs): # FIXME: EMPTY
  """
  Calculate definite integral.
  
  Parameters
  ----------
  y : array 
    Discrete function to differentiate.
  
  Returns
  -------
  yy : array 
    First derivative of y. Same lenght as y.
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Integral: '  
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  n = Kwarg('n', 1, kwargs) # 1 IF 1ST DERIVATIVE AND SO ON
  acc = Kwarg('acc', 2, kwargs) # ORDER OF NUMERICAL ACCURACY
  dx = Kwarg('dx', 1, kwargs)
  
  yy = np.zeros(len(y))
  
  #if n == 1:
    #if acc == 2:
 
      
      
  #print this_func + 'END'
  return yy


# -------------------------------------------------------------------------------
# SCHEME'S QC
# -------------------------------------------------------------------------------


def Check_Stability(logfile_object, dx, dt, v_max, kernel):
  """
  Check stability of the scheme (kernel).
  
  Parameters
  ----------
  logfile_object : file object 
    Already-opened file object 
    (return by open(file_name, ...))
    to write log into it.  
  dx : float 
    Size of the spatial grid cell [m].
  dt : float 
    Time step [s].
  v_max : float
    Max. velocity of the model.
  kernel : str 
    Type of kernel. Current capabilities:
     - 'low' (fullwave3D)
     - 'high' (fullwave3D)
  
  Returns
  -------
  C : float 
    Courant's number used later in 
    CFL criterion C < C_max where 
    C_max is scheme-specific.
  
  Notes
  -----
  The smaller the more stable.
  
  If C > 1 the fastest wave in the model will 
  be covering more than 1 grid cell per time step.
  
  FIXME: Understand why it leads to instability.
  
  """  
  this_func = this_lib + 'Check_Stability: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  # CALCULATE C
  C = Courant_Number(dx, dt, v_max)
  
  # CHOOSE MAX C FOR A GIVEN SCHEME
  if kernel == 'high':
    C_max = 0.38
  elif kernel == 'low':
    C_max = 0.5
  else:
    eprint(this_func + 'Error. Unknown kernel: ' + kernel + '\n')
    quit()
    
  # ADVICE FOR A USER WHAT TO CHANGE TO MAKE THE SCHEME STABLE
  if C > C_max:
    eprint(this_func + 'C_max = ' + str(C_max) + ', C = ' + str(C) + '\n')
    eprint(this_func + 'Error! Unstable for ' + kernel + ' kernel.\n')
    eprint(this_func + 'Possible solutions which will work independently: \n')
    eprint(this_func + '1. Decrease v_max of the model below: ' + str(C_max * dx / dt) + ' m/s\n')
    eprint(this_func + '2. Increase grid cell above: ' + str(dt * v_max / C_max) + ' m\n')
    eprint(this_func + '3. Decrease time step below: ' + str(C_max * dx / v_max) + ' s\n')
    quit()
  
  text = this_func + 'Courant number, C: ' + str(C) + '\n'
  print(text)
  logfile_object.write(text)
  
  #print this_func + 'END'
  return C


# -------------------------------------------------------------------------------


def Courant_Number(dx, dt, v_max):
  """
  Calculate the Courant's number describing
  stability of the scheme.
  
  Parameters
  ----------
  dx : float 
    Size of the spatial grid cell [m].
  dt : float 
    Time step [s].
  v_max : float
    Max. velocity of the model.
  
  Returns
  -------
  C : float 
    Courant's number used later in 
    CFL criterion C < C_max where 
    C_max is scheme-specific.
  
  Notes
  -----
  The smaller the more stable.
  
  If C > 1 the fastest wave in the model will 
  be covering more than 1 grid cell per time step.
  
  FIXME: Understand why it leads to instability.
  
  """
  this_func = this_lib + 'Courant_Number: '  
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  v_grid = dx / float(dt) # VEL. OF 1 GRID-CELL PER TIME-STEP IN PHYSIC. UNITS
  C = v_max / float(v_grid) 

  #print this_func + 'END'
  return C


# -------------------------------------------------------------------------------


def Check_Accuracy(logfile_object, dx, v_min, f_max, kernel):
  """
  Check accuracy of the scheme (kernel).
  
  Parameters
  ----------
  logfile_object : file object 
    Already-opened file object 
    (return by open(file_name, ...))
    to write log into it.  
  dx : float 
    Size of the spatial grid cell [m].
  v_min : float 
    Min. velocity of the model [m/s].
  f_max : float
    Max. frequency present in the source function.
  kernel : str 
    Type of kernel. Current capabilities:
     - 'low' (fullwave3D)
     - 'high' (fullwave3D)
  
  Returns
  -------
  f_allowed : float 
    Max. allowed frequency.
  
  Notes
  -----
  
  """    
  this_func = this_lib + 'Check_Accuracy: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  print((this_func + 'f_max: ' + str(f_max) + ' Hz'))
  
  ricker_fmax2fpeak_ratio = 2. # 1.6 IN Lombard et al. 2004 (TOO SMALL FOR FW3D?)
  
  # CHOOSE MIN NODES PER WAVELENGTH FOR A GIVEN SCHEME
  if kernel == 'high':
    min_nodes_per_wavelen = 3.6
  elif kernel == 'low':
    min_nodes_per_wavelen = 5
  else:
    eprint('Error! Unknown kernel: ' + kernel + '\n')
    quit()
  
  f_allowed = v_min / (min_nodes_per_wavelen * dx)
  
  if f_max > f_allowed:
    eprint(this_func + 'f_allowed = ' + str(f_allowed) + ', f_max = ' + str(f_max) + '\n')
    eprint(this_func + 'Error! Inaccurate for ' + kernel + ' kernel\n')
    eprint(this_func + 'Possible solutions which will work independently: \n')
    nf_max = v_min / (min_nodes_per_wavelen * dx)
    eprint(this_func + '1. Decrease f_max below: ' + str(nf_max) + ' Hz\n')
    eprint(this_func + '  (Decrease f_peak below: ' + str(nf_max / ricker_fmax2fpeak_ratio) + ' Hz)\n')
    eprint(this_func + '2. Decrease grid cell below: ' + str(v_min / (min_nodes_per_wavelen * f_max)) + ' m\n')
    eprint(this_func + '3. Increase v_min of the model above: ' + str(min_nodes_per_wavelen * f_max * dx) + ' m/s\n')
    quit()
  
  shortest_wavelength = v_min / f_max
  
  nodes_per_shortest_wavelength = shortest_wavelength / dx
  text = this_func + 'No. of nodes per (shortest) wavelength: ' + str(nodes_per_shortest_wavelength) + '\n'
  print(text)
  logfile_object.write(text)
  
  #print this_func + 'END'
  return f_allowed


# -------------------------------------------------------------------------------


def Max_Resolution(logfile_object, dx, vel, f_max):
  """
  Calculate max. resolution (size of the smallest resolvable 
  structure) of the model.
  
  Parameters
  ----------
  logfile_object : file object
    Already-opened file object 
    (return by open(file_name, ...))
    to write log into it.  
  dx : float 
    Size of the spatial grid cell [m].
  vel : float / 3D array  
    Velocity of the model. It's checked whether
    it is a single number or gridded data structure.
  f_max : float
    Max. frequency present in the source function.
  
  Returns
  -------
  resol_max : float 
    Max resolution.
  
  Notes
  -----
  
  """    
  this_func = this_lib + 'Max_Resolution: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  resolution_as_fraction_of_wavelength = 0.5 # MAYBE 0.25 (SEE Operto & Virieux)
  
  # CHECK IF VELOCITY IS A SINGLE NUMBER OF A GRIDDED MODEL
  gridded_vel = False
  
  try:
    float(vel)
    v_min = vel
  except TypeError:
    gridded_vel = True
  
  # BIT COMMON FOR BOTH VEL. DATA TYPES
  shortest_wavelength = v_min / f_max
  resol_max = shortest_wavelength * resolution_as_fraction_of_wavelength
  text = this_func + 'Max resolution (assuming v = v_min): ' + str(resol_max) + ' m, or ' + str(resol_max / dx) + ' nodes' + '\n'
  print(text)
  logfile_object.write(text)
  
  # BIT ONLY FOR VEL. MODEL
  if gridded_vel:
    x_mid = int(len(truevp_model) / 2)
    y_mid = int(len(truevp_model[0]) / 2)
    v_bottom = truevp_model[x_mid][y_mid][-1]
    resol_bottom = v_bottom / f_max * resolution_as_fraction_of_wavelength
    text = this_func + 'Resolution at the bottom of the model (center (x, y)): ' + str(resol_bottom) + ' m, or ' + str(resol_bottom / dx) + ' nodes' + '\n'
    print(text)
    logfile_object.write(text)
  
  #print this_func + 'END'
  return resol_max


# -------------------------------------------------------------------------------


def Propagation_Dists(logfile_object, dx, nx1, nx2, nx3, v_min, v_max, f_min, f_max):
  """
  Calculate propagation distance across 
  the model covered by the fastest waves.
  
  Parameters
  ----------
  logfile_object : file object 
    Already-opened file object 
    (return by open(file_name, ...))
    to write log into it.  
  dx : float 
    Size of the spatial grid cell [m].
  nx1 : int 
    In-line dim of the model. 
  nx2 : int 
    Cross-line dim of the model.  
  nx3 : int 
    Depth dim of the model.  
  f_min : float
    Min. frequency present in the source function.
  f_max : float
    Max. frequency present in the source function.
  v_min : float
    Min. vel. of the true model.
  v_max : float
    Max. vel. of the true model.   
  
  Returns
  -------
  dist : float 
    Distance in metres.
    
  Notes
  -----
  
  """    
  this_func = this_lib + 'Propagation_Dists: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')    
  
  ttimes = [10] #s
  shortest_wavelength = v_min / f_max
  
  try:
    longest_wavelength = v_max / f_min
  except ZeroDivisionError:
    eprint(this_func + 'Error. f_min must be > 0.\n')
    quit()
  
  for t in ttimes:
    text = ''
    text += this_func + 'Assuming t = ' + str(t) + ' s, the fastest wave will cover: \n'
    dist = v_max * t
    text += this_func + str(dist) + ' m \n'
    dist_per_shortest_wavelength = v_max * t / shortest_wavelength
    text += this_func + str(dist_per_shortest_wavelength) + ' shortest wavelengths \n'
    dist_per_longest_wavelength = v_max * t / longest_wavelength
    text += this_func + str(dist_per_longest_wavelength) + ' longest wavelengths \n'
    dist_in_nodes = v_min * t / dx
    text += this_func + str(dist_in_nodes) + ' nodes \n'
    dist_as_fraction_nx1 = dist_in_nodes / nx1
    text += this_func + str(dist_as_fraction_nx1) + ' model-sizes in X direction \n'
    dist_as_fraction_nx2 = dist_in_nodes / nx2
    text += this_func + str(dist_as_fraction_nx2) + ' model-sizes in Y direction \n'    
    dist_as_fraction_nx3 = dist_in_nodes / nx3
    text += this_func + str(dist_as_fraction_nx3) + ' model-sizes in Z direction \n'    
    print(text)
    logfile_object.write(text)

  #print this_func + 'END'
  return dist  


# -------------------------------------------------------------------------------
# VARIA
# -------------------------------------------------------------------------------


def Metres2nodes(r, dr):
  """
  
  Example
  -------
  if r = 100 and dr = 50:
   n = 2
  elif r = 101 and dr = 50:
   ERROR.
  
  """
  if float(r / dr).is_integer():
    return int(r / dr)
  else:
    eprint('lib_fwi_fs.py/Metres2nodes: ERROR! ' + str(r) + ' IS NOT A MULTIPLE OF ' + str(dr) + '\n')
    print(r / dr)
    quit()


# -------------------------------------------------------------------------------

