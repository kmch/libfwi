"""
Author: Kajetan Chrapkiewicz, 2018. 
All rights reserved. Ask for permision writing to K.Chrapkiewicz17@imperial.ac.uk.

This library provides framework procedures for...
At the moment fully implemented is:

1. ...

"""

## MODULES
#from pylab import *
import matplotlib # FOR LOG. XTICKS
import matplotlib.pyplot as plt

from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FuncFormatter
import matplotlib.colors as mcolors # FOR OWN COLORMAPS
from matplotlib import ticker
from mpl_toolkits.mplot3d import Axes3D
# COLORMAPS
import matplotlib.cm as cm
import cmocean

from lib_generic import *
from lib_generic_CONST import *
from lib_fwi_generic_CONST import *
from lib_fwi_project_CONST import *
##

## CONSTS
this_lib = 'lib_generic_PLOTT.py/'



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
# TOY DATA
# -------------------------------------------------------------------------------


def Toy_Data_3D(**kwargs):
  """
  In .vtr format.
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  NOTE: Meshgrid is limited to 2D - we can't
  use it here!
  
  """
  this_func = this_lib + 'Toy_Data_3D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from numpy import pi as PI
  
  nx1, nx2, nx3 = Kwarg('dims', dims_default, kwargs)
  vel = Kwarg('vel', 5000, kwargs)
  inhomo = Kwarg('inhomo', 'lateral', kwargs)
  vel_lambda = Kwarg('vel_lambda', 10, kwargs)
  vel_ampl = Kwarg('vel_ampl', 0.01*vel, kwargs)
  
  lamb_x = vel_lambda
  lamb_y = vel_lambda
  lamb_z = vel_lambda
  phi_x = 0
  phi_y = 0
  phi_z = 0
  
  # NOTE: HOMOG. MODEL
  model = np.ones(shape=(nx1, nx2, nx3), dtype=np.float32) * vel
  
  for x in range(nx1):
    value_x = np.sin(phi_x + 2 * PI * x / lamb_x) 
    
    for y in range(nx2):
      value_y = np.sin(phi_y + 2 * PI * y / lamb_y) 
      
      for z in range(nx3):
        value_z = np.sin(phi_z + 2 * PI * z / lamb_z)   
        
        if inhomo == 'homo':
          continue
        
        elif inhomo == 'vertical':
          raise ValueError('Vertical inhom. not yet implemented')
          
        elif inhomo == 'lateral':
          dm = vel_ampl * value_x * value_y
          #model[x]
          #raise ValueError('Lateral inhom. not yet implemented')
        
        elif inhomo == 'full':
          dm = vel_ampl * value_x * value_y * value_z
        
        #elif (inhomo == 'lateral') or (inhomo == 'full'):
        #  if inhomo == 'lateral':
        #    arg = np.sqrt(x**2 + y**2)
        #  
        #  else:
        #    arg = np.sqrt(x**2 + y**2 + z**2)
        #    
        #  model[x][y][z] = model[x][y][z] + vel_ampl * model[x][y][z] * np.sin(arg * 2 * np.pi / vel_lambda)
          
        else:
          raise ValueError('Unknown inhomogeneity type: ' + inhomo + '\n')  

        model[x][y][z] += dm

  #print this_func + 'END'
  return model


# -------------------------------------------------------------------------------


def Toy_Data_2D(mgrid, **kwargs):
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
  this_func = this_lib + 'Toy_Data_2D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_math_generic import Sine_2D
  
  X = Kwarg('X', 'none', kwargs)
  Y = Kwarg('Y', 'none', kwargs)
  
  #if X == 'none':
    #X = range(100)
  #if Y == 'none':
    #Y = range(100)
  #XX, YY = np.meshgrid(X, Y)
  XX, YY = mgrid
  # TRANSPOSE <-> MIMIC MESHGRID (FOR COMPATIBLITY WITH THE REST OF PACKAGE)
  XX = XX.T
  YY = YY.T
  
  #ZZ = Vec(Sine_2D, XX, YY, **kwargs)
  sin2d = lambda x, y : np.sin(x/10)
  ZZ = sin2d(XX, YY)
  ZZ = Array2D_Convert_To_3D(ZZ,**kwargs)
  
  #print this_func + 'END'
  return ZZ #XX, YY, ZZ


# -------------------------------------------------------------------------------


def Toy_Data_1D(**kwargs):
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
  this_func = this_lib + 'Toy_Data_1D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  X = np.array(list(range(1000)))
  Y = Vec(Sine, X, y0=0, A=10, k=100)

  #print this_func + 'END'
  return X, Y


# -------------------------------------------------------------------------------
# GENERIC FUNCTIONS
# -------------------------------------------------------------------------------


def Save(fname, **kwargs):
  """
  
  Notes
  -----
  save_as overwrites fname + suffix.
  
  """
  this_func = this_lib + 'Save: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
    
  save_as = Kwarg('save_as', None, kwargs)
  suffix = Kwarg('suffix', None, kwargs)
  
  if suffix:
    ext = Ext(fname)
    fname = Strip(fname) + '_' + suffix + '.' + ext
    if verbos > 0:
      print(this_func, fname)
  
  if not save_as:
    save_as = fname     
  
  #print 'save_as', save_as
  fname = Save_As(save_as, **kwargs)
   
  #print this_func + 'END'  
  return fname


# -------------------------------------------------------------------------------


def Save_As(fname, **kwargs):
  """
  
  """
  this_func = this_lib + 'Save_As: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
    
  core = Strip(fname)
  plt.title(core)
  fname = core + '.png'
  plt.savefig(fname, format='png')
  plt.close()
  
  if verbos > 0:
    print(this_func, 'Saved the file:', fname)
  
  #print this_func + 'END'  
  return fname


# -------------------------------------------------------------------------------


def Animate(fnames, anim_name, **kwargs):
  """
  
  """
  delay = Kwarg('delay', 50, kwargs)
  loop = Kwarg('loop', 0, kwargs)
  
  fnames_string = ''
  for fname in fnames:
    fnames_string += fname
  
  o, e = Bash2('convert -delay ' + str(delay) + ' -loop ' + str(loop) + ' ' + 
               fnames_string + ' ' + anim_name)
  
  o, e = Bash2('cp ' + anim_name + ' ~/Desktop/')


# -------------------------------------------------------------------------------
# GENERIC WRAPPERS
# -------------------------------------------------------------------------------


def Plot(fname_or_Z, **kwargs):
  """
  Just a framework to plot both files 
  and arrays with the same function.
  
  Parameters
  ----------
  fname_or_Z : see below 
    If string, it is assumed to 
    stand for a file name (incl.
    the path if file is outside './'),
    otherwise it must be an array.
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Plot: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  if isinstance(fname_or_Z, str):
    Plot_File(fname_or_Z, **kwargs)
  
  elif type(fname_or_Z) == type(np.array([])):
    Plot_Data(fname_or_Z, **kwargs)
  
  else:
    raise TypeError('First argument needs to be either ' + 
                    'a file-name or an array!\n')
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Files(fnames, func, **kwargs): # FIXME
  """
  Plot multiple files (e.g. time
  evolution).
  
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
  this_func = this_lib + 'Plot_Files: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  file_min = Kwarg('file_min', 0, kwargs) # START WITH 0 (PYTHONIC)
  file_max = Kwarg('file_max', len(fnames) , kwargs)
  file_step = Kwarg('file_step', 1, kwargs)
  #path = Kwarg('path', './', kwargs)
  #fnames = Get_Files(path, 
  
  i = -1 # TO START WITH 0
  for fname in fnames:
    i += 1
    if (i < file_min) or (i > file_max) or (i % file_step != 0):
      continue # SKIP
    
    func()
    #Plot_Wavefield_Overlain(fname, **kwargs)

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------  


def Plot_File(fname, **kwargs):
  """
  Just a short-name wrapper 
  to have less to type.
  
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
  this_func = this_lib + 'Plot_File: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START') 
  
  Plot_Gridded_Data_From_File(fname, **kwargs)
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
  

def Plot_Data(Z, **kwargs):
  """
  Just a short-name wrapper 
  to have less to type.
  
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
  this_func = this_lib + 'Plot_Data: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  Plot_Gridded_Data(Z, **kwargs)

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------  


def Plot_Gridded_Data_From_File(fname, **kwargs):
  """
  Framework function to plot 3D gridded data 
  from file of any reasonable format.
  # NOTE LEFT FOR COMPATIBILITY
  
  # FIXME: INCLUDE PATH AS PART OF FILE NAME
  OR Proj_path?
  THE FORMER, PROBABLY.
  Proj_path is only for saving location!
  
  Parameters
  ----------
  fname : str 
    File name. It must include extension proceeded by 
    a full stop (e.g. 'file.vtr'). It can include path 
    if needed (e.g. '/home/file.vtr').  
  **kwargs : keyword arguments, optional
      Just passing it down.
  
  Returns
  -------
  0
  
  Notes
  -----
  To plot 2D data (nx2=1) just set yslice=0 (one of the kwargs).
  The same trick works for seismic gathers as well.
  
  """  
  this_func = this_func = this_lib + 'Plot_Gridded_Data_From_File: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from os.path import isfile
  
  ext = fname.split(".")[-1]
  fname_core = fname[ :-(len(ext)+1)]
  kwargs['fname_core'] = fname_core # ADD FOR EXTRA PLOTTING (Plot_Trace)
  proj_path = Kwarg('proj_path', './', kwargs)
  # FIXME: IT ASSUMES proj_name-Rest
  #proj_name = Kwarg('proj_name', Split(fname_core, '-')[0], kwargs) 
  suffix = Kwarg('suffix', '', kwargs)
  save = Kwarg('save', False, kwargs)
  save_as = Kwarg('save_as', None, kwargs)
  
  #ax = Kwarg('ax', None, kwargs)
  #plt.sca(ax)
  
  #if not Exists(fname):
    #eprint(this_func + 'File ' + fname + ' not found.\n')
 
  if ext == 'vtr':
    from lib_fwi_generic_PLOTT import Plot_vtr
    Plot_vtr(fname, **kwargs)
  
  elif ext == 'ttr':
    from lib_fwi_generic_PLOTT import Plot_ttr
    Plot_ttr(fname, **kwargs)

  elif ext == 'su':
    from lib_fwi_generic_PLOTT import Plot_su
    Plot_su(fname, **kwargs)  
  
  elif ext == 'sgy':
    from lib_fwi_generic_PLOTT import Plot_sgy
    Plot_sgy(fname, **kwargs)
  
  else:
    raise ValueError('Unknown extension: ' + ext + '\n')
  
  # SAVE
  #if save and not evo:
  #  if save_as: # if not None
  #    new_fname = save_as
  #  else: # if None
  #    new_fname = fname_core + suffix + '.png'
  #  print this_func, 'Saving a new figure: ', new_fname 
  #  plt.savefig(new_fname)
  #  plt.close()
  
  # RETURN PLOT-HANDLES
  fig, ax = plt.gcf(), plt.gca()
  
  #print this_func + 'END'
  return fig, ax


# -------------------------------------------------------------------------------


def Plot_Gridded_Data(Z, **kwargs): #NOTE ESSENTIAL
  """
  Framework function to plot 3D gridded data 
  of various types. 
  # NOTE LEFT FOR COMPATIBILITY
  
  NOTE!!!
  It assumes to produce only a single plot
  of a single array.
  
  => DOESN'T (!) NEED TO BE MERGED WITH Plot_2D.
  
  Parameters
  ----------
  Z : array / list
    Gridded 3D data structure of (nx1*nx2*nx3) shape.
  **kwargs : keyword arguments, optional except for:
      data_type : str 
        Physical character of the gridded data. 
        At the moment possible types are: 
        'FS' - free surface plots
        'model' - for the model / wavefield (!) plots
      
  Returns
  -------
  0
  
  Notes
  -----
  At the moment only FS needs a special treatment. 
  Wavefield can be plotted in the same way as model.
  
  NOTE: change 'model' -> 'volume', 'FS' -> 'surface'?
  (needs changes in all notebooks).
  
  """    
  this_func = this_lib + 'Plot_Gridded_Data: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_math_signal import Data_Modify
  
  # CHECK KEYWORD ARGUMENTS
  interleave = Kwarg('interleave', [], kwargs)
  
  try:
    data_type = kwargs['data_type']
  except KeyError:
    data_type = 'model'
    if verbos > 2:
      eprint(this_func + "Warning. data_type not specified - assuming 'model'. \n")
  ##
  
  
  # NOTE: NORMALIZE
  Z = Data_Modify(Z, **kwargs) 
  if len(interleave) > 0: # NORMALIZE ALSO 2ND ARRAY
    interleave = Data_Modify(interleave, **kwargs) 
  ##
  
  
  
  # PLOT
  if data_type == 'model':
    from lib_fwi_generic_PLOTT import Plot_Model
    Plot_Model(Z, **kwargs)
    
  elif data_type == 'FS':
    from lib_fwi_generic_PLOTT import Plot_FS
    Plot_FS(Z, **kwargs)
  
  elif data_type == 'wavelet': # FIXME: TEMPORARILY?
    from lib_fwi_generic_PLOTT import Plot_Trace, Plot_Trace_Evolution
    evo = Kwarg('evo', False, kwargs)
    trace = Z[0][0]
    if evo:
      Plot_Trace_Evolution(trace, **kwargs)
    else:
      Plot_Trace(trace, **kwargs)
  
  elif data_type == 'data': # FIXME
    Plot_Wiggles(Z, **kwargs)
  
  else:
    raise ValueError('Unknown data_type: ' + data_type + '\n')
    
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
# 3D DATA VISUALIZATION
# -------------------------------------------------------------------------------


def Plot_3D(**kwargs):
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
  this_func = this_lib + 'Plot_3D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  
  figsize = Kwarg('figsize', figsize_default, kwargs)
  ax = Kwarg('ax', None, kwargs)
  #XX = Kwarg('XX', None, kwargs)
  #YY = Kwarg('YY', None, kwargs)
  surfs = Kwarg('surfs', [], kwargs)
  curves = Kwarg('curves', [], kwargs)
  scatts = Kwarg('scatts', [], kwargs)
  
  if not ax:
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')
  
  # VOXELS 
  # ...
  
  # SURFACES
  for surf in surfs:
    Plot_Surface(surf, ax=ax, **kwargs)
  
  # CURVES 
  # ...
  
  # POINTS
  rainbow = Rainbow(len(scatts))
  for scatt in scatts:
    Plot_Points(scatt, ax=ax, c=next(rainbow), **kwargs)
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Voxels(Z, **kwargs):
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
  this_func = this_lib + 'Plot_Voxels: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')



  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Isosurf(Z, **kwargs):
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
  this_func = this_lib + 'Plot_Isosurf: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')



  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Cube(Z, **kwargs): # BACKEND
  """
  Plot cube with 3 faces being X/Y/Z slices
  respectively. NOTE: works only for inhomogeneous 
  models (required by contourf).
  
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
  this_func = this_lib + 'Plot_Cube: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_fwi_generic import Slice_Data
  
  verbos = Kwarg('verbos', 0, kwargs)
  cmap = Kwarg('cmap', 'inferno', kwargs)
  offset = Kwarg('offset', 0, kwargs)
  minn = Kwarg('minn', np.min(Z), kwargs)
  maxx = Kwarg('maxx', np.max(Z), kwargs)    
  if cmap == 'seismic' or cmap == 'seismic0':
    absmax = max(abs(minn), abs(maxx)) # FIXME  
    minn = -absmax
    maxx = absmax
  
  # OVERWRITES ABOVE
  absmax = Kwarg('absmax', None, kwargs)
  if absmax:
    print('jee')
    minn = -absmax
    maxx = absmax
  else:
    print('yo')
    
  alpha = Kwarg('alpha', 1, kwargs)
  
  nx1, nx2, nx3 = len(Z), len(Z[0]), len(Z[0][0])
  
  # READ ARGS
  xyz = Kwarg('xyz', [nx1/2, nx2/2, nx3/2], kwargs)

  
  faces = []
  for slice_coord, coord_value in zip(['x', 'y', 'z'], xyz):
    if verbos > 2:
      print(this_func, 'slice_coord', slice_coord)
      print(this_func, 'coord_value', coord_value)
    sliced, shots = Slice_Data(Z, slice_coord=slice_coord, coord_value=coord_value, labels=False, **kwargs)
    faces.append(sliced)
  
  xface, yface, zface = faces
  if verbos > 3:
    for face in faces:
      print(this_func, 'face.shape', face.shape)
  
  fig = plt.figure(figsize=[16,8])
  ax = fig.gca(projection='3d')
  levels = np.linspace(np.min(Z),np.max(Z),100)
  
  # Z SLICE
  X = list(range(1, nx1 + 1))
  Y = list(range(1, nx2 + 1))
  X, Y = np.meshgrid(X, Y)  
  try:
    ax.contourf(X, Y, zface.T, zdir='z', offset=offset, levels=levels, cmap=cmap, vmin=minn, vmax=maxx, alpha=alpha)
  except ValueError:
    raise ValueError('Built-in function contourf cannot plot a homogeneous model!') 
  
  # X SLICE
  X = list(range(1, nx2 + 1))
  Y = list(range(1, nx3 + 1))
  X, Y = np.meshgrid(X, Y)
  X = X.T
  Y = Y.T
  ax.contourf(xface, X, Y, zdir='x', offset=nx1+offset, levels=levels, cmap=cmap, vmin=minn, vmax=maxx, alpha=alpha)
  
  # Y SLICE
  X = list(range(1, nx1 + 1))
  Y = list(range(1, nx3 + 1))
  X, Y = np.meshgrid(X, Y)
  X = X.T
  Y = Y.T  
  ax.contourf(X, yface, Y, zdir='y', offset=offset, levels=levels, cmap=cmap, vmin=minn, vmax=maxx, alpha=alpha)
  
  
  #
  ## setting 3D-axis-limits:    
  ax.set_xlim3d(1, nx1)
  ax.set_ylim3d(1, nx2)
  ax.set_zlim3d(1, nx3)
  ax.invert_zaxis()
  ax.set_aspect('equal')
  ax.set_xlabel('x')
  ax.set_ylabel('y')
  ax.set_zlabel('z')

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Surface(Z, **kwargs): # BACKEND
  """
  Framework function to plot 2D surface
  embedded in 3D space in a 3D view.
  
  Parameters
  ----------
  **kwargs : keyword arguments, optional except for:
  
  Returns
  -------
  0
  
  Notes
  -----
  This completely replaced Plot_FS*
  
  """  
  this_func = this_func = this_lib + 'Plot_Surface: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_generic_PLOTT import Set_3D_Window
  from lib_generic_CONST import alpha_default
  
  figsize = Kwarg('figsize', figsize_default, kwargs)
  try: 
    ax = kwargs['ax']
  except KeyError:
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')    
  
  stride = Kwarg('stride', 10, kwargs)
  color = Kwarg('color', None, kwargs)
  cmap = Kwarg('cmap', 'viridis_r', kwargs)
  alpha = Kwarg('alpha', alpha_default, kwargs)
  
  
  # NOTE: CONVERT .vtr -> meshgrid IF NEEDED
  if len(Z.shape) == 3:
    Z = [[j[0] for j in i] for i in Z]
  elif len(Z.shape) == 2:
    pass
  else:
    raise ValueError('Wrong len(Z.shape): ' + str(len(Z.shape))) 
    
  try:
    X = kwargs['X']
    Y = kwargs['Y']
  except KeyError:
    x = list(range(1, len(Z) + 1))
    y = list(range(1, len(Z[0]) + 1))
    X, Y = np.meshgrid(x, y)
  
  Z = np.array(Z)
  Z = Z.T
  
  if not color: # COLORMAP OR COLOR 
    ax.plot_surface(X, Y, Z, cmap=cmap, alpha=alpha, rstride=stride, cstride=stride)  
  else:
    ax.plot_surface(X, Y, Z, color=color, alpha=alpha, rstride=stride, cstride=stride)

  Format_Plot_3D(**kwargs)
  
  #print this_func + 'END'
  return 0
  

# - SLICES ----------------------------------------------------------------------


def Plot_Slice_Slider(**kwargs):
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
  this_func = this_lib + 'Plot_Slice_Slider: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  Plot_Slice_Single(**kwargs)

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Slice_Animation(**kwargs):
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
  File was read by parent function.
  
  """
  this_func = this_lib + 'Plot_Slice_Animation: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  # ------------------------------------------------------------------------------
  # READ KWARGS
  # ------------------------------------------------------------------------------
  try:
    vols = kwargs['vols']
    Z = vols[0]
  except KeyError:
    raise ValueError('You need to provide at least 1 volume to read its shape.')

  prefix = Kwarg('prefix', 'beta', kwargs)
  anim_name = Kwarg('anim_name', 'anim.gif', kwargs)
  #kwargs['minn'] = Kwarg('minn', np.min(Z), kwargs)
  #kwargs['maxx'] = Kwarg('maxx', np.max(Z), kwargs)
  # NOTE: delay IS SET AT THE BOTTOM BASED ON nsteps
  loop = Kwarg('loop', 0, kwargs)
  plot = Kwarg('plot', 1, kwargs) # YOU CAN TURN OFF PLOTTING JUST TO ANIMATE
  anim = Kwarg('anim', 1, kwargs) # YOU CAN TURN OFF ANIMATING JUST TO PLOT
  single_slice = Kwarg('single_slice', False, kwargs)
  slice_coord = Kwarg('slice_coord', None, kwargs)
  if slice_coord:
    slice_coords = [slice_coord]
    del kwargs['slice_coord']
    anim_name = Strip(anim_name) + '_' + slice_coord + '.gif'
  else:
    slice_coords = ['y', 'x', 'z'] # ALL THREE
  
  # NOTE: GLOBAL COLOR BOUNDS
  minn = np.min(Z)
  maxx = np.max(Z)
  
  all_fnames_str = ''
  nsteps_min = 1e7 # JUST A BIG NUMBER
  for slice_coord in slice_coords:  
  #slice_coord = Kwarg('slice_coord', slice_coord_default, kwargs)
    if slice_coord == 'x':
      n = Z.shape[0]
      y = int(Z.shape[1] / 2)
      z = int(Z.shape[2] / 2)
      xyz = ['dummy', y, z]
      dummy = 0
    elif slice_coord == 'y':
      n = Z.shape[1]
      x = int(Z.shape[0] / 2)
      z = int(Z.shape[2] / 2)
      xyz = [x, 'dummy', z]
      dummy = 1    
    elif slice_coord == 'z':
      n = Z.shape[2]
      x = int(Z.shape[0] / 2)
      y = int(Z.shape[1] / 2) 
      xyz = [x, y, 'dummy']
      dummy = 2    
    else:
      raise ValueError('Wrong slice_coord ' + slice_coord)
    
    # SLICING RANGE
    i_min = Kwarg('i_min', 0, kwargs)
    i_step = Kwarg('i_step', 1, kwargs)
    i_max = Kwarg('i_max', n, kwargs) #NOTE: USES n => MUST BE HERE   
    
    # ------------------------------------------------------------------------------
    # SLICE ALL ALONG
    # ------------------------------------------------------------------------------  
    if verbos > 0:
      print(this_func, 'Slicing with i_min, i_max, i_step', i_min, i_max, i_step)
    
    
    fnames_str = ''
    nsteps = 0
    for coord_value in range(i_min, i_max, i_step):
      #kwargs['coord_value'] = coord_value
      
      
      kwargs['cmap'] = 'gist_ncar'#'jet' #'rainbow' #'nipy_spectral' # 'gnuplot2' #'cubehelix' #'gist_earth' #
      if single_slice:
        xyz = None
        kwargs['slice_coord'] = slice_coord
      else:
        xyz = xyz
        xyz[dummy] = coord_value
      
      Plot_Slice(xyz=xyz, minn=minn, maxx=maxx, 
                 alpha=1.,**kwargs)
      
      fname = prefix + '_' + slice_coord + '_' + str(coord_value) + '.png'
      print(this_func, 'Plotting: ', fname)
      fnames_str += fname
      fnames_str += ' '
      plt.savefig(fname, format='png')
      plt.close()
      nsteps += 1
  
    all_fnames_str += fnames_str
    if nsteps < nsteps_min:
      nsteps_min = nsteps
    
  delay = 1000.0 / nsteps_min # VELOCITY EMPIRICALLY CHOSEN
  delay = Kwarg('delay', delay, kwargs) # 1 IS QUITE FAST ACTUALLY  
  
  o, e = Bash2('convert -delay ' + str(delay) + ' -loop ' + str(loop) + 
               ' ' + all_fnames_str + ' ' + anim_name)
  o, e = Bash2('cp ' + anim_name + ' ~/Desktop')
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Slice(**kwargs):
  """
  Framework function selecting 
  whether to plot 1 or 3 slices.
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  fig, ax
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Plot_Slice: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  xyz = Kwarg('xyz', [], kwargs)

  if len(xyz) != 0:
    Plot_Slices_XYZ(**kwargs)
  else:
    figsize = Kwarg('figsize', figsize_default, kwargs)
    fig = plt.figure(figsize=figsize)
    kwargs['fig'] = fig
    kwargs['ax'] = plt.gca()

    Plot_Slice_Single(**kwargs)
  
  #print this_func + 'END'
  return plt.gcf(), plt.gca()


# -------------------------------------------------------------------------------


def Plot_Slices_XYZ(**kwargs):
  """
  Plot 3 slices of volumetric data
  at a given int(point).
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
    
  
  title = Kwarg('title', None, kwargs)
  
  n1, n2, n3 = Z.shape
  x, y, z = Kwarg('xyz', [n1/2, n2/2, n3/2], kwargs)
  x, y, z = [int(i) for i in [x, y, z]]

  
  Returns
  -------
  0
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Plot_Slices_XYZ: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_fwi_generic_PLOTT import Plot_FS
  
  
  title = Kwarg('title', None, kwargs)
  layout = Kwarg('layout', '2+1', kwargs)
  figsize = Kwarg('figsize', figsize_default, kwargs)
  #print 'sizze', figsize
  #fig, ax, i = Subplots(nrow, ncol, figsize=figsize_default)
  fig = plt.figure(figsize=figsize)
  #fig = Kwarg('fig', plt.figure(), kwargs)
  #print 'fig', fig  
  #x, y, z = [int(i/2) for i in vols[0]]
  x, y, z = Kwarg('xyz', [0, 0, 0], kwargs)
  dims = Kwarg('dims', None, kwargs) # TO lim SR PLOTS ETC.
  data_type = Kwarg('data_type', 'model', kwargs)
  box = Kwarg('box', None, kwargs) # FOR  PHYSICAL UNITS VIA extent
  
  if box:
    kwargs['extent'] = [box[2], box[3], box[4], box[5]] # NOTE: FLIPPED
  
  #x, y, z = [int(i) for i in [x, y, z]]
  #
  #pad = Kwarg('pad', None, kwargs)
  #if pad:
  #  padx, pady, padz = pad, pad, pad
  #else:
  #  padx = Kwarg('padx', min(n1/2, 10), kwargs)
  #  pady = Kwarg('pady', min(n2/2, 10), kwargs)
  #  padz = Kwarg('padz', min(n3/2, 10), kwargs)
  #
  #if verbos > 2:
  #  print(this_func, 'xyz', xyz) 
  
  
  if layout == 'horiz':
    nrow = 1
    ncol = 3
    pos1, rowspan1, colspan1 = [0, 0], 1, 1  # NOTE: IT STARTS FROM 0
    pos2, rowspan2, colspan2 = [0, 1], 1, 1
    pos3, rowspan3, colspan3 = [0, 2], 1, 1
  
  elif layout == 'verti':
    nrow = 3
    ncol = 1
    pos1, rowspan1, colspan1 = [0, 0], 1, 1
    pos2, rowspan2, colspan2 = [1, 0], 1, 1
    pos3, rowspan3, colspan3 = [2, 0], 1, 1    
  
  elif layout == '2+1':
    nrow = 2
    ncol = 3
    pos1, rowspan1, colspan1 = [0, 0], 1, 2
    pos2, rowspan2, colspan2 = [1, 0], 1, 3
    pos3, rowspan3, colspan3 = [0, 2], 1, 1    

  else:
    raise ValueError('Wrong layout: ' + layout)
  
  #if title: #FIXME
    #fig.canvas.set_window_title(title)  
  
  if verbos > 1:
    kwargs['cbar'] = 0
    eprint(this_func + 'Switching off cbar. Some issues.')
    
  # ------------------------------------------------------------------------------
  # 1ST SUBPLOT (X-SLICE)
  # ------------------------------------------------------------------------------  
  
  #i = Subplot(fig, i)
  plt.subplot2grid((nrow, ncol), pos1, rowspan=rowspan1, colspan=colspan1, fig=fig)
  if dims:
    kwargs['xlim'] = [0, dims[1]]
    kwargs['ylim'] = [0, dims[2]]
  
  if box:
    kwargs['extent'] = [box[2], box[3], box[5], box[4]]  
  
  
  Plot_Slice_Single(ax=plt.gca(), slice_coord='x', coord_value=x, 
                    xslice_lines=[y, z], yflip=1, **kwargs)
  

  # ------------------------------------------------------------------------------
  # 2ND SUBPLOT (Y-SLICE)
  # ------------------------------------------------------------------------------    
  
  #i = Subplot(fig, i)
  plt.subplot2grid((nrow, ncol), pos2, rowspan=rowspan2, colspan=colspan2, fig=fig)
  if dims:
    kwargs['xlim'] = [0, dims[0]]
    kwargs['ylim'] = [0, dims[2]]  
  
  if box:
    kwargs['extent'] = [box[0], box[1], box[5], box[4]]
    
  Plot_Slice_Single(ax=plt.gca(), slice_coord='y', coord_value=y, 
                    yslice_lines=[x, z], yflip=1, **kwargs)
                

  # ------------------------------------------------------------------------------
  # 3RD SUBPLOT (Z-SLICE) - NOTE: SPECIAL BEHAVIOUR FOR SURFS, SHOW SLICE-LINES
  # ------------------------------------------------------------------------------                
  
  #i = Subplot(fig, i)
  plt.subplot2grid((nrow, ncol), pos3, rowspan=rowspan3, colspan=colspan3, fig=fig)
  if dims:
    kwargs['xlim'] = [0, dims[0]]
    kwargs['ylim'] = [0, dims[1]]
  
  
  kwargs['zslice_lines'] = [x, y] # NOTE!!!!!!!!!
  #xslice_x = range(kwargs['xlim'])
  
  #curves = []
  
  if box:
    kwargs['extent'] = box[ :-2]  
  
  if data_type == 'FS':
    if 'plot_type' in kwargs:
      del kwargs['plot_type']
    Z = kwargs['surfs'][0]
    Plot_FS(Z, ax=plt.gca(), plot_type='map', yflip=0, **kwargs)
    #pass
  else:
    Plot_Slice_Single(ax=plt.gca(), slice_coord='z', coord_value=z, yflip=0, **kwargs)
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Slice_Single(**kwargs):
  """
  Plot multi-layer 2D slice of various data originally 
  embedded in 3D space.
  
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
  this_func = this_lib + 'Plot_Slice_Single: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_fwi_generic import Slice_Points, Slice_Data
  from lib_generic_CONST import slice_coord_default  
  
  
  vols = Kwarg('vols', [], kwargs)
  surfs = Kwarg('surfs', [], kwargs)
  curves = Kwarg('curves', [], kwargs) # NOTE: A BIT TRICKIER AND NOT NEEDED
  #slice_lines = Kwarg('slice_lines', [], kwargs) # FOR 2+1, SHOW XY LINES ON Z-SLICE
  scatts = Kwarg('scatts', [], kwargs)
  slice_coord = Kwarg('slice_coord', slice_coord_default, kwargs)
  interleave = Kwarg('interleave', [], kwargs)
  
  
  # AS FIRST FOR labels TO BE OVERWRITTEN 
  scatts2d = []
  for scatt in scatts:
    Z, labels = Slice_Points(scatt, **kwargs)
    scatts2d.append(Z)
    
  vols2d = []
  for vol in vols:
    #NOTE: BY DEFAULT THEY WILL BE SLICED AT THEIR MIDDLES -
    
    Z, labels = Slice_Data(vol, **kwargs)
    if len(interleave) != 0 :
      Z2, labels = Slice_Data(interleave, **kwargs)
      Z = Array_Interleave(Z, Z2, **kwargs)
      # NOTE: Z2 AND Z ARE SWAPPED TO APPEND WITHOUT A NEW if
      #Z2, Z = Array_Interleave_Split(Z, Z2, **kwargs)
      #kwargs['alpha'] = .5 # OTHERWISE WHITE GAPS WILL COVER THE OTHER
      #vols2d.append(Z2)
      if len(Z) == 0:
        eprint(this_func + 'Got empty array. ' +
               'skipping this plot...\n')
        return 1
    
    vols2d.append(Z)
  
  surfs2d = []
  for surf in surfs:
    if slice_coord == 'z':
      labels = ['No depth-slice available for surfaces', 'in-line node', 'cross-line node']
      if verbos > 2:
        eprint(this_func + 'Skipped slicing surf in Z coord.\n')
    else:
      try:
        del kwargs['data_type']
      except KeyError:
        pass
      Z, labels = Slice_Data(surf, data_type='FS', **kwargs)
      surfs2d.append(Z)
  
  try:
    del kwargs['surfs'] # VOLS->SURFS a.k.a. VOLS2D
  except KeyError:
    pass
  if 'curves' in kwargs:
    del kwargs['curves']
  
  
  # NOTE: NAMING CONFUSION
  surfs = vols2d
  curves += surfs2d
  Plot_2D(surfs=surfs, curves=curves, scatts2d=scatts2d, labels=labels, **kwargs)
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
# 2D DATA VISUALIZATION
# -------------------------------------------------------------------------------


def Plot_2D(**kwargs): #NOTE: ESSENTIAL
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
  this_func = this_lib + 'Plot_2D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_generic_CONST import zorder_back, zorder_front, zorder_mid
  
  
  ax = Kwarg('ax', None, kwargs)
  XX = Kwarg('XX', None, kwargs)
  YY = Kwarg('YY', None, kwargs)
  surfs = Kwarg('surfs', [], kwargs)
  curves = Kwarg('curves', [], kwargs)
  scatts2d = Kwarg('scatts2d', [], kwargs)
  xslice_lines = Kwarg('xslice_lines', None, kwargs)  
  yslice_lines = Kwarg('yslice_lines', None, kwargs)
  zslice_lines = Kwarg('zslice_lines', None, kwargs)
  
  try:
    if not ax:
      fig, ax = plt.subplots(1, 1) # NECESSARY FOR OVERLAYS
    else:
      del kwargs['ax']
  except ValueError:
    raise ValueError('You probably provided ax of subplots, use plt.gca() instead.')

  # SURFACES
  #if len(surfs) == 2:
    #cmaps = ['inferno_r', 'seismic']
    #del kwargs['cmap']
  for i, surf in enumerate(surfs):
    if i == 0:
      kwargs['alpha'] = 1
    if i == 1:
      kwargs['cmap'] = 'seismic' #'Greys'
      kwargs['cbar'] = 0
      kwargs['alpha'] = .3 
    if i > 1:
      raise NotImplementedError('Colormaps for more layers needed.')
    
    Plot_Map(surf, ax=ax, format_plot=0, zorder=zorder_back, **kwargs)

  # LINES
  # NOTE: RAINBOW TOO?
  for curve in curves: # NOTE: ACTUALLY IT MUST BE y(x) NOT ANY CURVE
    print(curve)
    # ALLOW FOR E.G. VERTICAL CURVES
    if len(curve) == 2:
      eprint('len(curve) == 2\n')
      X = curve[0]
      Y = curve[1]
    else:
      X = None
      Y = curve
      
    Plot_1D(Y, X=X, ax=ax, format_plot=0, zorder=zorder_mid, **kwargs)
  
  # NOTE!!!!!!!!!!!!!
  if len(surfs) > 0:
    nx = surfs[0].shape[0]
    ny = surfs[0].shape[1]  
    
    if xslice_lines:
      xline_x = list(range(nx))
      xline_y = np.zeros(nx) + xslice_lines[1]
      plt.plot(xline_x, xline_y, lw=3, ls=':', c='k') 
    
    if yslice_lines:    
      xline_x = list(range(nx))
      xline_y = np.zeros(nx) + yslice_lines[1]
      plt.plot(xline_x, xline_y, lw=3, ls=':', c='k')    
    
    if zslice_lines:
      xline_x = list(range(nx))
      xline_y = np.zeros(nx) + zslice_lines[1]
      plt.plot(xline_x, xline_y, lw=3, ls=':', c='k')
      
      yline_x = np.zeros(ny) + zslice_lines[0] # VERTICAL
      yline_y = list(range(ny))
      plt.plot(yline_x, yline_y, lw=3, ls=':', c='k')
    
  
  # POINTS
  colors = Rainbow(len(scatts2d), **kwargs)
  for scatt in scatts2d:
    if len(scatts2d) > 1:
      kwargs['marker_color'] = next(colors)    
    Scatter_2D(scatt, ax=ax, format_plot=0, zorder=zorder_front, **kwargs)
  
  Format_Plot_2D(**kwargs)
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Map(Z, **kwargs): # BACKEND
  """
  Framework function to plot 2D projection 
  of a 2D surface embedded in 3D space.
  
  Parameters
  ----------
  **kwargs : keyword arguments, optional except for:
  
  Returns
  -------
  0
  
  Notes
  -----
  
  """  
  this_func = this_func = this_lib + 'Plot_Map: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from pylab import pcolormesh
  
  # CHECK KEYWORD ARGUMENTS
  ax = Kwarg('ax', plt.gca(), kwargs)
  format_plot = Kwarg('format_plot', True, kwargs)
  aspect = Kwarg('aspect', 'equal', kwargs)
  cbar = Kwarg('cbar', 1, kwargs)
  elev = Kwarg('elev', 45, kwargs)
  azim = Kwarg('azim', 45, kwargs)
  stride = Kwarg('stride', 10, kwargs)
  color = Kwarg('color', 'dodgerblue', kwargs)
  cmap = Kwarg('cmap', 'viridis_r', kwargs)
  alpha = Kwarg('alpha', .5, kwargs)
  zorder = Kwarg('zorder', 0, kwargs) # BACK/FRONT POSITION ON CANVAS
  minn = Kwarg('minn', None, kwargs)
  maxx = Kwarg('maxx', None, kwargs)  
  shade = Kwarg('shade', False, kwargs)
  extent = Kwarg('extent', None, kwargs) #[1, len(Z)+1, 1, len(Z[0])+1], kwargs)
  xaxis = Kwarg('xaxis', None, kwargs) # E.G. FOR ASSIGNING CHANNEL NUMBERS
  y_unit = Kwarg('y_unit', None, kwargs) # FOR RESCALING yaxis (CONVERT TIME)
  interp = Kwarg('interp', None, kwargs)
  
  if not interp:
    if verbos > 1:
      print((this_func + "imshow-interpolation disabled. To switch it on, " + 
                       "choose interp=['bilinear','bicubic','sinc', 'spline16'," +
                       "'spline36', 'hanning', 'hamming', 'hermite', 'kaiser'," + 
                       "'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell'," + 
                       "'sinc','lanczos']"))
  
  Z = np.array(Z)
  Z = Z.T # IMPLIED BY THE SHAPE OF THE meshgrid's OUTPUT
  #pcolormesh(X, Y, Z, cmap=cmap, alpha=alpha)
  
  if shade:
    from matplotlib.colors import LightSource
    ls = LightSource(azdeg=315, altdeg=45)
    vert_exag = 1 # VERTICAL EXAGGERATION
    Z = ls.hillshade(Z, vert_exag=vert_exag)
  
  imshow_kwargs = {'cmap': cmap,
                   'alpha': alpha,
                   'zorder': zorder,
                   'interpolation': interp}
  
  if xaxis:
    #print Z.shape[1]
    #print len(xaxis)
    
    xticks = list(range(Z.shape[1]))
    xlabls = xaxis
    
    step = 10
    xticks = xticks[::step]
    xlabls = xlabls[::step]
    
    plt.gca().set_xticks(xticks)
    plt.gca().set_xticklabels(xlabls)
    #xextent = 
    #yextent = [0, Z.shape[1]]
  
  if y_unit:
    xmax = Z.shape[1]
    ymax = Z.shape[0] * y_unit
    extent = [0, xmax, ymax, 0]
    imshow_kwargs['extent'] = extent
    
    
  
  
  if extent:
    imshow_kwargs['extent'] = extent
  
  # CORRECT SHIFT OF THE DIVERGING COLORMAP
  if cmap == 'seismic' or cmap == 'cmo.topo_r' or cmap == 'cmo.topo': # OR OTHER DIRVERGING
    from math import copysign # copysign(x, y) - Return x with the sign of y
    if (not minn) and (not maxx):
      minn = np.min(Z)
      maxx = np.max(Z)
      if minn == 0.0:
        minn = -0.0 # SIGNED ZERO (PLATFORM DEPENDENT) - OTHERWISE WRONG BEHAVI
            
      if abs(minn) > abs(maxx):
        a = abs(minn)
        #signum_a = copysign(1, minn)
      else:
        a = abs(maxx)
        #signum_a = copysign(1, maxx)
      #print a
      #minn = copysign(a, minn)
      #maxx = copysign(a, maxx)
      #print minn, maxx
      
      # IT IS AS SIMPLE AS THAT: FOR COLORMAPS
      # DIVERGING AROUND ZERO WE NEED OPPOSITE SIGNS
      minn = -a 
      maxx = a 
  
  # OVERWRITES ABOVE
  absmax = Kwarg('absmax', None, kwargs)
  if absmax:
    minn = -absmax
    maxx = absmax
  
  if minn and maxx:
    imshow_kwargs['vmin'] = minn 
    imshow_kwargs['vmax'] = maxx
    
  plt.imshow(Z, **imshow_kwargs)
  
  plt.gca().invert_yaxis() # NOTE
  
  if format_plot:
    Format_Plot_2D(**kwargs)
 
  if cbar:
    #plt.colorbar(pad=-0.2,  shrink=0.9, aspect=30)
    plt.colorbar()
    #Colorbar() 
    #plt.gcf().colorbar(im, ax=ax)
  
    
  #print this_func + 'END'
  return 0
  

# -------------------------------------------------------------------------------


def Scatter_2D(XY, **kwargs):
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
  this_func = this_lib + 'Scatter_2D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  
  ax = Kwarg('ax', plt.gca(), kwargs)
  marker = Kwarg('marker', 'o', kwargs)
  marker_size = Kwarg('marker_size', 50, kwargs) # DEFINE BY ITS AREA!
  marker_color = Kwarg('marker_color', 'red', kwargs)
  marker_edge = Kwarg('marker_edge', 'k', kwargs)
  
  zorder = Kwarg('zorder', 0, kwargs)
  
  
  X = [i[0] for i in XY]
  Y = [i[1] for i in XY]
  
  ax.scatter(X, Y, marker=marker, s=marker_size, c=marker_color, 
             edgecolors=marker_edge, zorder=zorder)
  
  return 0
  
  
# -------------------------------------------------------------------------------


def Plot_Wiggles(Z, **kwargs):
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
  It should be merged with Plot_Trace.
  
  Taking the spectrum should be separated from Plot_Trace 
  or called from external.
  
  """
  this_func = this_lib + 'Plot_Wiggles: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  #print Z.shape
  gap = Kwarg('gap', 10, kwargs) # GAP BETWEEN TRACE
  #from lib_fwi_generic_PLOTT import Plot_Trace
  
  t = np.arange(Z.shape[-1])
  
  for i, trace in enumerate(Z):
    trace = trace[0]
    trace += i * gap
    zero_axis = np.ones(len(t)) * i * gap
    Plot_2_Series(t, zero_axis, trace, orient='verti',
                  c1='k', c2='w', c_line='k', lw=.1, **kwargs)
    
    #Plot_Trace(trace)
    
  #plt.gca().invert_yaxis() # DISABLED SINCE IT IS FLIPPED BY ANOTHER FUNCTION
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
# 1D DATA VISUALIZATION
# -------------------------------------------------------------------------------

   
def Plot_2_Series(x, y1, y2, **kwargs): # BACKEND
  """
  Plot 2 time-series with red/blue filling between them.
  It can deal with both horizontal and vertical
  orientations.
  
  Parameters
  ----------

  
  Returns
  -------
  0
  
  Notes
  -----
  It assumes that first samples of both series
  correspond to the same x.
  
  """
  this_func = this_lib + 'Plot_2_Series: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_generic import Equalize_Series_Lengths
  
  # READ ARGS
  xlim = Kwarg('xlim', None, kwargs)
  ylim = Kwarg('ylim', None, kwargs)
  c1 = Kwarg('c1', 'r', kwargs)
  c2 = Kwarg('c2', 'b', kwargs)
  alpha = Kwarg('alpha', 1, kwargs)
  lw = Kwarg('lw', 1, kwargs)
  interpolate = Kwarg('interpolate', True, kwargs)
  xlabel = Kwarg('xlabel', 'X', kwargs)
  ylabel = Kwarg('ylabel', 'Y', kwargs)
  l1 = Kwarg('l1', None, kwargs)
  l2 = Kwarg('l2', None, kwargs)
  orient = Kwarg('orient', 'horiz', kwargs)
  c_line = Kwarg('c_line', 'c2c1', kwargs)
  
  
  # FIXME: ACTUALLY https://matplotlib.org/gallery/lines_bars_and_markers/multicolored_line.html
  if c_line == 'c2c1':
    c_line1 = c2
    c_line2 = c1
  else:
    c_line1 = c_line #'grey'
    c_line2 = c_line #'grey'
    
  
  # PREPARE LISTS
  if (len(x) != len(y1)) or (len(x) != len(y2)):
    eprint(this_func + 'Error. All x, y1, y2 arrays must have the same length. Consider using Equalize_Series_Lengths before')
    quit()
  
  y1 = np.array(y1) # OTHERWISE 'dimensions are inconsistent'
  y2 = np.array(y2) ##
  
  # PLOT
  ax = plt.gca()  

  if orient == 'horiz':
    ax.plot(x, y1, color=c_line1, lw=lw, label=l1)
    ax.plot(x, y2, color=c_line2, lw=lw, label=l2)    
    ax.fill_between(x, y1, y2, where=y2 >= y1, facecolor=c1, interpolate=interpolate, alpha=alpha)
    ax.fill_between(x, y1, y2, where=y2 <= y1, facecolor=c2, interpolate=interpolate, alpha=alpha)
  elif orient == 'verti':
    ax.plot(y1, x, color=c_line1, lw=lw, label=l1)
    ax.plot(y2, x, color=c_line2, lw=lw, label=l2)    
    ax.fill_betweenx(x, y1, y2, where=y2 >= y1, facecolor=c1, interpolate=interpolate, alpha=alpha)
    ax.fill_betweenx(x, y1, y2, where=y2 <= y1, facecolor=c2, interpolate=interpolate, alpha=alpha)
  else:
    raise ValueError('Wrong orient: ' + orient)
  
  # FORMAT
  #plt.xlim(xlim)
  #plt.ylim(ylim)
  #plt.xlabel(xlabel)
  #plt.ylabel(ylabel)
  #plt.grid()
  
  if l1 and l2:    
    plt.legend(loc='upper right', frameon=False, prop={'size' : 15})
  
  Format_Plot_2D(**kwargs)
  
  #print this_func + 'END'
  return 0
  

# -------------------------------------------------------------------------------


def Plot_1D(Y, **kwargs):
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
  this_func = this_lib + 'Plot_1D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  
  ax = Kwarg('ax', plt.gca(), kwargs)
  X = Kwarg('X', None, kwargs)
  lw = Kwarg('lw', 3, kwargs)
  #curves = Kwarg('curves', [], kwargs)
  scatts = Kwarg('scatts', [], kwargs)
  zorder = Kwarg('zorder', 0, kwargs)
  extent = Kwarg('extent', None, kwargs)
  
  if extent:
    x1, x2, y1, y2 = extent # NOTE: y1, y2 ARE DUMMY
    X = np.linspace(x1, x2, len(Y))
  else:
    X = list(range(len(Y)))
    
  ax.plot(X, Y, zorder=zorder, lw=lw)
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
# VARIA (~BILINEAR?)
# -------------------------------------------------------------------------------


def Plot_Tile_3D(Q, **kwargs):
  """
  
  Parameters
  ----------
  Q : list 
    List of x1, x2, y1, y2 = 
    [Q[0][0], Q[2][0], Q[0][1], Q[2][1]]
  **kwargs : keyword arguments, optional
  
  Returns
  -------
  
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Plot_Tile_3D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  grid_x, grid_y, grid_z = Prepare_Gridded_Z('bilinear', Q)  
  Plot_Surface(grid_z, X=grid_x, Y=grid_y, **kwargs)
  plt.figure()
  Plot_Map(grid_z, X=grid_x, Y=grid_y, **kwargs)
  
  #print this_func + 'END'
  return grid_x, grid_y, grid_z


# -------------------------------------------------------------------------------


def Plot_Square(x1, x2, y1, y2, **kwargs): # BACKEND
  """
  Plot a 2D box (square) given its 2 vertices
  (along the diagonal).
  
  Parameters
  ----------
  x1 : float 
    Min. value of X-coord.
  x2 : float 
    Max. value of X-coord.
  y1 : float 
    Min. value of Y-coord.
  y2 : float 
    Max. value of Y-coord.
  **kwargs : keyword arguments (optional)
    Current capabilities:
    - c - color
    - ls - line style 
    - lw - line width
    - alpha - opacity
  
  Returns
  -------
  0
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Plot_Square: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  # CHECK KEYWORD ARGUMENTS
  try:
    c = kwargs['c']
  except KeyError:
    c = 'b'
  
  try:
    ls = kwargs['ls']
  except KeyError:
    ls = '.-'
  
  try:
    lw = kwargs['lw']
  except KeyError:
    lw = '3'
    
  try:
    alpha = kwargs['alpha']
  except KeyError:
    alpha = 0.5    
  
  # PLOT
  plt.plot([x1, x2], [y1, y1], ls, c=c, lw=lw, alpha=alpha)
  plt.plot([x1, x2], [y2, y2], ls, c=c, lw=lw, alpha=alpha)
  plt.plot([x1, x1], [y1, y2], ls, c=c, lw=lw, alpha=alpha)
  plt.plot([x2, x2], [y1, y2], ls, c=c, lw=lw, alpha=alpha)
  
  plt.gca().set_aspect('equal')
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Prepare_Gridded_Z(interp, **kwargs):
  """
  FIXME: FOR BILINEAR?
  
  Parameters
  ----------
  
  
  Returns
  -------
  
  
  Notes
  -----
  # FIXME: WRONG SHAPE?? OF grid_z??
  
  """
  this_func = this_lib + 'Prepare_Gridded_Z: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  from scipy.interpolate import griddata
  
  try:
    grid_x = kwargs['X']
    grid_y = kwargs['Y']
    
  except KeyError:
    try:
      Q = kwargs['Q']
    except KeyError:
      eprint(this_func + 'Error. You must provided either X, Y or Q.\n')
      quit()
    x1, x2, y1, y2 = [Q[0][0], Q[2][0], Q[0][1], Q[2][1]] # CAUSE Q22 = Q[2]
    no_samples = 20j # OR STEP (E.G. 0.1) IF WITHOUT j
    grid_x, grid_y = np.mgrid[x1 : x2 : no_samples, y1 : y2 : no_samples] # SIMILAR TO np.meshgrid()
    
  
  xy_points = np.array([np.array([i[0], i[1]]) for i in Q])
  z_values = np.array([i[2] for i in Q])
  
  if interp == 'nearest':
    grid_z = griddata(xy_points, z_values, (grid_x, grid_y), method='nearest')
  
  elif interp == 'linear':
    grid_z = griddata(xy_points, z_values, (grid_x, grid_y), method='linear')
  
  elif interp == 'cubic':
    grid_z = griddata(xy_points, z_values, (grid_x, grid_y), method='cubic')
    
  elif interp == 'bilinear':
    from lib_math_interp import Interpolate_Bilinear
    try:
      Q = kwargs['Q']
    except KeyError:
      eprint(this_func + 'Error. You must provided either X, Y or Q.\n')
      quit()    
    grid_z = []
    for gx, gy in zip(grid_x, grid_y):
      gz = []
      for x, y in zip(gx, gy):    
        z = Interpolate_Bilinear(x, y, Q)
        gz.append(z)
      grid_z.append(np.array(gz))
    grid_z = np.array(grid_z)    
  
  elif interp == 'dist':
    from lib_math_interp import Interpolate_Bilinear
    from lib_math_generic import Dist
    try:
      Q = kwargs['Q']
    except KeyError:
      eprint(this_func + 'Error. You must provided either X, Y or Q.\n')
      quit()    
    try:
      G = kwargs['G']
    except KeyError:
      eprint(this_func + 'Error. Point to measure the dist from is not specified \n')
      quit()
      
    grid_z = []
    for gx, gy in zip(grid_x, grid_y):
      gz = []
      for x, y in zip(gx, gy):    
        z = Interpolate_Bilinear(x, y, Q)          
        z = Dist(G, np.array([x, y, z]))
        gz.append(z)
      grid_z.append(np.array(gz))
    grid_z = np.array(grid_z)           
  
  else:
    eprint(this_func + 'Error. This interp ("' + interp + '") needs to be added\n')
    quit()
      
  #print this_func + 'END'
  return grid_x, grid_y, grid_z


# -------------------------------------------------------------------------------
# FORMATTING THE CANVAS
# -------------------------------------------------------------------------------


def Set_3D_Window(XX, YY, ZZ, **kwargs): # FIXME: MERGE WITH Format_Plot_3D?
  """
  Create and set parameters of a window 
  for 3D plots.
  
  Parameters
  ----------
  XX : 2D array 
    Array of x-coordinates. Output of meshrid.
  YY : 2D array 
    Array of y-coordinates. Output of meshrid.    
  ZZ : 2D array 
    Array of z-values. Output of function
    def foo(x, y) called as foo(XX, YY).
  **kwargs : keyword arguments (optional)
    fig_size : list
      Figure's dimensions [xdim, ydim].
      Default: [5, 5]
    zflip : bool 
      If true, flip the Z axis. 
      Default: False.
    xlabel : str 
      Label of the X axis.
      Default: 'X [nodes]'
    ylabel : str 
      Label of the Y axis.
      Default: 'Y [nodes]'
    zlabel : str 
      Label of the Z axis.
      Default: 'Z [nodes]'
      
  Returns
  -------
  fig : figure 
    Current figure.
  ax : axes
    Current axes.

  Notes
  -----
  Sets axes' limits as extremes of each array.
  
  FIXME: for some reason adding anything to ax
  invalidates (uniewaznia) zlimits...
  
  """
  this_func = this_lib + 'Set_3D_Window: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from mpl_toolkits.mplot3d import Axes3D
  
  try:
    fig_size = kwargs['fig_size']
  except KeyError:
    fig_size = [5, 5]
  
  try:
    zflip = kwargs['zflip']
  except KeyError:
    zflip = False
  
  try:
    xlabel = kwargs['xlabel']
  except KeyError:
    xlabel = "X [nodes]"

  try:
    ylabel = kwargs['ylabel']
  except KeyError:
    ylabel = "Y [nodes]"  

  try:
    zlabel = kwargs['zlabel']
  except KeyError:
    zlabel = "Z [nodes]"

  # CREATE THE WINDOW
  fig = plt.figure(figsize=fig_size)
  ax = fig.add_subplot(111, projection='3d')
  
  # SCALE THE AXES
  ax.auto_scale_xyz([np.min(XX), np.max(XX)], 
                    [np.min(YY), np.max(YY)], 
                    [np.min(ZZ), np.max(ZZ)]) 
  
  # LABEL THE AXES
  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
  ax.set_zlabel(zlabel)
 
  # FLIP Z AXIS
  if zflip:
    ax.invert_zaxis()

  #print this_func + 'END'
  return fig, ax


# -------------------------------------------------------------------------------


def Format_Plot_3D(**kwargs): # NOTE: MERGE WITH 2D (CLASSES?)
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
  this_func = this_lib + 'Format_Plot_3D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  fig = Kwarg('fig', plt.gcf(), kwargs)
  ax = Kwarg('ax', plt.gca(), kwargs)
  zflip = Kwarg('zflip', False, kwargs)
  azim = Kwarg('azim', 45, kwargs)
  elev = Kwarg('elev', 45, kwargs)
  aspect = Kwarg('aspect', 'auto', kwargs)
  xlim = Kwarg('xlim', None, kwargs)
  ylim = Kwarg('ylim', None, kwargs)  
  
  # ADJUST THE 3D VIEW 
  ax.view_init(azim=azim, elev=elev)  
  plt.gca().set_aspect(aspect)
  
  if xlim:
    plt.xlim(xlim)
  
  if ylim:
    plt.ylim(ylim)  
  
  if zflip:
    ax.invert_zaxis()

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Format_Plot_2D(**kwargs):
  """
  Set parameters of the plot.
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    fig : figure 
      Figure to format.
      Default: plt.gcf().
    ax : axes
      Axes to format.
      Default: plt.gca().
    fig_size : list
      Figure's dimensions [xdim, ydim].
      Default: [5, 5]
    yflip : bool 
      If true, flip the Y axis. 
      Default: False.
    xlabel : str 
      Label of the X axis.
      Default: 'X [nodes]'
    ylabel : str 
      Label of the Y axis.
      Default: 'Y [nodes]'
      
  Returns
  -------
  fig : figure 
    Current figure.
  ax : axes
    Current axes.

  Notes
  -----
  
  """
  this_func = this_lib + 'Format_Plot_2D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from matplotlib.ticker import FormatStrFormatter
  
  
  fig = Kwarg('fig', plt.gcf(), kwargs)
  ax = Kwarg('ax', plt.gca(), kwargs)
  l1 = Kwarg('l1', None, kwargs)
  l2 = Kwarg('l2', None, kwargs)
  yflip = Kwarg('yflip', False, kwargs)
  aspect = Kwarg('aspect', 'auto', kwargs)
  xlim = Kwarg('xlim', None, kwargs)
  ylim = Kwarg('ylim', None, kwargs)
  title = Kwarg('title', None, kwargs)
  xlabel = Kwarg('xlabel', None, kwargs)
  ylabel = Kwarg('ylabel', None, kwargs)
  labels = Kwarg('labels', None, kwargs)
  xticks = Kwarg('xticks', None, kwargs)
  yticks = Kwarg('yticks', None, kwargs)
  grid = Kwarg('grid', None, kwargs)
  
  for lim, axis in zip([plt.xlim(), plt.ylim()], [ax.xaxis, ax.yaxis]):
    minn, maxx = lim
    
    #if (abs(minn) > 1e-2) and (abs(maxx) < 1e3):
      #axis.set_major_formatter(FormatStrFormatter('%5.2f'))
    #else:
      #axis.set_major_formatter(FormatStrFormatter('%.1e'))
      
  plt.setp(ax.xaxis.get_majorticklabels(), rotation=15)
  #ax.tick_params(labelrotation=45) # NOT AVAILABLE
  #ax.xticks(rotation=45)
  #ax.set_yticklabels(rotation=45)
  
  plt.gca().set_aspect(aspect)
  
  if labels and not (title or xlabel or ylabel): # NOTE: IT CAN BE OVERWRITTEN
    title, xlabel, ylabel = labels
  
  if title:
    plt.gca().set_title(title)
  
  if xlabel:  
    plt.xlabel = xlabel

  if ylabel:  
    plt.ylabel = ylabel     
  
  if xlim:
    plt.xlim(xlim)
  
  if ylim:
    plt.ylim(ylim)
  
  if xticks:
    plt.xticks(xticks)

  if yticks:
    plt.yticks(yticks)
  
  if grid:
    plt.grid()
  
  # NOTE: IT OVERRIDES ylim
  if yflip:
    plt.gca().invert_yaxis()
    if verbos > 1:
      print(this_func, 'yflip: True')
  
  #print this_func + 'END'
  return fig, ax
  

# -------------------------------------------------------------------------------


def Subplots(nx, ny, **kwargs):
  """
  Create an axes of 
  nx * ny subplots.
  
  Parameters
  ----------
  n : int 
    Number of colors in the spectrum.
  
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  i : int 
    Iterator over subplots starting 
    from top-left (=> i=1).
    See Subplot(...).
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Subplots: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  figsize = Kwarg('figsize', figsize_default, kwargs) #[14,8]
  title = Kwarg('title', None, kwargs)
  
  fig, ax = plt.subplots(nx, ny, figsize=figsize)
  i = 1 
  
  if title:
    fig.suptitle(title) # MAIN TITLE COMMON FOR ALL SUBPLOTS
  
  plt.tight_layout()
  
  #print this_func + 'END'
  return fig, ax, i


# -------------------------------------------------------------------------------


def Subplot(fig, i, **kwargs):
  """
  Plot a subplot on a previously 
  created subplots-axes.
  
  Parameters
  ----------
  fig : figure 
    Figure's handle.
  i : int 
    As in plt.subplot(nx, ny, i) 
    of the PREVIOUS subplot.
  
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  i : int 
    Updated i.
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Subplot: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  # GET nx, ny INSTEAD OF PASSING IT AS ARGUMENTS
  nx, ny = fig.axes[0].get_subplotspec().get_topmost_subplotspec().get_gridspec().get_geometry()
  
  plt.subplot(nx, ny, i)
  i += 1
  
  #print this_func + 'END'
  return i


# -------------------------------------------------------------------------------


def Rainbow(n, **kwargs):
  """
  Create an iterator for rainbow colors.
  
  Parameters
  ----------
  n : int 
    Number of colors in the spectrum.
  
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  Usage: plot(..., c=next(colors))
  
  """
  this_func = this_lib + 'Rainbow: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  colors = iter(cm.rainbow(np.linspace(0, 1, n)))

  #print this_func + 'END'
  return colors


# -------------------------------------------------------------------------------


def Zoom_To_Point_2D(x, y, **kwargs):
  """
  Set plot limits to zoom in a point (x,y).
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    fig : figure 
      Figure to format.
      Default: plt.gcf().
    ax : axes
      Axes to format.
      Default: plt.gca().    
    r : float 
      Radius of zooming window.
      Default: 5.
    yflip : bool 
      Flip Y axis if true.
      Default: False.
      
  Returns
  -------
  fig : figure 
    Current figure.
  ax : axes
    Current axes.

  Notes
  -----
  
  """
  this_func = this_lib + 'Zoom_To_Point_2D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  # CHECK KEYWORD ARGUMENTS
  try:
    fig = kwargs['fig']
  except KeyError:
    fig = plt.gcf()

  try:
    ax = kwargs['ax']
  except KeyError:
    ax = plt.gca()
    
  try:
    r = kwargs['r']
  except KeyError:
    r = 5

  try:
    yflip = kwargs['yflip']
  except KeyError:
    yflip = False
  
  # ZOOM
  plt.xlim(x-r, x+r)
  if yflip:
    plt.ylim(y-r, y+r)
  else:  
    plt.ylim(y-r, y+r)
  
  #print this_func + 'END'
  return fig, ax


# -------------------------------------------------------------------------------


def Display_Every_Nth_Tick(ax, n):
  """
  Display only every nth tick of both X and Y axis.
  
  Parameters
  ----------
  n : int 
    Step of ticks-display.
  **kwargs : keyword arguments (optional)
    fig : figure 
      Figure to format.
      Default: plt.gcf().
    ax : axes
      Axes to format.
      Default: plt.gca().    
      
  Returns
  -------
  fig : figure 
    Current figure.
  ax : axes
    Current axes.

  Notes
  -----
  
  """
  this_func = this_lib + 'Display_Every_Nth_Tick: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  # CHECK KEYWORD ARGUMENTS
  try:
    fig = kwargs['fig']
  except KeyError:
    fig = plt.gcf()

  try:
    ax = kwargs['ax']
  except KeyError:
    ax = plt.gca()  
  
  # FORMAT
  i = 1
  for tic in ax.xaxis.get_major_ticks():
    if (i % n) != 0:
      tic.label1On = tic.label2On = False
    i += 1
  i = 1
  for tic in ax.yaxis.get_major_ticks():
    if (i % n) != 0:
      tic.label1On = tic.label2On = False
    i += 1
  
  #print this_func + 'END' 
  return fig, ax


# -------------------------------------------------------------------------------


def Force_Aspect(**kwargs):
  """
  Force aspect to actual proportions FIXME: d-c
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    fig : figure 
      Figure to format.
      Default: plt.gcf().
    ax : axes
      Axes to format.
      Default: plt.gca().    
    aspect : float 
      Default: 1.
      
  Returns
  -------
  fig : figure 
    Current figure.
  ax : axes
    Current axes.

  Notes
  -----
  
  """
  this_func = this_lib + 'Force_Aspect: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  # CHECK KEYWORD ARGUMENTS
  try:
    fig = kwargs['fig']
  except KeyError:
    fig = plt.gcf()

  try:
    ax = kwargs['ax']
  except KeyError:
    ax = plt.gca()  
  
  try:
    aspect = kwargs['aspect']
  except KeyError: 
    aspect = 1
  
  # FORMAT  
  im = ax.get_images()
  extent =  im[0].get_extent()
  ax.set_aspect(abs((extent[1] - extent[0]) / (extent[3] - extent[2])) / aspect)
  
  #print this_func + 'END'
  return fig, ax


# -------------------------------------------------------------------------------


def Pick_My_ColorMap(name):
  """
  Pick one of the custom, already 
  defined colormaps.
  
  Parameters
  ----------
  name : str 
    Name of the colormap. 
      
  Returns
  -------
  my_cmap : cm
    Colormap, formally a 'LinearSegmentedColormap'.

  Notes
  -----
  
  """
  this_func = this_lib + 'Pick_My_ColorMap: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  if name == 'R_Matrix': # Colormap used for plots of resolution matrices output from CPS package.
    
    colors = [(255, 255, 255), (0, 0, 255), (0, 246, 255), (255, 222, 0), (255, 0, 0)]
    #, (255, 0, 0), (255, 0, 0), (255, 0, 0)] # , (25, 0, 116)
    
    c1 = 'white'
    c2 = '#00035b' #'#004577'#'#00555a' #'#005f6a' #'petrol' #'#0e87cc' #water lbue#'#047495' #'sea blue''blue' #midnightblue'
    c3 = '#41fdfe' #'brightcyan' #'cyan'
    c4 = '#9cef43' #(kiwi) #lime
    c5 = '#f1da7a' #(sandy) #'yellow'
    c6 = '#ffb07c' #'peach' #'orange'
    c7 = '#fd8d49' #orangeish'orangered'
    c8 = 'red'
    c9 = 'maroon'
    c10 = 'black'
    
    c = mcolors.ColorConverter().to_rgb
    my_cmap = Make_Colormap(
        [c(c1), 0.01, c(c1), c(c2), 0.1, c(c2), c(c3), 0.2, c(c3), c(c4), 0.3, c(c4), 
         c(c5), 0.4, c(c5), c(c6), 0.5, c(c6), c(c7), 0.6, c(c7), c(c8), 0.7, c(c8), 
         c(c9), 0.8, c(c9), c(c10), 0.9, c(c10)])
  
  else:
    eprint(this_func + 'Error.  Unknown name of the cmap: ' + name + '\n')
    quit()
  
  #print this_func + 'END'
  return my_cmap


# -------------------------------------------------------------------------------


def Make_Colormap(seq):
  """
  Make a custom colormap provided 
  a sequence of RGB values and scaling.
  
  Parameters
  ----------
  seq : list
    A sequence of floats and RGB-tuples. 
    The floats should be increasing
    and in the interval (0,1).
      
  Returns
  -------
  my_cmap : cm
    Colormap, formally a 'LinearSegmentedColormap'.

  Notes
  -----
  
  """
  this_func = this_lib + 'Make_Colormap: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START') 
  
  seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
  
  cdict = {'red': [], 'green': [], 'blue': []}
  
  for i, item in enumerate(seq):
    if isinstance(item, float):
      r1, g1, b1 = seq[i - 1]
      r2, g2, b2 = seq[i + 1]
      cdict['red'].append([item, r1, r2])
      cdict['green'].append([item, g1, g2])
      cdict['blue'].append([item, b1, b2])
      
  my_cmap = mcolors.LinearSegmentedColormap('CustomMap', cdict)
  
  #print this_func + 'END'
  return my_cmap


# -------------------------------------------------------------------------------


def Shifted_Colormap(name, minn, maxx):
  """
  Assign a name to a shifted colormap.
  
  Parameters
  ----------
  name : str 
    It should consist of a name 
    of the original colormap (e.g. 'viridis') 
    followed by 0 (=> 'viridis0').
  minn : float 
    Min. amplitude of the data.
  maxx : float
    Max. amplitude of the data.
      
  Returns
  -------
  0

  Notes
  -----
  
  """
  this_func = this_lib + 'Shifted_Colormap: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  shift = True
  
  if name == 'seismic0':
    mapp = matplotlib.cm.seismic
  
  elif name == 'bwr0':
    mapp = matplotlib.cm.bwr
     
  elif name == 'BrBG0':
    mapp = matplotlib.cm.BrBG
  
  elif name == 'RdYlBu0':
    mapp = matplotlib.cm.RdYlBu
  
  else:
    eprint(this_func + 'Unknown colormap: ' + name + ' proceeding without a shift. \n')
    shift = False 
    #cmap = name
  
  if shift:
    mapp = Shift_Colormap(mapp, minn, maxx, start=0, stop=1.0, name=name) 
  
  #print this_func + 'END' 
  return 0


# -------------------------------------------------------------------------------


def Shift_Colormap(cmap, minn, maxx, start=0, stop=1.0, name='shiftedcmap'):
  """
  Function to offset the "center" of a colormap. 
  
  Parameters
  ----------
  cmap : cm 
    The matplotlib colormap to be altered
  minn : float 
    Min. amplitude of the data.
  maxx : float
    Max. amplitude of the data.
  start : float
    Offset from lowest point in the colormap's range.
    Should be between 0.0 and `midpoint` (see notes below).
    Default: 0.0 (no lower ofset). 
  stop : float
    Offset from highets point in the colormap's range.
    Should be between `midpoint`  (see notes below) and 1.0.
    Default: 1.0 (no upper ofset)
  name : str 
    Default: 'shiftedcmap'.
  
  Returns
  -------
  newcmap : cm
    Colormap after shifting. 
  
  Notes
  -----
  From stackoverflow.
  
  'Useful for data with a negative min and positive max 
  and you want the middle of the colormap's dynamic 
  range to be at zero.'
  
  midpoint : float
    The new center of the colormap. Defaults to 
    0.5 (no shift). Should be between 0.0 and 1.0. In
    general, this should be  1 - vmax/(vmax + abs(vmin))
    For example if your data range from -15.0 to +5.0 and
    you want the center of the colormap at 0.0, `midpoint`
    should be set to  1 - 5/(5 + 15)) or 0.75  
  
  """
  this_func = this_lib + 'Shift_Colormap: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')   
  
  from mpl_toolkits.axes_grid1 import AxesGrid
  
  print(this_func, 'cmap, minn, maxx, start, stop, name', cmap, minn, maxx, start, stop, name)
  
  
  if (maxx == 0) and (minn == 0):
    midpoint = 0.5
  else:
    midpoint = 1 - maxx / (maxx + abs(minn))
  
  print(this_func, 'midpoint', midpoint)
  
  cdict = {
      'red': [],
      'green': [],
      'blue': [],
      'alpha': []
  }

  # regular index to compute the colors
  reg_index = np.linspace(start, stop, 257)

  # shifted index to match the data
  shift_index = np.hstack([
      np.linspace(0.0, midpoint, 128, endpoint=False), 
      np.linspace(midpoint, 1.0, 129, endpoint=True)
  ])

  for ri, si in zip(reg_index, shift_index):
      r, g, b, a = cmap(ri)

      cdict['red'].append((si, r, r))
      cdict['green'].append((si, g, g))
      cdict['blue'].append((si, b, b))
      cdict['alpha'].append((si, a, a))

  newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
  plt.register_cmap(cmap=newcmap)

  #print this_func + 'END'
  return newcmap


# -------------------------------------------------------------------------------


def Colorbar(**kwargs):
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
  this_func = this_lib + 'Colorbar: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  ax = Kwarg('ax', plt.gca(), kwargs)
  
  from mpl_toolkits.axes_grid1 import make_axes_locatable
  
  divider = make_axes_locatable(ax)
  # create an axes on the right side of ax. The width of cax will be 5%
  # of ax and the padding between cax and ax will be fixed at 0.05 inch.
  cax = divider.append_axes("right", size="5%", pad=0.15)
  cbar = plt.colorbar(cax=cax, format='%.0e') # NOTE: IT WILL CRASH FOR FS! NEEDS cbar=0 IN THIS CASE
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
















# -------------------------------------------------------------------------------
# 3D-POINTS VISUALIZATION # FIXME: OBSOLETE?
# -------------------------------------------------------------------------------


def Plot_Points(points, **kwargs):
  """
  Framework function to plot 3D points 
  of the form [x, y, z].
  
  Parameters
  ----------
  points : list 
    List of 3D tuples [x, y, z].
  **kwargs : keyword arguments, optional except for:
    plot_type : str 
      Current capabilities:
       - 'scatter'
       - 'slice'
  
  Returns
  -------
  0
  
  Notes
  -----
  It checks if they are co-planar.
  
  """  
  this_func = this_func = this_lib + 'Plot_Points: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  ## CHECK IF THE POINTS ARE CO-PLANAR (2D)
  #ndims = 3
  #const_dim = -1
  #for i, coord in enumerate(range(ndims)):
  #  coord_value = points[0][coord]
  #  const = True
  #  for point in points:
  #    if point[coord] != coord_value:
  #      const = False 
  #      break
  #  if const:
  #    const_dim = i
  #    break
  
  # IF YES, PLOT THEM EASILY
  #if const_dim == 0:
  #  plt.scatter([i[1] for i in points], [i[2] for i in points])
  #  plt.title('x=const=' + str(coord_value))
  #  plt.xlabel('y')
  #  plt.ylabel('z')
  #
  #elif const_dim == 1:
  #  plt.scatter([i[0] for i in points], [i[2] for i in points])
  #  plt.title('y=const=' + str(coord_value))
  #  plt.xlabel('x')
  #  plt.ylabel('z')
  #  
  #elif const_dim == 2:
  #  plt.scatter([i[0] for i in points], [i[1] for i in points])
  #  plt.title('z=const=' + str(coord_value))
  #  plt.xlabel('x')
  #  plt.ylabel('y')
    
  # IF NOT, PLOT THEM IN 3D
  #if False:
  #  raise('Plotting of co-planar points not yet implemented')
  #else:
  
  Plot_Points_3D(points, **kwargs)
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Points_3D(points, **kwargs): 
  """
  Framework function to plot 3D points 
  of the form [x, y, z].
  
  Parameters
  ----------
  points : list 
    List of 3D tuples [x, y, z].
  **kwargs : keyword arguments, optional
    plot_type : str 
      Current capabilities:
       - 'scatter' (default)
       - 'slice'
       - 'projection'
  
  Returns
  -------
  0
  
  Notes
  -----
  To plot 2D data (nx2=1) just set yslice=0 (one of the kwargs).
  
  """  
  this_func = this_func = this_lib + 'Plot_Points_3D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  plot_type = Kwarg('plot_type', 'scatter', kwargs, 1)
    
  if plot_type == 'scatter': 
    Plot_Points_3D_Scatter(points, **kwargs)

  elif plot_type == 'slice': 
    Plot_Points_3D_Slice(points, **kwargs)
    
  elif plot_type == 'proj': 
    Plot_Points_3D_Projection_2D(points, **kwargs)
    
  else:
    raise ValueError('Unknown plot-type: ' + plot_type + '\n')
    
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Points_3D_Scatter(points, **kwargs): # BACKEND
  """
  
  Parameters
  ----------
  points : list 
    List of 3D tuples [x, y, z].  
  **kwargs : keyword arguments, optional
  
  Returns
  -------
  0
  
  Notes
  -----
  NOTE WE NEGLECT THE AMPLITUDE OF 4D POINTS HERE
  
  """  
  this_func = this_func = this_lib + 'Plot_Points_3D_Scatter: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  try: 
    ax = kwargs['ax']
  except KeyError:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')    
    
  c = Kwarg('c', 'k', kwargs)
  s = Kwarg('s', 10, kwargs)  
  alpha = Kwarg('alpha', 1, kwargs)
  

  X = [i[0] for i in points]
  Y = [i[1] for i in points]
  Z = [i[2] for i in points]
  # NOTE WE NEGLECT THE AMPLITUDE OF 4D POINTS HERE

  ax.scatter(X, Y, Z, c=c, s=s, alpha=alpha)
  Format_Plot_3D(**kwargs)
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Points_3D_Slice(points, **kwargs): # BACKEND
  """
  Framework function to plot
  
  Parameters
  ----------
  points : list 
    List of 3D tuples [x, y, z].  
  **kwargs : keyword arguments, optional
      Just passing it down.
  
  Returns
  -------
  0
  
  Notes
  -----
  To plot 2D data (nx2=1) just set yslice=0 (one of the kwargs).
  
  """  
  this_func = this_func = this_lib + 'Plot_Points_3D_Slice: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  
  ax = Kwarg('ax', plt.gca(), kwargs)  
  
  #npoints = []
  #
  #try:
  #  xslice = kwargs['xslice']
  #  for p in points:
  #    x = p[0]
  #    y = p[1]
  #    z = p[2]
  #    if x == xslice:
  #      npoints.append([y, z])
  #      
  #  plt.title('Slice for in-line node no. ' + str(xslice))
  #  plt.xlabel('cross-line node')
  #  plt.ylabel('depth node')
  #except KeyError:
  #  try:
  #    yslice = kwargs['yslice']
  #    for p in points:
  #      x = p[0]
  #      y = p[1]
  #      z = p[2]
  #      if y == yslice:
  #        npoints.append([x, z])      
  #    
  #    plt.title('Slice for cross-line node no. ' + str(yslice))
  #    plt.xlabel('in-line node')
  #    plt.ylabel('depth node')
  #  except KeyError:
  #    try:
  #      zslice = kwargs['zslice']
  #      for p in points:
  #        x = p[0]
  #        y = p[1]
  #        z = p[2]
  #        if z == zslice:
  #          npoints.append([x, y])
  #          
  #      plt.title('Slice for depth node no. ' + str(zslice))
  #      plt.xlabel('in-line node')
  #      plt.ylabel('cross-line node')
  #    except KeyError:
  #      eprint(this_func + 'Error. No slice-coordinate (x / y / z) defined for plot_type=slice\n')
  #      quit()
  #
  #plotx = [i[0] for i in npoints]
  #ploty = [i[1] for i in npoints]
  #
  #
  ## FORMAT THE PLOT
  #try: 
  #  c = kwargs['c']
  #except KeyError:
  #  c = 'k'
  #
  #try: 
  #  s = kwargs['s']
  #except KeyError:
  #  s = 10
  #
  #try: 
  #  alpha = kwargs['alpha']
  #except KeyError:
  #  alpha = 1.   
  
  from lib_fwi_generic import Slice_Points
  
  points2d = Slice_Points(points)
  pointsx = [i[0] for i in points2d]
  pointsy = [i[1] for i in points2d]
  ax.scatter(pointsx, pointsy)#, c=c, s=s, alpha=alpha)
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Points_3D_Projection_2D(points, **kwargs): # BACKEND
  """
  Plot 3D points projected on a chosen plane.
  
  Parameters
  ----------
  points : list 
    List of 3D tuples [x, y, z].  
  **kwargs : keyword arguments, optional
  
  Returns
  -------
  0
  
  Notes
  -----
  
  """  
  this_func = this_func = this_lib + 'Plot_Points_3D_Projection_2D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  try: 
    plane = kwargs['plane']
  except KeyError:
    eprint(this_func + "Warning. plane not specified - assuming 'XY'. \n")
    plane = 'XY'
  
  # FORMAT THE PLOT
  try: 
    c = kwargs['c']
  except KeyError:
    c = 'k'
  
  try: 
    s = kwargs['s']
  except KeyError:
    s = 10

  try: 
    alpha = kwargs['alpha']
  except KeyError:
    alpha = 1. 
    
  label = Kwarg('label', None, kwargs)
    
  x, y, z = Split_Tuples(points)
  #print this_func, x,y,z
  
  if plane == 'XY':
    X = x
    Y = y
    plt.xlabel('x [nodes]')
    plt.ylabel('y [nodes]')  

  elif plane == 'XZ':
    X = x
    Y = z
    plt.xlabel('x [nodes]')
    plt.ylabel('z [nodes]') 
  
  elif plane == 'YZ':
    X = y
    Y = z
    plt.xlabel('y [nodes]')
    plt.ylabel('z [nodes]') 

  # PLOT
  if c != 'rainbow':
    plt.scatter(X, Y, c=c, s=s, alpha=alpha, label=label)
  
  else:
    colors = iter(cm.rainbow(np.linspace(0, 1, len(X))))
    for i in range(len(X)):
      plt.scatter(X[i], Y[i], c=next(colors), s=s, alpha=alpha)

  #print this_func + 'END'
  return 0



# -------------------------------------------------------------------------------

