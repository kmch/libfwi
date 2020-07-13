"""
Author: Kajetan Chrapkiewicz, 2018. 
All rights reserved. Ask for permision writing to K.Chrapkiewicz17@imperial.ac.uk.

This library provides framework procedures for...
At the moment fully implemented is:

1. ...

"""

## MODULES
import numpy as np
#import matplotlib.pyplot as plt

## MY MODULES
from lib_generic import *
#from lib_fwi_generic import *
from lib_generic_CONST import verbos_func #*

from fullwavepy.generic.decor import timer
##


## CONSTS
this_lib = 'lib_io_fullwave.py/'

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
# CONVERSIONS BETWEEN FORMATS (APIs)
# -------------------------------------------------------------------------------


def Convert_xyz2vtr(file_name, dr, **kwargs):
  """
  Read an .xyz file and save it as .vtr. 
  
  Parameters
  ----------
  file_name : string
    File name. It should include .xyz extension. It can include path if needed.
  dr : int 
    Size of the grid cell. All X and Y coordinates present in the file must be 
    its multiples and their units must match.
  
  **kwargs : keyword arguments, optional
    z_origin : int
      z-coordinate of the origin of a new, shifted frame of coordinates.
      For all z: z_new = z_old + z_origin.
  
  Returns
  -------
  0
  
  Notes
  -----
  Assumptions about .xyz file:
  - sea level is at z=0, land is at z>0
  NOTE: READY TO RELEASE!
  """
  this_func = this_lib + 'Convert_xyz2vtr: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_generic_CONST import very_big
  
  # CHECK KEYWORD ARGUMENTS
  z_origin = Kwarg('z_origin', 0, kwargs)    
  z_clip_min = Kwarg('z_clip_min', -very_big, kwargs)
  z_clip_max = Kwarg('z_clip_max', very_big, kwargs)
  z_scale = Kwarg('z_scale', -1, kwargs)
  
  print(this_func, 'Assuming xyz contains elevation not depth.')
  
  # READ THE FILE    
  x, y, z = Read_xyz(file_name, dr)
  
  # CLIP POINTS BELOW SEA LEVEL, I.E CONVERT DATA TO THE FREE-SURFACE TYPE
  z = np.clip(z, z_clip_min, z_clip_max) 
  
  # SCALE (E.G. FLIP FOR z_scale=-1)
  z *= z_scale
  
  # TRANSFORM COORDINATES TO NODES 
  z /= float(dr)
  
  # SHIFT THE ORIGIN OF THE COORDINATES FRAME
  z += z_origin
  
  # CONFORM TO .vtr CONVENTION (MORE SENSIBLE AND COMMON ANYWAY)
  z = z.T 
  z = [[[j] for j in i] for i in z] # WE NEED A TRACE (NOT A SINGLE NUMBER) IN THE FASTEST DIMENSION
  z = np.array(z) # BACK TO ARRAY
  
  # GET DIMENSIONS
  nx = z.shape[0]
  ny = z.shape[1]
  nz = 1
  
  # STRIP THE EXTENSION ('core' WOULD NOT WORK BECAUSE OF PATHS CONTAINING WHITE SPACE
  file_name = file_name[ :-len('.xyz')] 
  
  # SAVE AS .vtr FILE
  Save_vtr(file_name, nx, ny, nz, z) # NOT -FreeSurf BECAUSE WE MIGHT NOT CLIP TO SEA SURFACE

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Convert_sgy2vtr(fname, **kwargs): # FIXME
  """
  # NOTE:
  # MODELS NEED TO BE HANDLED BY convert_sgy2vtr_3D.sh  
  
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
  this_func = this_lib + 'Convert_sgy2vtr: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  from lib_su import SU_Get_Traces_No
  
  ntraces = SU_Get_Traces_No(fname)
  print(this_func, 'ntraces', ntraces)
  #o = Bash('sgy2vtr.sh ' + fname + ' ' + str(ntraces))  
  #print this_func, o
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Convert_Plain_Txt_2_Model_3D(nx, ny, plain_list): # NOTE: D-C
  """
  Model 3D is in .vtr-like format.
  Plain txt is with z as fast index (inner loop),
  x as a slowest one (outer loop) 
  
  # FREQUENTLY USED!
  """
  
  this_func = this_lib + 'Convert_Plain_Txt_2_Model_3D: '
  print(this_func + 'START')

  model = []
  i = 0
  for x in np.arange(nx):
    xtraces = []
    for y in np.arange(ny):
      xtraces.append(plain_list[i])
    model.append(xtraces)  

  print(this_func + 'END')
  return model


# -------------------------------------------------------------------------------


def Convert_Synthetic_sgy(proj_name, **kwargs): # FIXME: D-C
  """
  Split and convert to vtr.
  
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
  this_func = this_lib + 'Convert_Synthetic_sgy: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  
  path = Kwarg('path', '../output/', kwargs)
  
  key_station = Kwarg('key_station', key_station_proteus, kwargs)
  key_line = Kwarg('key_line', key_line_proteus, kwargs)
  fname = Kwarg('fname', path + proj_name + '-Synthetic.sgy', kwargs)

  o = Bash('su_supergather_split_stations_n_lines.sh ' + 
           key_station + ' ' + key_line + ' ' + fname)
  if verbos > 1:
    print(this_func, o)
  #Supergather_sgy_Split(fname, **kwargs)
  #fnames = Get_Files(...)
  #Convert_sgy2vtr(fnames)

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
# I/O SEGY
# -------------------------------------------------------------------------------


#def Save_sgy(**kwargs): # EMPTY


# -------------------------------------------------------------------------------


def Read_sgy(fname, **kwargs):
  """
  Converts to vtr.
  
  """
  this_func = this_lib + 'Read_sgy: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
    
  nx = Kwarg('nx', None, kwargs)
  dim = Kwarg('dim', 3, kwargs)
  fname_vtr = fname[:-len('.sgy')] + '.vtr'
    
  if dim == 3:
    if nx:
      if verbos > 0:
        print(this_func, 'Converting sgy2vtr for file: ' + fname)
      o, e = Bash2('convert_sgy2vtr_3D.sh ' + str(int(nx)) + ' ' + fname)
    
    if Exists(fname_vtr, **kwargs):
      nx, ny, nz, Z = Read_vtr(fname_vtr, **kwargs)
      eprint(this_func + 'Plotting ' +  fname_vtr  + 
             '. Provide nx to convert sgy (nx of ' + fname_vtr + 
             ' is ' + str(nx) + ').\n')
      
    else:
      eprint(this_func + 'File in .vtr format: ' + fname_vtr + ' not found. Skipping this figure. You can provide nx to automatically convert sgy 2 vtr\n')  
      Z = []
  
  #elif dim == 2:
    #if verbos > 0:
      #print this_func, 'Converting sgy2vtr for file: ' + fname
    #o, e = Bash2('convert_sgy2vtr_2D.sh ' + fname)
    
  
  
  else:
    raise ValueError('Wrong no. of dimensions: ' + str(dim))
  
  #print this_func + 'END'  
  return Z
  

# -------------------------------------------------------------------------------


#def Save_geo(...)


# -------------------------------------------------------------------------------


def Read_geo(file_name, dx, **kwargs):
  """
  Read Fullwave's .geo geometry files.
  
  Parameters
  ----------
  file_name : str    
    File name. It should include  extension. 
    It can include path if needed.
  dx : float 
    Size of the grid cell in metres.

  Returns
  -------
  dims : list
    Dimension of the grid as in the geom.
    files.
  records : dict 
    records[id] = [x, y, z]
  
  Notes
  -----
  Merge with Read_pgy (DRY)?
  
  NOTE: THESE ARE NOT PHYSICAL DIMENSION, ONLY GRID DIMS (BUT IN metres)
  
  """
  this_func = this_lib + 'Read_geo: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  # CONVERT TO NODES (YES BY DEFAULT FOR COMPATIBILITY)
  convert = Kwarg('convert', True, kwargs) 
  
  if verbos > 1:
    print(this_func, 'File to read: ', file_name)
  
  content = Read_File(file_name)
  header = content[0]
  data = content[1: ]
  
  dims = [int(float(i) / dx) for i in header[1: ]]
  
  records = {}
  for row in data:
    # NOTE: x AND z ARE NOT SWAPPED IN THIS FORMAT (IN CONTRAST TO .pgy)
    # NOTE : '1 +' ENSURES CONSISTENCY WITH .pgy FILES  
    #print 'row', row
    
    if convert:
      records[row[0]] = [1 + float(row[1]) / dx, 1 + float(row[2]) / dx, 1 + float(row[3]) / dx]  
    else: #NOTE: THESE ARE NOT PHYSICAL DIMENSION, ONLY GRID DIMS (BUT IN metres)
      records[row[0]] = [float(row[1]), float(row[2]), float(row[3])]    
    #print 'records[row[0]]', records[row[0]]
  
  #print this_func + 'END'
  return dims, records#points
  

# -------------------------------------------------------------------------------
# I/O FW3D
# -------------------------------------------------------------------------------


def Save_pgy(fname, xyz_list, dims, **kwargs):
  """
  Save Fullwave's .pgy geometry files.
  
  Parameters
  ----------
  file_name : str    
    File name. It should include  extension. 
    It can include path if needed.

  Returns
  -------
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Save_pgy: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  nx1, nx2, nx3 = dims
  
  n = len(xyz_list) # NO. OF POINTS
  sn   = '{:10}'.format(n)
  snx3 = '{:15}'.format(nx3)
  snx2 = '{:15}'.format(nx2)
  snx1 = '{:15}'.format(nx1)
  
  f = open(fname, 'w')
  f.write(sn + snx3 + snx2 + snx1 + '\n')
  for i, xyz in enumerate(xyz_list):
    si   = '{:10}'.format(i + 1)
    sz   = '{:15.8f}'.format(xyz[2])
    sy   = '{:15.8f}'.format(xyz[1])
    sx   = '{:15.8f}'.format(xyz[0])  
    #print si, sz, sy, sx
    f.write(si + sz + sy + sx + '\n')
    i += 1
  f.close()  
  
  #print this_func + 'END'
  return 0
  

# -------------------------------------------------------------------------------


def Read_pgy(file_name, **kwargs):
  """
  Read Fullwave's .pgy geometry files.
  
  Parameters
  ----------
  file_name : str    
    File name. It should include  extension. 
    It can include path if needed.

  Returns
  -------
  nx1, nx2, nx3 : int
    Dimension of the grid as in the geom.
    files.
  points : list 
    List of points [x, y, z]
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Read_pgy: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  
  
  if verbos > 2:
    print(this_func, 'File to read: ', file_name)
  
  content = Read_File(file_name)
  header = content[0]
  data = content[1: ]
  nx3, nx2, nx1 = [int(float(i)) for i in header[1: ]]
  
  records = {}
  for row in data:
    # NOTE: x AND z ARE SWAPPED IN THIS FORMAT
    records[row[0]] = [float(row[3]), float(row[2]), float(row[1])]

  #print this_func + 'END'
  return nx1, nx2, nx3, records
  

# -------------------------------------------------------------------------------


def Save_vtr(file_name_core, nx1, nx2, nx3, model, **kwargs):
  """
  Save gridded data structure ('model') 
  as a .vtr file.
  
  Parameters
  ----------
  file_name_core : string
    File name without extension.
    It can include path if needed.
  nx1, nx2, nx3 : int / float
    Dimensions of the 'model' 3D array (see below) 
    from the slowest (nx1) to the 
    fastest (nx3) index.
  model : array / list
    Gridded 3D data structure, 
  
  Returns
  -------
  0
  
  Notes
  -----
  The order of indices is the same as in 'fullwave3d' code by M. Warner et al.
  
  It might not be compatibile with 'wavefor2d' code by L. Guasch.
  
  """  
  this_func = this_lib + 'Save_vtr: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from scipy.io import FortranFile
  
  file_name = file_name_core + '.vtr'
  
  if verbos > 2:
    print(this_func + 'File to save: ', file_name)
    
  model = np.array(model)
  
  nx1 = int(nx1)
  nx2 = int(nx2)
  nx3 = int(nx3)
  
  f = FortranFile(file_name, 'w')
  
  if nx2 == 1:
    ndims = 2
    f.write_record(np.array([1, ndims, 0], dtype = np.int32))
    f.write_record(np.array([nx3, nx1], dtype = np.int32))
  elif nx2 > 1:
    ndims = 3
    f.write_record(np.array([1, ndims, 0], dtype = np.int32))  
    f.write_record(np.array([nx3, nx2, nx1], dtype = np.int32)) #NOT COMPATIBLE WITH LLUIS CODE
  else:
    eprint(this_func + 'Error. nx2 cannot be < 1: ' + str(nx2) + '\n')
    quit()
    
  for x in range(int(nx1)):
    for y in range(int(nx2)):
      trace = model[x, y, :]
      f.write_record(np.array(trace, dtype = np.float32))
  
  f.close()  
  
  # DEBUGG
  if verbos > 3:
    print(this_func + 'ndims: ', ndims)
    print(this_func + 'nx3, nx2, nx1: ', nx3, nx2, nx1)

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------

@timer
def Read_vtr(file_name, **kwargs):
  """
  Read a .vtr file.
  
  Parameters
  ----------
  file_name : str    
    File name. It should include  extension. 
    It can include path if needed.
  **kwargs : keyword arguments (optional)
    Current capabilities:
    mode : int 
      
  
  Returns
  -------
  nx1, nx2, nx3 : int
    Dimensions of the model (see below) read from the header.
  model : array
    Gridded 3D data structure of shape that should be consistent 
    with nx1, nx2, nx3 (see above). It should be such as long as 
    the file was written with the Save_vtr function.
  
  Notes
  -----
  .vtr is a Fullwave3D's native binary format. 
  
  """
  this_func = this_lib + 'Read_vtr: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  from scipy.io import FortranFile
  
  
  
  if verbos > 1:
    print(this_func, 'File: ', file_name)
  
  # CHECK KEYWORD ARGUMENTS
  try:
    mode = kwargs['mode'] 
  except KeyError:
    mode = 1
    
  if mode == 1:
    f = FortranFile(file_name, 'r')
      
    ncomponents, ndims, dummy_zero = f.read_ints(dtype = np.int32)
    #print this_func, 'ncomponents, ndims, dummy_zero', ncomponents, ndims, dummy_zero
    
    if ndims == 3: 
      nx3, nx2, nx1 = f.read_ints(dtype = np.int32)
    elif ndims == 2:
      nx3, nx1 = f.read_ints(dtype = np.int32)
      nx2 = 1
    elif ndims == 1:
      nx3, nx1 = f.read_ints(dtype = np.int32)
      nx2 = 1
      
    else:
      eprint(this_func + 'Error. Wrong ndims: ' + str(ndims) + '\n')
      quit()
    
    model = []
    for x in range(int(nx1)):
      xslices = []
      for y in range(int(nx2)):
        trace = f.read_reals(dtype = np.float32)
        xslices.append(trace)
      model.append(xslices)    
    
  elif mode == 2: # FIXME: DELETE?
    try:
      f = FortranFile(file_name, 'r')
      header1 = f.read_ints(dtype = np.int32)
      nx3, nx1 = f.read_ints(dtype = np.int32)
      nx2 = 1
      model = [] #np.zeros((nx1, nx2))
      for x in np.arange(nx1):
         model.append(f.read_reals(dtype = np.float32))
      model = np.transpose(model) # PLOTTING REQUIRES THIS
      
    except ValueError:
      f.close() # OTHERWISE IT WOULD READ THE NEXT LINE (3RD)
      f = FortranFile(file_name, 'r')
      header1 = f.read_ints(dtype = np.int32)
      nx3, nx2, nx1 = f.read_ints(dtype = np.int32)
      
      model = [] #np.zeros((nx1, nx2))
      for x in np.arange(nx1):
        xtraces = []
        for y in np.arange(nx2):
          xtraces.append(f.read_reals(dtype = np.float32))
        model.append(xtraces)
   
  elif mode == 3: # FIXME: DELETE?
    # FROM fsprep
    f = FortranFile(file_name, 'r')
    ncomponents, ndims, dummy_zero = f.read_ints(dtype = np.int32)    
    if ndims == 3: 
      nx3, nx2, nx1 = f.read_ints(dtype = np.int32)
    
    elif ndims == 2:
      nx3, nx1 = f.read_ints(dtype = np.int32)
      nx2 = 1
    
    else:
      eprint(this_func + 'Error. Wrong ndims: ' + str(ndims) + '\n')
      quit()
    
    model = []
    for x in range(int(nx1)):
      xslices = []
      for y in range(int(nx2)):
        trace = f.read_reals(dtype = np.float32)
        print(this_func, 'trace: ', trace)
        quit()
        xslices.append(trace)
      model.append(xslices)       
     
  f.close()
  
  # CONVERT LIST->ARRAY
  model = np.array(model)
  
  # DEBUGG
  nx1_mod = model.shape[0]
  nx2_mod = model.shape[1]
  nx3_mod = model.shape[2]
  
  if (nx1_mod != nx1 or nx2_mod != nx2 or nx3_mod != nx3) and (verbos > 0):
    eprint(this_func + 'Shape of the model inconsistent with the header\n')
    eprint(this_func + 'nx1 (model): ' + str(nx1_mod) + ', nx1 (header): ' + str(nx1) + '\n')
    eprint(this_func + 'nx2 (model): ' + str(nx2_mod) + ', nx2 (header): ' + str(nx2) + '\n')
    eprint(this_func + 'nx3 (model): ' + str(nx3_mod) + ', nx3 (header): ' + str(nx3) + '\n')
  
  
  # HANDY REPORT
  if verbos > 2:
    print(this_func + 'min value: ', np.min(model))
    print(this_func + 'max value: ', np.max(model))
    print(this_func + 'nx1, nx2, nx3: ', nx1, nx2, nx3)
  
  #print this_func + 'END'
  return nx1, nx2, nx3, model


# -------------------------------------------------------------------------------


def Save_ttr(proj_name, nshots, nrecs, dt, nsamp): # FIXME
  from scipy.io import FortranFile
  # CURRENTLY AIMS AT Observed-0000.ttr FILE ONLY
  # FIXME NOT FINISHED
  file_name = proj_name + '-Observed-0000.ttr'
  f = FortranFile(file_name, 'w')
  f.write_record(np.array([nshots, nrecs, nsamp, dt * nsamp], dtype = np.int32)) 
  for index1 in (list(range(nshots)) + 1):
    print(index1)
  
  
  f.close()


# -------------------------------------------------------------------------------


def Read_ttr(file_name): # FIXME
  print('lib_fwi_generic.py/read_ttr: START')
  #import struct
  #
  #
  #with open(file_name, "rb") as f:
  #  
  #  byte = f.read(4) # READ
  #  i = 1
  #  while i < 5:#byte != "": # SINCE read OUTPUTS A STRING!
  #    print struct.unpack('i', byte)[0]
  #    byte = f.read(4)
  #    i += 1
  #    #byte = f.read(1)
  #  #while byte != "": # SINCE read OUTPUTS A STRING!
  #  while i < 20:
  #    print struct.unpack('f', byte)[0]
  #    byte = f.read(4)
  #    i += 1
  #  #  byte = f.read(4)
  #  #  print struct.unpack('f', byte)[0]
  #    
  #quit()
  #byte = f.read(1)
  #with open(file_name, "rb") as f:
  #  print(struct.unpack('i', f.read(4)))
  #  quit()
  #  while byte != "":
  #      # Do stuff with byte.
  #      byte = f.read(1)
  #      print byte
  #      quit()
  #quit()
  
  
  #from scipy.io import FortranFile
  #f = FortranFile(file_name, 'r')
  #
  #model = []
  #i = 0
  #i_start = i
  #
  #
  #while True:
  #  trace = [] # EMPTY JUST IN CASE
  #  try:
  #    trace = f.read_reals(dtype = np.float32)
  #  except TypeError:
  #    eprint('EOF reached.\n')
  #    break
  #  if i > i_start: # SKIP FIRST LINE (INDEX OF WHICH IS i_start = 0) WHICH IS ONLY 4-FLOATS LONG
  #    model.append(trace)
  #  i += 1
  ##model = np.transpose(model)
  #print len(model)
  #
  #f.close()
  with open(file_name,'r') as f:
    header = np.fromfile(f, dtype=np.int32, count=4)
    time = np.fromfile(f, dtype=np.float32, count=1)
    data = np.fromfile(f, dtype=np.float32)
  
  
  print(header)
  print(time)
  print(data[:7])
  #quit()
  print('lib_fwi_generic.py/read_ttr: END')
  return data[4:]#model
  

# -------------------------------------------------------------------------------
# OTHERS
# -------------------------------------------------------------------------------


#def Save_xyz(...)


# -------------------------------------------------------------------------------


def Read_xyz(file_name, dr, **kwargs):
  """
  Read file in a .xyz format (used by GMT) and output gridded X, Y, Z ready to plot etc.
  
  Parameters
  ----------
  file_name : file name 
    It should include .xyz extension. It can include path if needed.
  dr : size of the grid cell
    All X and Y coordinates present in the file must be its multiples and their units must match.
  
  Returns
  -------
  X, Y, Z : arrays
    Description below.
  
  X : array 
    (nx, ny) shaped array being an x-component of the meshgrid(x, y),
    where x, y are lists of coordinates read from the file_name,
    and nx, ny are their respective lengths.
  Y : array
    (nx, ny) shaped array being an y-component of the meshgrid(x, y),
    where x, y are lists of coordinates read from the file_name,
    and nx, ny are their respective lengths.
  Z : array
    (nx, ny) shaped array being reshaped list of x-coordinates.
  
  Notes
  -----
  FIXME: only 3D now
  
  It assumes that the file was written in the following way:
  
  >>> for y in np.linspace(ymax, ymin, ny): # (np.a)range DOES NOT WORK WHEN y_start > y_end
  >>>   for x in range(xmin, xmax, nx):
  >>>      f.write(x, y, z[x][y])
  
  which is an output of:
  
  >>> gmt grd2xyz santorini_merged_local_50.grd -R' '$xmin/$xmax/$ymin/$ymax > output.xyz
  
  where  'santorini_merged_local_50.grd' is the bathymetry file prepared by MP.
  
  NOTE: y is sorted in a decreasing order which is INCOMPATIBLE with plotters etc.; the function reverts this order
  
  NOTE: treating x as the fastest index is a default option in np.meshgrid too:
  
  >>> x = range(xmin, xmax)
  >>> y = range(ymin, ymax)
  >>> X, Y = np.meshgrid(x, y)
  >>> Z = fs.reshape(X.shape)

  """
  this_func = this_lib + 'Read_xyz: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  from lib_generic import Read_File, Get_Chunks, Flatten_List
  from lib_math_num import Metres2nodes
  
  c = Read_File(file_name)

  # SPLIT [[x1, y1, z1], ..., [xN, yN, zN]] INTO: [x1, ..., xN], [y1, ..., yN], [z1, ..., zN]
  x, y, z = list(zip(*c)) 
  
  
  # GET nx, ny
  x = list(map(float, x)) # = [float(i) for i in x]
  y = list(map(float, y)) #
  #z = map(float, z) #  
  xmin, xmax = min(x), max(x)
  ymin, ymax = min(y), max(y)
  #zmin, zmax = min(z), max(z)
  print('xmin, xmax', xmin, xmax)
  print('ymin, ymax', ymin, ymax)
  #print 'zmin, zmax', zmin, zmax
  xdiff = xmax - xmin # THIS WORKS IN ALL CASES PROVIDED xmax >= xmin
  ydiff = ymax - ymin ##
  print('xdiff', xdiff)
  print('ydiff', ydiff)
  nx = Metres2nodes(xdiff + dr, dr) # '+ dr' BECAUSE E.G. FOR xmin=50, xmax=150 WE WANT 3 NODES, NOT 2 
  ny = Metres2nodes(ydiff + dr, dr)
  print('nx, ny', nx, ny)
  
  
  # CORRECT THE SORTING ORDER OF Y-COORDINATE
  yslices = Get_Chunks(c, nx)
  yslices_list = []
  for yslice in yslices:
    yslices_list.append(yslice)
  
  yslices_list = yslices_list[ : : -1] # REVERT THE LIST
  new_c = Flatten_List(yslices_list)
    
  
  # EXTRACT Z-COORD
  fs = [float(i[2]) for i in new_c]
  
  
  # PREPARE THE MESHGRID
  x = np.linspace(xmin, xmax, nx)
  y = np.linspace(ymin, ymax, ny)
  X, Y = np.meshgrid(x, y)
  
  
  # TRANSFORM Z-COORD LIST INTO A GRIDDED ARRAY CORRESPONDING TO X AND Y 
  fs = np.array(fs)
  Z = fs.reshape(X.shape)
  
  
  # DEBUGG
  a = 5
  print(X[0][:a])
  print(Y[0][:a])
  print(Z[0][:a])

  #print this_func + 'END'
  return X, Y, Z


# -------------------------------------------------------------------------------

