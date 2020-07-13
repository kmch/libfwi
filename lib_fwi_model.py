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
#from lib_fwi_generic import Save_vtr, Read_vtr
from lib_io_fullwave import *
from lib_fwi_project_CONST import dims_default
##


## CONSTS
this_lib = 'lib_fwi_model.py/'

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
# CREATE A NEW MODEL
# -------------------------------------------------------------------------------


def Model_Create(fname, **kwargs):
  """
  Create a model, i.e. create new or copy
  an existing one.
  
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
  this_func = this_lib + 'Model_Create: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  #pars_gen = Kwarg('pars_gen', {}, kwargs)
  #pars_mod = Kwarg('pars_mod', {}, kwargs)
  
  
  nx1, nx2, nx3 = Kwarg('dims', dims_default, kwargs)# pars_gen)
  vel = Kwarg('vel', 5000, kwargs)#pars_mod)
  inhomo = Kwarg('inhomo', 'homo', kwargs)#pars_mod)
  vel_lambda = Kwarg('vel_lambda', 10, kwargs)#pars_mod)
  vel_ampl = Kwarg('vel_ampl', 0.1, kwargs) #pars_mod)
  
  model = np.ones(shape=(nx1, nx2, nx3), dtype=np.float32) * vel
  
  for x in range(nx1):
    for y in range(nx2):
      for z in range(nx3):
  
        if inhomo == 'homo':
          continue # ALREADY CORRECT VALUE
        
        elif inhomo == 'vertical':
          print() 
          
        elif (inhomo == 'lateral') or (inhomo == 'full'):
          if inhomo == 'lateral':
            arg = np.sqrt(x**2 + y**2)
          
          else:
            arg = np.sqrt(x**2 + y**2 + z**2)
            
          model[x][y][z] = model[x][y][z] + vel_ampl * model[x][y][z] * np.sin(arg * 2 * np.pi / vel_lambda)
          
        else:
          eprint(this_func + 'Error. Unknown inhomogeneity type: ' + inhomo + '\n')
          quit()
  
  Save_vtr(fname[ :-len('.vtr')], nx1, nx2, nx3, model, **kwargs)  # NOTE: 
  
  #print this_func + 'END'
  return model


# -------------------------------------------------------------------------------
# CHECK THE EXISTING MODEL (EXTRACT SOME INFO)
# -------------------------------------------------------------------------------


def Model_Check_Sea_Level(model, v_air, v_h2o, **kwargs):
  """
  In a 3D velocity model with an air-layer on top of it
  find the coordinates of the sea-level.
  
  Parameters
  ----------
  model : array 
    Gridded 3D array with Z coordinate being the fastest index 
    (v(x,y,z) = model[x-1][y-1][z-1])
  v_air : float 
    The velocity of the air-layer. It is assumed to be homogeneous and
    different from any non-air value.
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  Returns
  -------
  z_sea_xy : array
    Array of indices of the first water-node
    (within uncertainty) for each XY coord.
  
  Notes
  -----
  Python counts from 0!!!
  
  It gives the index of the first water-node.
  
  Z axis is assumed to point downwards => air layer starts 
  at z=1 nodes.
  """
  this_func = this_lib + 'Model_Check_Sea_Level: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  v_h2o_uncertainty = Kwarg('v_h2o_uncertainty', 0, kwargs)
  v1 = v_h2o - v_h2o_uncertainty
  v2 = v_h2o + v_h2o_uncertainty
  
  print(this_func, 'Searching for water velocity between ' + str(v1) + ' and ' + str(v2) + ' m/s')  
  
  if v_air == v_h2o:
    raise ValueError('v_air equal to v_h2o')
  
  nx, ny, nz = model.shape
  z_sea_xy = np.zeros((nx, ny)) - 1
  xy_not_sea = []
  for ix in range(nx):
    for iy in range(ny):
      found = False
      for iz in range(nz):
        v_mod = model[ix][iy][iz]
        if (v_mod > v1) and (v_mod < v2):
          z_sea_xy[ix][iy] = iz
          found = True
          break
      if not found:
        xy_not_sea.append([ix, iy])
  
  eprint(this_func + 'Length of list of XY indices for which water layer was not found within uncertainty: ' + str(len(xy_not_sea)) + '\n')
  #print this_func, 'min, max model', np.min(model), np.max(model)  
  z_sea_xy = Array2D_Convert_To_3D(z_sea_xy)
  
  #print this_func + 'END'
  return z_sea_xy


# -------------------------------------------------------------------------------


def Model_Check_Highest_Topo(model, v_air, **kwargs):
  """
  In a 3D velocity model with an air-layer on top of it
  find the coordinates of the highest-altitude,
  non-air point.
  
  Parameters
  ----------
  model : array 
    Gridded 3D array with Z coordinate being the fastest index 
    (v(x,y,z) = model[x-1][y-1][z-1])
  v_air : float 
    The velocity of the air-layer. It is assumed to be homogeneous and
    different from any non-air value.
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  Returns
  -------
  topo_min : vector
    Vector of coordinates of the highest peak:
    [x, y, z]
  
  Notes
  -----
  Z axis is assumed to point downwards => air layer starts 
  at z=1 nodes.
  """
  this_func = this_lib + 'Model_Check_Highest_Topo: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  topo_min = np.zeros(3)
  z_topo_min = big # HIGHEST TOPOGR. PEAK MINIMIZES (!) Z BECAUSE VERTICAL AXIS IS FLIPPED
  for ix in range(len(model)):
    for iy in range(len(model[0])):
      for iz in range(len(model[0][0])):
        if model[ix][iy][iz] > 0: # FIXME: IT ASSUMES vel=0 ABOVE FS
          break
      z_topo = iz
      if z_topo < z_topo_min:
        z_topo_min = z_topo
        topo_min = np.array([ix, iy, iz])
  
  print('topo_min (=highest peak) coordinates (x, y, z):', topo_min)
  
  #print this_func + 'END'
  return topo_min


# -------------------------------------------------------------------------------
# MODIFY THE EXISTING MODEL
# -------------------------------------------------------------------------------


def Model_Modify(prefix, model, **kwargs):
  """
  
  Parameters
  ----------
  
  Returns
  -------
  
  Notes
  -----
  """
  this_func = this_lib + 'Model_Modify_Extend: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  pars_mod = Kwarg('pars_mod', {}, kwargs)
  modify = Kwarg('modify', None, pars_mod)
  
  #if modify == 'extend':
    #Model_Modify_Extend(prefix, model)
  
  if not modify:
    pass
  
  else:
    raise NotImplementedError('modify: ' + modify)
  
  return model


# -------------------------------------------------------------------------------


def Model_Modify_Set_Air_Vel(model, old_v_air, mode, **kwargs):
  """
  Change air-layer velocity of the 3D model. 
  
  Parameters
  ----------
  model : array 
    Gridded 3D array with Z coordinate being the fastest index 
    (v(x,y,z) = model[x-1][y-1][z-1])
  old_v_air : float 
    The old velocity of the air-layer. It is assumed to be homogeneous and
    different from any non-air value.
  mode : str 
    May be 'homo' / 'extrapol'
  **kwargs : keyword arguments (optional)
    Current capabilities:
    - new_v_air : float 
      The new velocity of the air-layer. It should be homogeneous and
      different from any non-air value.
  
  Returns
  -------
  model : array 
    Corrected version of the input model.
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Model_Modify_Set_Air_Vel: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  if mode == 'homo':
    try:
      new_v_air = kwargs['new_v_air']
    except KeyError:
      eprint(this_func + 'Error. For homogeneous air-layer you must provide new_v_air\n')
      quit()
  
  for x in range(len(model)):
   for y in range(len(model[0])):
     for z in range(len(model[0][0])):
       if model[x][y][z] != old_v_air:
         break
     
     if mode == 'homo':
       model[x][y][ :z] = new_v_air
     
     elif mode == 'extrapol':
       model[x][y][ :z] = model[x][y][z]
     
     else:
       eprint(this_func + 'Error. Unknown mode:' + mode + '\n')
       quit()

  #print this_func + 'END'
  return model


# -------------------------------------------------------------------------------


def Model_Modify_Fill_Gaps(model, v_air, proj_name, **kwargs):
  """
  Correct the 3D model in case it has air-velocity cells 
  below the real free surface. This is needed e.g. when 
  air layer was derived from coarser free surface. 
  
  Parameters
  ----------
  model : array 
    Gridded 3D array with Z coordinate being the fastest index 
    (v(x,y,z) = model[x-1][y-1][z-1])
  v_air : float 
    The velocity of the air-layer. It is assumed to be homogeneous and
    different from any non-air value.
  proj_name : str 
    Project name needed to read proj_name + '-FreeSurf.vtr'.
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  Returns
  -------
  model : array 
    Corrected version of the input model.
  
  Notes
  -----
  2D array fs_z stores Z coordinates.
  
  """
  this_func = this_lib + 'Model_Modify_Fill_Gaps: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  proj_path = Kwarg('proj_path', './', kwargs)
  
  fname = Get_Files(proj_path, proj_name + '-FreeSurf.vtr')[0] 
  nx1_fs, nx2_fs, nx3_fs, fs_z = Read_vtr(fname)
  
  for ix in range(len(model)):
    for iy in range(len(model[0])):
      # ITERATE FROM THE BOTTOM UP
      for iz in np.arange(len(model[0][0])-1, -1, -1):
        z_mod = iz + 1 # MAKE z OUT OF GRID-ARRAY INDEX
        z_fs = fs_z[ix][iy][0]
        # FOR NODES BELOW THE FREE SURFACE...
        if z_mod >= z_fs:
          # ...CHECK IF THE VELOCITY IS EQUAL TO THE VELOCITY OF AIR  
          if model[ix][iy][iz] == v_air:
            # IF SO, USE THE VALUE FROM A CELL BELOW
            # NOTE: AS WE ITERATE FROM THE BOTTOM UP IT'S ALWAYS POSSIBLE
            model[ix][iy][iz] = model[ix][iy][iz+1]
            
  #print this_func + 'END'
  return model


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


def Model_Modify_Extend(prefix, model):
  """
  
  Parameters
  ----------
  
  Returns
  -------
  
  Notes
  -----
  """
  this_func = this_lib + 'Model_Modify_Extend: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  
  from lib_fwi_generic import Set_Extra_Nodes
  etop, eleft, efront, ebot, eright, eback = Set_Extra_Nodes(prefix + '-Runfile.key')
  
  print(this_func, 'etop, eleft, efront, ebot, eright, eback', etop, eleft, efront, ebot, eright, eback)
  
  n1, n2, n3 = model.shape
  
  en1 = int(n1 + eleft + eright)
  en2 = int(n2 + efront + eback)
  en3 = int(n3 + etop + ebot)
  
  nm = np.zeros((en1, en2, en3))
  #test = nm[eleft:eleft+n1, efront:efront+n2, etop:etop+n3]
  #print this_func, test.shape
  #print this_func, nm.shape #[1:3][1:4][1:2].shape
  nm[eleft:eleft+n1, efront:efront+n2, etop:etop+n3] = model
  
  return nm
  
  
# -------------------------------------------------------------------------------


def cut_any_model(file_name, xbnds, ybnds, zbnds):
  """
  
  Parameters
  ----------
  proj_name : str
    Project name assumed to be provided 
    as a first argument of the FWI code
    being a prefix of all FWI input files
  nx1, nx2, nx3 : int 
  
  velocity
  
  
  
    
  Returns
  -------
  
  Notes
  -----
  
  """  
  this_func = 'lib_fwi_generic/cut_any_model: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  ext = file_name.split('.')[-1]
  if ext == 'vtr':
    nx1, nx2, nx3, model = Model_Cut_vtr(file_name, xbnds, ybnds, zbnds)
  
  else:
    eprint(this_func + 'ERROR! UNKONWN FILE FORMAT:' + ext + '\n')
    quit()    
  
  #print this_func + 'END'
  return nx1, nx2, nx3, model
  
  
# -------------------------------------------------------------------------------


def Model_Cut_vtr(file_name, xbnds, ybnds, zbnds):
  this_func = 'lib_fwi_generic/Model_Cut_vtr: '
  print(this_func + 'START')  
  
  #from lib_fwi_generic import Read_vtr
  
  nx1, nx2, nx3, model = Read_vtr(file_name)
  print('nx1, nx2, nx3: ', nx1, nx2, nx3)
  print('lens: ', len(model), len(model[0]), len(model[0][0]))
  
  nx1, nx2, nx3, model = Model_Cut(nx1, nx2, nx3, model, xbnds, ybnds, zbnds)
	
  print(this_func + 'END')
  return nx1, nx2, nx3, model  
  
  
# -------------------------------------------------------------------------------


def Model_Cut(nx1, nx2, nx3, model, xbnds, ybnds, zbnds):  
  this_func = 'lib_fwi_generic/Model_Cut: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START') 
  
  #bnds = ''
  if len(xbnds) != 0:
    xmin, xmax = xbnds
    nnx1 = xmax - xmin + 1 # 3 - 1 = 2 BUT NX = 3
  else:
    xmin = -big 
    xmax = big
    nnx1 = len(model) #nx1
    nx1 = nnx1
    
  if len(ybnds) != 0:
    ymin, ymax = ybnds
    nnx2 = ymax - ymin + 1
  else:
    ymin = -big 
    ymax = big
    nnx2 = len(model[0]) #nx2
    nx2 = nnx2
    
  if len(zbnds) != 0:
    zmin, zmax = zbnds
    nnx3 = zmax - zmin + 1
  else:
    zmin = -big 
    zmax = big
    nnx3 = len(model[0][0]) 
    nx3 = nnx3

  nmodel = []
  for x in range(1, nx1 + 1):
    #print 'x = ', x
    xslice = []
    for y in range(1, nx2 + 1):
      #print 'y = ', y
      yslice = []
      for z in range(1, nx3 + 1):
        if (z >= zmin) and (z <= zmax): 
          yslice.append(model[x - 1][y - 1][z - 1])
      if (y >= ymin) and (y <= ymax):
        xslice.append(yslice)
    if (x >= xmin) and (x <= xmax):
      nmodel.append(xslice)
  
  #print this_func + 'END'
  return nnx1, nnx2, nnx3, nmodel     


# -------------------------------------------------------------------------------


#def Model_Resample(...)


# ------------------------------------------------------------------------------- 


def Model_Add_Anomaly(model, chckr, **kwargs):
  """
  Add checkerboard to the model of the same shape.
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  model : array
    Updated.
  
  
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Model_Add_Anomaly: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  ampl_type = Kwarg('ampl_type', 'relative', kwargs)
  
  if ampl_type == 'relative':
    Add_Checker = lambda model, chckr :  model[x][y][z] + chckr[x][y][z] * model[x][y][z]
  
  elif ampl_type == 'absolute':  
    Add_Checker = lambda model, chckr :  model[x][y][z] + chckr[x][y][z]
    
  if model.shape != chckr.shape:
    eprint(this_func + 'Error. Shapes must agree!\n')
    quit()
  
  nx1, nx2, nx3 = model.shape
  
  for x in range(nx1):
    for y in range(nx2):
      for z in range(nx3):
        model[x][y][z] = Add_Checker(model, chckr)
          
  return model


# ------------------------------------------------------------------------------- 


def Model_Modify_Add_Checker(fname, **kwargs):
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
  this_func = this_lib + 'Model_Modify_Add_Checker: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  proj_name = Kwarg('proj_name', 'test', kwargs)
  fname_new = Kwarg('fname_new', fname, kwargs)
  
  ext = Ext(fname_new)
  fname_chckr = fname_new[ :-(len(ext)+1)] + '_chckr.' + ext
  
  model = Model_Read(fname, **kwargs)
  
  # CREATE THE CHECKERBOARD
  dims = model.shape
  XYZ, chckr = Checker_Create(proj_name, dims=dims, **kwargs)
  Model_Save(fname_chckr, chckr)
  
  # MODIFY THE MODEL
  model = Model_Add_Anomaly(model, chckr, **kwargs)
  Model_Save(fname_new, model)

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
# PREPARE SPHERICAL ANOMALY(IES)
# -------------------------------------------------------------------------------


def Spheres_Create(proj, **kwargs):
  """
  
  """
  dims = Kwarg('dims', proj.dims, kwargs)
  
  nx, ny, nz = dims
  xmid, ymid, zmid = int(nx/2), int(ny/2), int(nz/2)
  
  centres = Kwarg('centres', [[xmid, ymid, zmid]], kwargs)
  radii = Kwarg('radii', [10], kwargs)
  amps = Kwarg('amps', [0.05], kwargs)  
  
  #squared_radii = [i**2 for i in radii]
  
  if not (len(radii) == len(centres) and len(amps) == len(centres)):
    raise IOError('All lists must be of the same length.')
  
  Z = np.zeros(dims)
  
  for x in range(nx):
    for y in range(ny):
      for z in range(nz):
        for cen, r, amp in zip(centres, radii, amps):
          x0, y0, z0 = cen
          dist = np.sqrt((x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2)
          if dist <= r:
            Z[x, y, z] += amp
        
  return Z


# -------------------------------------------------------------------------------
# PREPARE A CHECKERBOARD ANOMALY
# -------------------------------------------------------------------------------


def Checker_Create(proj, **kwargs):
  """
  Create a checkerboard pattern of any dimension.
  
  Parameters
  ----------
  dims  : list
    List of lengths of data in each direction
  dx    : list 
    List of space between checkers in each direction;
    =anom_size/2 for dense case
  pads  : list 
    List of 2-tuples of shape 
    [pad_one_end, pad_other_end] 
    for each direction
  sizes : list 
    List of anomaly sizes in each direction.
  ampl  : float 
    Amplitude anomaly as a fraction of 
    a background model
    (one for all checkers).
  mode  : str 
    Shape of checkers: 'rect', 'gauss' etc.
    (one for all checkers)

  Returns
  -------
  nx1, nx2, nx3 : int
    Dimension of the grid as in the geom.
    files.
  points : list 
    List of points [x, y, z]
  
  Notes
  -----
  Works only for regular grids.
  
  """
  this_func = this_lib + 'Checker_Create: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')    

  from lib_math_signal import Dirac_Comb_1D, Filter1D_Gauss_Unnorm, Filter1D_Rect_Unnorm, fwhm2std

  dims = proj.dims
  if verbos > 3:
    print(this_func, 'Inherited dims after', proj.name + ': ', dims)  
  sizes = Kwarg('sizes', 20*np.ones(len(dims)).astype(int), kwargs)
  sizes = sizes.astype(int)
  pads = Kwarg('pads', [[int(size / 2), int(size/ 2)] for size in sizes], kwargs)  
  shape = Kwarg('shape', 'rect', kwargs)
  sparsity = Kwarg('sparsity', 'dense', kwargs)
  ampl = Kwarg('ampl', 0.05, kwargs)
  


  dxs = Checker_Create_Set_Spacing(sizes, shape, sparsity)
  
  if verbos > 2:
    print(this_func, 'sizes', sizes, 'pads', pads, 'dxs', dxs)
  
  # ITERATE OVER ALL DIMENSIONS
  X, A = [], []
  for (nxi, dxi, padi, sizei) in zip(dims, dxs, pads, sizes):
    
    padi = Checker_Create_Adjust_Padding(sizei, padi)
    
    # PREPARE A DIRAC COMB
    Xi, Ai = Dirac_Comb_1D(nxi, dxi, padi, ampl)
    
    
    # PREPARE A SINGLE-CHECKER WAVELET TO CONVOLVE WITH DIRAC COMB
    if shape == 'rect':
      Ai = Filter1D_Rect_Unnorm(Ai, sizei)

    elif shape == 'gauss':
      #FIXME: STH WRONG WITH Indexing
      stdi = fwhm2std(sizei) 
      Ai = Filter1D_Gauss_Unnorm(Ai, stdi)
    
    else:
      raise ValueError('Unknown shape: ' + shape)
    
    if verbos > 0:
      plt.figure()
      plt.plot(Xi, Ai)
      
    X.append(Xi)
    A.append(Ai)
  
  # PREPARE GRIDDED DATA
  XYZ, chckr = Checker_Create_Set_Arrays(dims, X, A, ampl)

  #plt.legend()
  
  #print this_func + 'END'  
  return chckr


# -------------------------------------------------------------------------------


def Checker_Create_Set_Spacing(sizes, shape, sparsity, **kwargs):
  """
  
  
  Parameters
  ----------
  sizes 
  
  shape 
  
  sparsity 
  
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  dxs : list 
    List spacing for each dimension
    respectively.
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Checker_Create_Set_Spacing: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  if sparsity == 'dense': # NOTE: TUNE IT
    gap_factor = 2.2
  else:
    gap_factor = 5 
    
  
  dxs = []
  for size in sizes:
    
    if shape == 'rect':
      dx = size * gap_factor
    
    elif shape == 'gauss':
      dx = size * gap_factor
      
    else:
      raise ValueError('Unknown shape: ' + shape)
  
    dxs.append(dx)
    
  #print this_func + 'END'
  return dxs


# -------------------------------------------------------------------------------


def Checker_Create_Adjust_Padding(sizei, padi, **kwargs):
  """
  
  
  Parameters
  ----------
  
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  padi : list 
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Checker_Create_Adjust_Padding: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  min_pad = int(sizei / 2)
  for j, padij in enumerate(padi):
    if padij < min_pad:
      eprint(this_func + 'Increasing p' + str(j) + ' padding to prevent checker-cutting\n')
      padi[j] = min_pad
    
  #print this_func + 'END'
  return padi


# -------------------------------------------------------------------------------


def Checker_Create_Set_Arrays(dims, X, A, ampl, **kwargs):
  """
  
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  XYZ : array 
  
  chckr : array
  
  
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Checker_Create_Set_Arrays: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  # PREPARE GRIDDED DATA
  if len(dims) == 1:
    chckr = np.ones(len(X[0]))
    for i, A3_i in enumerate(A[0]):
      chckr[i] *= A3_i
      #print i, A3_i, ampl, chckr[i]
    XYZ = np.meshgrid(X[0])
  
  elif len(dims) == 2:
    chckr = np.ones((len(X[0]), len(X[1])))
    for i, A1_i in enumerate(A[0]):
      for j, A3_j in enumerate(A[1]):
        chckr[i][j] *= A1_i * A3_j / float(ampl) # NORMALIZATION!
    
    XYZ = np.meshgrid(X[0], X[1])
  
  elif len(dims) == 3:
    chckr = np.ones((len(X[0]), len(X[1]), len(X[2])))
    for i, A1_i in enumerate(A[0]):
      for j, A2_j in enumerate(A[1]):
        for k, A3_k in enumerate(A[2]):
          chckr[i][j][k] *= A1_i * A2_j * A3_k / float(ampl)**2 # DOUBLE-CHECK IT!
    
    XYZ = np.meshgrid(X[0], X[1], X[2])
  
  else:
    eprint(this_func + 'Error. No. of dimension = ' + dims + ' not yet implemented.\n')
    quit()

  #print this_func + 'END'
  return XYZ, chckr


# ------------------------------------------------------------------------------- 
# GENERIC 
# -------------------------------------------------------------------------------


def Model_Read(fname, **kwargs):
  """
  Read any model.
  
  Parameters
  ----------
  fname : str 
    Name of the file including 
    the path if needed.
  
  Returns
  -------
  model : array 
    Model array.
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Model_Read: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  ext = Ext(fname)
  
  if ext == 'vtr':
    nx1, nx2, nx3, model = Read_vtr(fname)
    
  else:
    eprint(this_func + 'Error. Extension: ' + ext + ' cannot be read yet.\n')
    quit()
  
  #print this_func + 'END'
  return model 


# -------------------------------------------------------------------------------


def Model_Save(fname, model, **kwargs):
  """
  Save any model.
  
  Parameters
  ----------
  fname : str 
    Name of the file including 
    the path if needed.
  
  Returns
  -------
  model : array 
    Model array.
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Model_Save: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  ext = Ext(fname)
  
  if ext == 'vtr':
    nx1, nx2, nx3 = model.shape
    name_core = fname[ :-len('.vtr')]
    Save_vtr(name_core, nx1, nx2, nx3, model)
    
  else:
    eprint(this_func + 'Error. Extension: ' + ext + ' not implemented yet.\n')
    quit()
  
  #print this_func + 'END'
  return model 


# -------------------------------------------------------------------------------

