"""
Author: Kajetan Chrapkiewicz, 2018. 
All rights reserved. Ask for permision writing to K.Chrapkiewicz17@imperial.ac.uk.

This library provides framework procedures for...
At the moment fully implemented is:

1. ...

"""

## MODULES
import numpy as np
from lib_generic import *
from lib_fwi_project_CONST import dims_default


## CONSTS
this_lib = 'lib_fwi_fs.py/'


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
# GENERATE A FREE SURFACE FILE
# -------------------------------------------------------------------------------


def FS_Create(proj_name, **kwargs):
  """
  
  Parameters
  ----------
  proj_name : str
    Project name assumed to be provided 
    as a first argument of the FWI code
    being a prefix of all FWI input files
  
  Returns
  -------
  fs : array 
    Gridded FS.
  
  Notes
  -----
  
  """
  this_func = this_lib + 'FS_Create: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')   
  
  from lib_io_fullwave import Save_vtr
  from lib_io_fullwave_CONST import interfix_fs
 
  pars_gen = Kwarg('pars_gen', {}, kwargs)
  dims = Kwarg('dims', dims_default, pars_gen)
  nx1, nx2, nx3 = dims  
  pars_FS = Kwarg('pars_FS', {}, kwargs)
  fs_type = Kwarg('fs_type', 'sine', pars_FS)
  
  fname = proj_name + interfix_fs + '.vtr'
  fname_core = fname[ :-len('.vtr')]
  
  # SELECT THE MODE
  if fs_type == 'sine':
    fs = FS_Create_Sine(nx1, nx2, **kwargs)
  
  elif fs_type == 'flat':
    fs = FS_Create_Flat(nx1, nx2, **kwargs)
  
  elif fs_type == 'rotated':
    fs = FS_Create_Rotated(nx1, nx2, **kwargs)
    
  else:
    raise NotImplementedError('fs_type: ' + fs_type)
    
  Save_vtr(fname_core, nx1, nx2, nx3, fs, **kwargs)       
    
  #print this_func + 'END'  
  return fs


# -------------------------------------------------------------------------------


def FS_Create_Sine(nx1, nx2, **kwargs):
  """
  Prepare a free surface.
  
  Parameters
  ----------
  
  Returns
  -------
  fs_z : 3D array
    Gridded.
    
  Notes
  -----
  
  """
  this_func = this_lib + 'FS_Create_Sine: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')   
  
  from lib_math_generic import Sine_2D_Old  
  
  try:
    nx3 = kwargs['nx3']
  except KeyError:
    nx3 = 1
  
  try:
    lamb_x = kwargs['lamb_x']
  except KeyError:
    lamb_x = nx1 / 2.
  
  try:
    lamb_y = kwargs['lamb_y']
  except KeyError:
    lamb_y = nx2 / 2.
    
  try:
    phi_x = kwargs['phi_x']
  except KeyError:
    phi_x = np.pi / 2
    
  try:
    phi_y = kwargs['phi_y']
  except KeyError:
    phi_y = np.pi / 2

  try:
    z_min = kwargs['z_min']
  except KeyError:
    z_min = 5
    if verbos > 2:
      print(this_func, 'Default z_min chosen: ', z_min)
  
  try:
    z_max = kwargs['z_max']
  except KeyError:
    z_max = 10
    if verbos > 2:
      print(this_func, 'Default z_max chosen: ', z_max)
  
  
  # CALCULATE VALUES
  fs_z = np.zeros((nx1, nx2, nx3))
  
  x0 = nx1 / 2. - 1 # FIXME: D-C
  y0 = nx2 / 2. - 1
  
  for x in range(nx1):
    for y in range(nx2):
      for z in range(nx3):
        fs_z[x][y][z] = Sine_2D_Old(x-x0, y-y0, lamb_x, lamb_y, phi_x, phi_y, z_min, z_max)  
  
  #print this_func + 'END'  
  return fs_z


# -------------------------------------------------------------------------------


def FS_Create_Flat(nx1, nx2, **kwargs):
  """
  Prepare a free surface.
  
  Parameters
  ----------
  
  Returns
  -------
  model : 3D array
    Gridded.
    
  Notes
  -----
  
  """
  this_func = this_lib + 'FS_Create_Flat: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')   

  nx3 = Kwarg('nx3', 1, kwargs)  
  pars_FS = Kwarg('pars_FS', {}, kwargs)
  z_fs = Kwarg('z_fs', 5, pars_FS)
  
  model = np.ones((nx1, nx2, nx3)) * z_fs
  
  #print this_func + 'END'   
  return model


# -------------------------------------------------------------------------------


def FS_Create_Rotated(nx1, nx2, **kwargs):
  """
  Prepare a free surface.
  
  Parameters
  ----------
  
  Returns
  -------
  model : 3D array
    Gridded.
    
  Notes
  -----
  
  """
  this_func = this_lib + 'FS_Create_Rotated: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')   

  nx3 = Kwarg('nx3', 1, kwargs)  
  theta = Kwarg('theta', 30, kwargs)
  tilt_axis = Kwarg('tilt_axis', 'x', kwargs)
  z0 = Kwarg('z0', 5, kwargs) # z0 = a*0 + b = b 
  
  # CALC. THE SLOPE OF LINEAR FUNCTION
  a = np.tan(theta * np.pi / 180.) # CHECKED
  b = z0
  model = np.zeros((nx1, nx2, nx3))
  
  for x in range(nx1):
    for y in range(nx2):
      for z in range(nx3):
        model[x][y][z] = a * x + b  
    
  
  #print this_func + 'END'   
  return model


# -------------------------------------------------------------------------------
# HANDLE GHOST DATA
# -------------------------------------------------------------------------------


def Index_Ghost_From_Coords(ghosts, **kwargs):
  """
  Return an array-index/indices of the ghost(s) 
  of given coordinates.
  
  Parameters
  ----------
  ghosts : list 
  
  **kwargs : keyword arguments, optional
    x
    y
    z
  
  
  Returns
  -------
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Index_Ghost_From_Coords: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')   
  
  # CHECK KEYWORD ARGUMENTS
  try:
    x = [kwargs['x']]
  except KeyError:
    x = np.arange(np.min([g[0] for g in ghosts]), np.max([g[0] for g in ghosts]) + 1)

  try:
    y = [kwargs['y']]
  except KeyError:
    y = np.arange(np.min([g[1] for g in ghosts]), np.max([g[1] for g in ghosts]) + 1)

  try:
    z = [kwargs['z']]
  except KeyError:
    z = np.arange(np.min([g[2] for g in ghosts]), np.max([g[2] for g in ghosts]) + 1)
    
  print(this_func, 'x', x)
  print(this_func, 'y', y)
  print(this_func, 'z', z)
  
  indices = []
  for ix in x:
    for iy in y:
      for iz in z:
        xyz = [ix, iy, iz]
        print(this_func, 'xyz', xyz)
        index = Index_Ghost_XYZ(ghosts, xyz)
        if index > 0:
          indices.append(index)
  
  #print this_func + 'END'  
  return indices


# -------------------------------------------------------------------------------


def Index_Ghost_XYZ(ghosts, xyz, **kwargs):
  """
  Return an array-index/indices of the ghost(s) 
  of given coordinates.
  
  Parameters
  ----------
  ghosts : list 
  
  xyz : vector
    x, y, z = xyz
  
  Returns
  -------
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Index_Ghost_XYZ: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START') 
  
  x, y, z = xyz
  
  found = False
  indx = -1
  
  for i, g in enumerate(ghosts):
    xg, yg, zg = g[0], g[1], g[2]
    if x == xg and y == yg and z == zg:
      found = True
      indx = i
      break
  
  if not found:
    eprint(this_func + 'Ghost: (' + str(x) + ',' + str(y) + ',' + str(z) + ') not found.\n')  
  
  #print this_func + 'END'  
  return indx


# -------------------------------------------------------------------------------


def Read_GhostData_File_txt(fname, **kwargs):
  """
  Read all the information contained 
  in a GhostData.txt file.
  
  Parameters
  ----------
  fname : str    
    File name. It should include .txt extension. 
    It can include path if needed.
  
  Returns
  -------
  ghosts : list
  intersects : list
  ficts : list
  auxs : list
  weights : list
  
  Notes
  -----
  It works! 
  
  """
  this_func = this_lib + 'Read_GhostData_File_txt: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')    
  
  ct = Read_File(fname)
  header = ct[0]
  data = ct[1: ]
  fict_no, auxs_no = [int(i) for i in header]
  
  ghosts_no, intersects_no = 1, 1
  ghosts, intersects, ficts, auxs, weights = [], [], [], [], []
  
  j = 0
  while j < len(data):
    ghosts.append(data[j])
    j += 1
    intersects.append(data[j])
    j += 1
    
    fcs = []
    for f in range(fict_no):
      fcs.append(data[j])
      j += 1
    ficts.append(fcs)
    
    fcs_auxs = []
    for f in range(fict_no):
      f_auxs = []
      for ax in range(auxs_no):
        f_auxs_x = []
        for ay in range(auxs_no):
          f_auxs_x.append(data[j])
          j += 1
        f_auxs.append(f_auxs_x)
      fcs_auxs.append(f_auxs)  
    auxs.append(fcs_auxs)
    
    fcs_weights = []
    for f in range(fict_no):
      f_weights = []
      for w in range(2 * auxs_no):
        f_weights.append(data[j])
        j += 1
      fcs_weights.append(f_weights)
    weights.append(fcs_weights)
   
  ghosts = [[float(i) for i in j] for j in ghosts]
  intersects = [[float(i) for i in j] for j in intersects]
  ficts = [[[float(i) for i in j] for j in k] for k  in ficts]
  auxs = [[[[[float(i) for i in j] for j in k] for k in l] for l in m] for m in auxs]  
  weights = [[[[float(i) for i in j] for j in k] for k in l] for l in weights]
  
  #print this_func + 'END' 
  return ghosts, intersects, ficts, auxs, weights


# -------------------------------------------------------------------------------
# GENERATE SR FOR FREE-SURFACE TESTS
# -------------------------------------------------------------------------------


def prepare_flat_geometry(xm, ym, zm, dx, fs_depth, rs_depth_ratio, sou_depth, ns, rec_dxs, steps):
  from lib_math_num import Metres2nodes # FS BEFORE ROTATION FIXME ONLY ONE LINE OF RECEIVERS AT A TIME (IN  y= ny/2)  
  
  #print xm, ym, zm, dx, fs_depth, rs_depth_ratio, sou_depth, ns, rec_dxs, steps
  
  nx = Metres2nodes(1000 * xm, dx)
  if ym == 0:
    ny = 1
  else:
    ny = Metres2nodes(1000 * ym, dx)
  nz = Metres2nodes(1000 * zm, dx)
  
  rec_depth = sou_depth * rs_depth_ratio
  
  if ny == 1:
    y_coord = 1
  else:
    y_coord = ny / 2 # INT
  
  x_n, y_n = ns # NO. OF RECEIVERS
  rec_dx = rec_dxs[0]
  xleft = (nx - (2 * x_n * rec_dx)) / 2 # 2* CAUSE x_n DESCRIBES ONLY 1 SIDE
  xright = xleft 
  if xleft < 1:
    eprint('ERROR! wavefor2d NEEDS AT LEAST 1-NODE PADDING\n')
    quit()
  
  x_step, y_step = steps
        
  # SOURCE & RECEIVERS
  sou, rec = prepare_flat_source_n_rec(fs_depth, rec_depth, sou_depth, rec_dx, x_n, y_coord, xleft)
  
  if nz < sou[2] or nz < max([r[2] for r in rec]):
    eprint('ERROR! nz TOO SMALL EVEN FOR FLAT GEOMETRY, nz:' + str(nz) + '\n')
    eprint('sou[1]: ' + str(sou[1]) + '\n')
    eprint('max r[2]: ' + str(max([r[2] for r in rec])) + '\n')
    quit()
  
  # FREE SURFACE
  fs_x = np.arange(1, rec[-1][0] + xright + 1)
  x_start = 1
  x_end = nx
  y_start = 1
  y_end = ny
  fs = prepare_flat_fs_grd(x_start, nx, x_step, y_start, ny, y_step, fs_depth)
  
  print('Model dimensions: ', nx, ny, nz)
  return nx, ny, nz, rec, fs, sou, rec_depth


# -------------------------------------------------------------------------------


def prepare_flat_source_n_rec(fs_depth, rec_depth, sou_depth, rec_dx, x_n, y_coord, xleft):
  rec_x = np.arange(1, (2 * x_n + 1) * rec_dx + 1, rec_dx)
  rec_x += xleft
  #rec_y = np.arange(1, (2 * y_n + 1) * rec_dy + 1, rec_dy)
  #rec_y += yleft
  rec_y = np.zeros(len(rec_x)) + y_coord # IN THE MIDDLE OF Y AXIS
  #print np.meshgrid(rec_x, rec_y)
  rec_z = np.zeros(len(rec_x)) + rec_depth + fs_depth # OTHERWISE THEY AIN'T IN-NODE 
  
  if sou_depth < 0:
    sou_depth = 2 * rec_depth + fs_depth
  else:
    sou_depth += fs_depth
    
  sou = (rec_x[len(rec_x) / 2], y_coord, sou_depth) # INDICES START FROM 0
  rec = list(zip(rec_x, rec_y, rec_z))
  return sou, rec


# -------------------------------------------------------------------------------


def save_fs_3d_pgy(fname, nx1, nx2, nx3, X, Y, Z):
  Z = np.ravel(Z) # SAVE 3D DATA TO .pgy FILE
  X = [float(i) for i in X]
  Y = [float(i) for i in Y]
  Z = [float(i) for i in Z]
  
  
  n = len(X)# nx1 * nx2 # NO. OF POINTS
  sn   = '{:10}'.format(n)
  snx3 = '{:15}'.format(nx3)
  snx2 = '{:15}'.format(nx2)
  snx1 = '{:15}'.format(nx1)
  
  f = open(fname + '.pgy', 'w')
  f.write(sn + snx3 + snx2 + snx1 + '\n')
  for i in range(len(X)):
    si   = '{:10}'.format(i + 1)
    sz   = '{:15.8f}'.format(Z[i])
    sy   = '{:15.8f}'.format(Y[i])
    sx   = '{:15.8f}'.format(X[i])  
    #print si, sz, sy, sx
    f.write(si + sz + sy + sx + '\n')
    i += 1
  f.close()


# -------------------------------------------------------------------------------




##############################################################################################









# TOPO 
def slice_xyz(prefix, X, Y, Z, x_min, x_max, y_min, y_max):
  f = open(prefix + str(x_min) + '-' +  str(x_max) + '-' + str(y_min) + '-' + str(y_max)+ '.xyz', 'w')
  NX, NY, NZ = [], [], [] 
  for i in range(len(Z)):
    
    x = X[i]
    y = Y[i]
    z = Z[i]

    if ((y >= y_min) and (y <= y_max) and (x >= x_min) and (x <= x_max)):
      NX.append(x)
      NY.append(y)
      NZ.append(z)
      f.write(str("{:10.4f}".format(x)) + '  ' + str("{:10.4f}".format(y)) + '   ' + str("{:10.4f}".format(z)) + '\n')
  f.close()
  return NX, NY, NZ    





#
def prepare_auxs_interp_geometry(ficts, auxs):
  print('lib_fwi_fs/prepare_auxs_interp_geometry: START')
  from lib_math import Dist
  
  xx = []
  for i in range(len(ficts)):
    f = np.array(ficts[i][ :3])
    x_f_auxs = []
    for j, aux in enumerate(auxs[i]):
      a = np.array(aux[ :3])
      if j < len(auxs[i]) / 2: # FIRST HALF OF AUXS ARE LEFT TO THE Y-AXIS
        x_f_auxs.append(-Dist(f, a))
      else:
        x_f_auxs.append(Dist(f, a))
    xx.append(x_f_auxs + [0.]) # 0. IS FOR FICT

  print('lib_fwi_fs/prepare_auxs_interp_geometry: END') 
  return xx

#
def prepare_fict_interp_geometry(ghost, isect, ficts):
  print('lib_fwi_fs/prepare_fict_interp_geometry: START') 
  from lib_math import Dist
  
  g = np.array(ghost)[ :3]
  i = np.array(isect)
  dr = i - g 
  rpoint = i + dr

  f_x = []
  for fic in ficts:
    f = np.array(fic[ :3])
    f_x.append(dist(f, i))
  i_x = Dist(i, i)
  g_x = -Dist(g, i) # NOTE THE SIGN
  r_x = Dist(rpoint, i)
  
  print('lib_fwi_fs/prepare_fict_interp_geometry: END')  
  return [g_x, i_x, r_x] + f_x 


  






# PREPARE 3D GRIDDED FREE SURFACE IN X DIRECTION
def prepare_X_sin_fs_grd(fs, z_min, amp, lambd, nx, ny):
  X, Y, Z = fs
  X = [float(x) for x in X] # ABSOLUTELY CRUCIAL, OTHERWISE 5-POINT PYRAMID!
  Z = Z.transpose() # TO MODIFY IN X DIRECTION AND KEEP CONSTANT IN Y DIRECTION
  for i, y in enumerate(Y):
    Z[i] = [(z_min + amp) + amp * np.cos(2 * np.pi * (x - nx / 2 - 1) / lambd) * np.cos(2 * np.pi * (y - ny / 2 - 1) / lambd) for x in X]
  Z = Z.transpose() 
  #plt.plot(Z)
  #plt.show()
  # CONVERT BACK TO SHAPE COMPATIBLE WITH OTHER READ/WRITE FUNCTION
  return X, Y, Z

# PREPARE 3D GRIDDED FREE SURFACE GAUSSIAN IN X DIRECTION
def prepare_X_gauss_fs_grd(fs, z_min, z_max, mu, sigma):
  from lib_math_signal import Gaussian_1D_Unnormed
  X, Y, Z = fs
  X = [float(x) for x in X] # ABSOLUTELY CRUCIAL, OTHERWISE 5-POINT PYRAMID!
  Z = Z.transpose() # TO MODIFY IN X DIRECTION AND KEEP CONSTANT IN Y DIRECTION
  for i in range(len(Y)):
    Z[i] = [z_max + (z_min - z_max) * Gaussian_1D_Unnormed(x, mu, sigma) for x in X]
  Z = Z.transpose() # CONVERT BACK TO SHAPE COMPATIBLE WITH OTHER READ/WRITE FUNCTION
  return X, Y, Z
  
# PREPARE 3D GRIDDED FREE SURFACE TILTED IN X DIRECTION
def prepare_X_tilted_fs_grd(fs, alpha, y0):
  X, Y, Z = fs
  Z = Z.transpose() # TO MODIFY IN X DIRECTION AND KEEP CONSTANT IN Y DIRECTION
  for i in range(len(Y)):
    Z[i] = Straight_Line(X, alpha, y0)
  Z = Z.transpose() # CONVERT BACK TO SHAPE COMPATIBLE WITH OTHER READ/WRITE FUNCTION
  return X, Y, Z
  




  
  
#
def rotate_flat_geometry(project, theta, fs_flat, fs_depth, rec_depth, sou, rec, nz):
  from lib_fwi_fs import plot_flat_n_rotated
  from lib_math import Dist
  from lib_generic_CONST import epsi
  
  X, Y, Z = fs_flat
  Z = np.ravel(Z)
  if len(Y) == 1:
    Y = np.zeros(len(X)) + Y[0]
  else:
    eprint('ERROR! 3D NOT YET HERE\n')
    quit()
  fs = list(zip(X, Y, Z))

  #rot_axis = prepare_rot_axis(rec, int(rec_depth), int(fs_depth))
  rot_axis = (sou[0], sou[2]) # FIXME: 2D ONLY
  #print rot_axis, theta
  #print fs
  theta, rot_fs, rot_sou, rot_rec = rotate_all(fs, sou, rec, rot_axis, theta)
  #print 'Rotated receivers:\n', rot_rec

  
  rot_fs = extend_fs(fs, rot_fs) # STRETCH!
  
  rot_fs, rot_rec, rot_sou = shift_all(rot_fs, rot_rec, rot_sou, fs_depth)
  
  if abs(Dist(rot_sou, rot_rec[0]) - Dist(sou, rec[0])) > epsi:
    eprint('ERROR. DISTANCE |SR| CHANGED THROUGH ROTATION\n')
    print(Dist(rot_sou, rot_rec[0]), Dist(sou, rec[0]))
    quit()
  
  return rot_rec, rot_sou, rot_fs

# DETERMINE AXIS AND ANGLE OF ROTATION
def prepare_rot_axis(rec, rec_depth, fs_depth):
  from lib_math_num import rad2deg
  centre = 'first'#'middle'
  # MUST COINCIDE WITH 1 OF THE RECEIVERS (FACTUAL OR 'EXTRAPOLATED' + i*rec_dx)
  if centre == 'first':
    rot_axis = (1, int(fs_depth)) 
    #rot_axis = (rec[0][0], int(rec_depth) + int(fs_depth))   
  elif centre == 'middle': # EQUAL DISTANCE FROM BOTH MARGINS
    rot_axis = (rec[len(rec) / 2][0], int(rec_depth) + int(fs_depth))
  else:
    eprint('ERROR. UNKNOWN CENTRE OF ROTATION\n')
    quit()
  print('Rotation axis: ', rot_axis)
  
  return rot_axis

# ROTATE FS, S AND Rs
def rotate_all(fs, sou, rec, rot_axis, theta):
  from lib_math import Rotate_2D , Dist
  rot_rec = [Rotate_2D(rot_axis, r, theta) for r in rec] 
  rot_fs = [Rotate_2D(rot_axis, f, theta) for f in fs] 
  #print rot_fs
  rot_sou = Rotate_2D(rot_axis, sou, theta)
  
  #print dist(rot_rec[0], rot_sou)
  #if dist(rec[len(rec) / 2], sou) != dist(rot_rec[len(rec) / 2], rot_sou):
  #  eprint('ERROR. DISTANCE |SR| CHANGED THROUGH ROTATION\n')
  #  print (dist(fs[0], rec[len(rec) / 2]), ' vs ',
  #        dist(rot_fs[0], rot_rec[len(rec) / 2]))
  #  quit()
  return theta, rot_fs, rot_sou, rot_rec

# STRETCH FS TO MODEL BOUNDARIES MAKING USE OF ITS LINEARITY  
def extend_fs(fs, rot_fs):
  from lib_math import Find_Linear_Function
  if len(rot_fs[0]) > 2:
    eprint('WARNING! 3D->2D TO FIX\n')
    rot_fs = [[r[0], r[2]] for r in rot_fs]
    y = 1.0
  
  a, b = Find_Linear_Function(rot_fs[0], rot_fs[1])
  new_fs = []
  for i in range(len(fs)):
    x = i + 1
    z = a * x + b
    new_fs.append((x, y, z))
  return new_fs

# PUT EVERYTHING BACK INTO MODEL DOMAIN
def shift_all(rot_fs, rot_rec, rot_sou, fs_depth):
  rot_rec_x, rot_rec_y, rot_rec_z = [r[0] for r in rot_rec], [r[1] for r in rot_rec], [r[2] for r in rot_rec]
  rot_fs_x, rot_fs_y, rot_fs_z = [r[0] for r in rot_fs], [r[1] for r in rot_fs], [r[2] for r in rot_fs]
  
  #shift =  np.ceil(-min(rot_fs_z) + fs_depth)
  shift = -min(rot_fs_z) + fs_depth
  
  rot_fs_z = [z + shift  for z in rot_fs_z]
  rot_rec_z = [z + shift  for z in rot_rec_z]
  rot_sou = np.add(rot_sou, (0, 0, shift)) # OTHERWISE IT WOULD BE APPENDED
  
  rot_fs = list(zip(rot_fs_x, rot_fs_y, rot_fs_z))
  rot_rec = list(zip(rot_rec_x, rot_rec_y, rot_rec_z))
  
  return rot_fs, rot_rec, rot_sou








# SET NO. OF AUXILIARY NODES PER FICTITIOUS POINT FOR 2D/3D
def set_auxs_per_fict(auxs_no, ny):
  if ny == 1:
    return auxs_no
  elif ny > 1:
    return auxs_no ** 2
  else:
    eprint('ERROR! set_auxs_per_fict: WEIRD ny\n')
    quit()

 
  
  
  
  
  
  
# CUT FS TO LEAVE ONLY AREA WITHIN radius FROM (x, y) <- INDICES (MAKE IT MORE GENERAL)
def cut_fs(X, Y, Z, x, y, radius):
  xmin = x - radius + 1
  if xmin < 0:
    xmin = 0  
  xmax = x + radius 

  
  X = X[xmin : xmax]
  
  if len(Y) > 1:
    ymin = y - radius + 1
    if ymin < 0:
      ymin = 0    
    ymax = y + radius   
    Y = Y[ymin : ymax]
    Z = [i[ymin : ymax] for i in Z[xmin : xmax]]  
  else:
    Z = Z[xmin : xmax]
  return X, Y, Z

# READ FS GRID FILE
def read_fs_grd(fname):
  this_func = this_lib + 'read_fs_grd: '
  print(this_func + 'START')
  
  content = Read_File(fname)
  try:
    x_start, nx, x_step, y_start, ny, y_step = content[0]
  except ValueError: 
    eprint(this_func + 'Probably header line is missing.\n')
    #print len(content)
    #print len(content[0])
    print(('First line of the file: ' + line1 + '\n'))
    #eprint('Trying to figure the header values another way\n') # FIXME: IT NEEDS PASSING MORE ARGUMENTS => MODIFICATION OF THE CALLERS
    quit()
  
  # CAST STRINGS TO CORRECT TYPES
  x_start, nx, y_start, ny = [int(i) for i in [x_start, nx, y_start, ny]]
  x_step, y_step = [float(i) for i in [x_step, y_step]]
  ##
  X = np.arange(x_start, x_start + x_step * nx, x_step)
  if ny > 1:
    Y = np.arange(y_start, y_start + y_step * ny, y_step)
  else:
    Y = [1]
  data = content[1: ]
  data = np.ravel(data) # FLATTEN THE LIST
  data = [float(i) for i in data]
  
  x_slices = Get_Chunks(data, ny)
  Z = []
  for x_slice in x_slices:
    Z.append(x_slice)
  return X, Y, Z
  
# SAVE 2D FS TO .pgy FILE 
def save_fs_2d_pgy(fname, nx1, nx2, nx3, X, Z):
  Z = np.ravel(Z)
  X = [float(i) for i in X]
  Z = [float(i) for i in Z]
  
  n = len(X)# nx1 * nx2 # NO. OF POINTS
  sn   = '{:10}'.format(n)
  snx3 = '{:15}'.format(nx3)
  snx2 = '{:15}'.format(nx2)
  snx1 = '{:15}'.format(nx1)
  
  f = open(fname + '.pgy', 'w')
  f.write(sn + snx3 + snx2 + snx1 + '\n')
  for i in range(len(X)):
    si   = '{:10}'.format(i + 1)
    sz   = '{:15.8f}'.format(Z[i])
    sy   = '{:15.8f}'.format(1.0)
    sx   = '{:15.8f}'.format(X[i])  
    #print si, sz, sy, sx
    f.write(si + sz + sy + sx + '\n')
    i += 1
  f.close()

# PREPARE .pgy FREE SURFACE FILE
def prepare_flat_fs_2d_pgy(fname, nx1, nx2, nx3):
  n = nx1 * nx2 # NO. OF POINTS
  sn   = '{:10}'.format(n)
  snx3 = '{:15}'.format(nx3)
  snx2 = '{:15}'.format(nx2)
  snx1 = '{:15}'.format(nx1)
  
  f = open(fname, 'w')
  f.write(sn + snx3 + snx2 + snx1 + '\n')
  i = 1
  for y in range(1, nx2 + 1):
    for x in range(1, nx1 + 1):
      si   = '{:10}'.format(i)
      sz   = '{:15.8f}'.format(1.)
      sy   = '{:15.8f}'.format(y)
      sx   = '{:15.8f}'.format(x)  
      #print si, sz, sy, sx
      f.write(si + sz + sy + sx + '\n')
      i += 1
  f.close()

# READ .pgy FREE SURFACE FILE
def read_fs_2d_pgy(fname):
  nx1, nx2, nx3, x, y, z = Read_pgy(fname)
  yslices = Get_Chunks(z, nx1)
  fs = []
  for yslice in yslices:
    fs.append(yslice)
  fs = np.array(fs)
  fs = np.transpose(fs)
  return nx1, nx2, nx3, fs

















## RUBBISH


# RETURN AN INDEX OF A GHOST OF GIVEN COORDINATES
def pick_ghost_by_x(ghosts, x):
  print('lib_fwi_fs/pick_ghost_by_x: START')  
  
  i_list = []
  
  for i, g in enumerate(ghosts):
    xg, yg, zg = g[0], g[1], g[2]
    if x == xg:
      i_list.append(i)
  
  print('lib_fwi_fs/pick_ghost_by_x: END')  
  return i_list 
  

  

# READ FS GRID FILE
def read_fs_grd_old(fname):
  content = Read_File(fname)
  x_start, nx, x_step, y_start, ny, y_step = content[0]
  # CAST STRINGS TO CORRECT TYPES
  x_start, nx, y_start, ny = [int(i) for i in [x_start, nx, y_start, ny]]
  x_step, y_step = [float(i) for i in [x_step, y_step]]
  ##
  X = np.arange(x_start, x_start + x_step * nx, x_step)
  if ny > 1:
    Y = np.arange(y_start, y_start + y_step * ny, y_step)
  else:
    Y = [1]
  data = content[1: ]
  #print data
  data = np.ravel(data)
  
  data = [float(i) for i in data]
  #print 'after'
  #print data
  
  if ny > 1:
    y_slices = Get_Chunks(data, nx)
    Z = []
    for y_slice in y_slices:
      Z.append(y_slice)
  elif ny == 1:
    Z = data
  else:
    eprint('ERROR. STRANGE VALUE OF ny')
  
  print(this_func + 'END')
  return X, Y, Z



# PREPARE 3D FLAT (HORIZONTAL) FREE SURFACE GRID FILE
def prepare_flat_fs_grd(x_start, nx, x_step, y_start, ny, y_step, z_value):
  X = np.arange(x_start, x_step * nx + 1, x_step)
  Y = np.arange(y_start, y_step * ny + 1, y_step) 
  Z = np.ones((len(X), len(Y))) * z_value
  return X, Y, Z
