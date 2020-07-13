"""
Author: Kajetan Chrapkiewicz, 2018. 
All rights reserved. Ask for permision writing to K.Chrapkiewicz17@imperial.ac.uk.

This library provides framework procedures for...
At the moment fully implemented is:

1. ...

"""

## MODULES
import matplotlib.pyplot as plt
import numpy as np
from lib_generic import *
from lib_io_fullwave import *
from lib_generic_PLOTT import *
##


## CONSTS
this_lib = 'lib_fwi_fs_PLOTT.py/'

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
# GHOST DATA
# -------------------------------------------------------------------------------


def Plot_Ghost_Data(fname, **kwargs):
  """
  Plot various subsets of ghost data.
  
  Parameters
  ----------
  fname : str 
    It must include extension proceeded by a full stop (e.g. 'file.txt').
    It can include path if needed (e.g. '/home/file.txt'). 
  **kwargs : keyword arguments, optional
      Just passing it down.
  
  Returns
  -------
  0
  
  Notes
  -----
  
  """  
  this_func = this_func = this_lib + 'Plot_Ghost_Data: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  from lib_fwi_fs import Read_GhostData_File_txt
  from lib_generic_PLOTT import Plot_Gridded_Data_From_File, Plot_Gridded_Data
  
  # READ THE DATA-FILE
  ghosts, intersects, ficts, auxs, weights = Read_GhostData_File_txt(fname)
  
  # PLOTT GHOSTS
  Plot_Ghosts(ghosts, **kwargs)
  
  # PLOTT GIFA
  #Plot_GIFA(ghosts, intersects, ficts, auxs, **kwargs)
  
  # PLOTT FICT. INTERPOLATION 
  
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Ghosts(ghosts, **kwargs):
  """
  Plot ghost nodes.
  
  Parameters
  ----------
  ghosts : list 
    Raw list output from Read_GhostData_File_txt.
  **kwargs : keyword arguments, optional
      Just passing it down.
  
  Returns
  -------
  0
  
  Notes
  -----
  
  """  
  this_func = this_func = this_lib + 'Plot_Ghosts: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  import matplotlib.cm as cm
  from lib_generic_PLOTT import Plot_Gridded_Data, Plot_Points
  
  # CHECK KEYWORD ARGUMENTS
  try:
    fs_z = kwargs['fs_z']
  except KeyError:
    #eprint(this_func + 'No FS provided - plotting only wavefield.\n')
    fs_z = []  
  
  # READ GHOSTS
  ghosts_xyz = [g[ :3] for g in ghosts] # NOTE: REDUNDANT
  ghosts_xyz_all_lvls = []
  nglvls = 3
  for lvl in range(1, nglvls + 1):
    ghosts_xyz_lvl = []
    for g in ghosts:
      if g[4] == lvl:
        ghosts_xyz_lvl.append(g[ :3])
    ghosts_xyz_all_lvls.append(ghosts_xyz_lvl)
  
  # 2D SLICE
  plt.figure()
  if len(fs_z) > 0: # PLOTT FS IF PROVIDED
    Plot_Gridded_Data(fs_z, data_type='FS', plot_type='slice', **kwargs)  
  
  # BECAUSE GRIDDED DATA STRUCTURES (FS, MODEL, ETC.) ARE COUNTED FROM 0
  try:
    kwargs['xslice'] += 1
  except KeyError:
    pass

  try:
    kwargs['yslice'] += 1
  except KeyError:
    pass

  try:
    kwargs['zslice'] += 1
  except KeyError:
    pass  
  
  # 2D SLICE
  #colors = iter(cm.rainbow(np.linspace(0, 1, nglvls)))
  #for lvl in range(nglvls):
  #  Plot_Points(ghosts_xyz_all_lvls[lvl], plot_type='slice', c=next(colors), **kwargs)
     
  # ALL IN 3D
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  colors = iter(cm.rainbow(np.linspace(0, 1, nglvls)))
  for lvl in range(nglvls):
    Plot_Points(ghosts_xyz_all_lvls[lvl], plot_type='scatter', ax=ax, c=next(colors), alpha=.5, azim=45, elev=40)
  plt.gca().invert_zaxis()
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_GIFA(ghosts, intersects, ficts, auxs, **kwargs):
  """
  Plot a ghost node and its intersection, 
  fict. points, aux. nodes.
  
  Parameters
  ----------
  ghosts : list 
    Raw list output from Read_GhostData_File_txt.
  intersects : list
    Raw list output from Read_GhostData_File_txt.
  ficts : list
    Raw list output from Read_GhostData_File_txt.
  auxs : list
    Raw list output from Read_GhostData_File_txt.
  **kwargs : keyword arguments, optional
      Just passing it down.
  
  Returns
  -------
  0
  
  Notes
  -----
  
  """  
  this_func = this_func = this_lib + 'Plot_GIFA: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  from lib_generic_CONST import big
  from lib_fwi_fs import Index_Ghost_From_Coords
  import matplotlib.cm as cm
  from lib_generic_PLOTT import Plot_Gridded_Data, Plot_Points
  
  # CHECK KEYWORD ARGUMENTS
  try:
    fs_z = kwargs['fs_z']
  except KeyError:
    eprint(this_func + 'No FS provided - plotting only wavefield.\n')
    fs_z = []  
  
  xmin, ymin = big, big
  xmax, ymax = -big, -big
  for g in ghosts:
    if g[0] < xmin:
      xmin = g[0]
    if g[0] > xmax:
      xmax = g[0]

    if g[1] < ymin:
      ymin = g[1]
    if g[1] > ymax:
      ymax = g[1]  
  
  xmid = int((xmin + xmax) / 2)
  ymid = int((ymin + ymax) / 2)
  
  indices = Index_Ghost_From_Coords(ghosts, x=xmid, y=ymid)
  #indices += Index_Ghost_From_Coords(ghosts, x=xmin, y=ymid)
  #indices += Index_Ghost_From_Coords(ghosts, x=xmid, y=ymin)
  #indices += Index_Ghost_From_Coords(ghosts, x=xmax, y=ymid)
  #indices += Index_Ghost_From_Coords(ghosts, x=xmid, y=ymax)
  
  j = -1
  for g, i, f, a  in zip(ghosts, intersects, ficts, auxs):
    j += 1
    if not j in indices:
      continue 
    
    g = g[ :3]
    x, y, z = g
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.title('Ghost (' + str(x) + ',' + str(y) + ','  + str(z) + ')')  
    
    azim = 45
    elev = 45
    Plot_Points([g], plot_type='scatter', ax=ax, c='b', s=40, azim=azim, elev=elev)
    Plot_Points([i], plot_type='scatter', ax=ax, c='r', s=40, azim=azim, elev=elev)
    #normal = np.array(g[:3]) - np.array(i[:3])
    #normal = Normed(normal)
    
    Plot_Points([ff[ :3] for ff in f], plot_type='scatter', ax=ax, c='green', s=40, azim=azim, elev=elev)
    points = []
    for aa in a:
      for aaa in aa:
        for aaaa in aaa:
          points.append(aaaa[ :3])
    Plot_Points(points, plot_type='scatter', ax=ax, c='gray', s=40, alpha=.5, azim=azim, elev=elev)
    
    plt.gca().invert_zaxis()
    
    if len(fs_z) > 0: # PLOTT FS IF PROVIDED
      Plot_Gridded_Data(fs_z, data_type='FS', plot_type='surf', ax=plt.gca(), alpha=0.1, **kwargs)  
   
    #plt.gca().auto_scale_xyz([10, 20], [10, 20], [15, 5])
    #break
     
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------

















##########################################################################################


# MAKE INTERACTIVE 3D PLOTT OF GHOSTS GIVEN COORDINATES 1D ARRAYS (X, Y) AND 2D DATA ARRAY (Z)   
def plot_ghosts_3D(X_fs, Y_fs, Z_fs, X_gh, Y_gh, Z_gh):
  XX_fs, YY_fs = np.meshgrid(X_fs, Y_fs)
  Z_fs = np.array(Z_fs)
  Z_fs = Z_fs.reshape(XX_fs.shape)
  
  XX_gh, YY_gh = np.meshgrid(X_gh, Y_gh)
  Z_gh = np.array(Z_gh)
  Z_gh = Z_gh.reshape(XX_gh.shape)
  
  fig = plt.figure()
  ax = fig.add_subplot(111, projection = '3D')
  
  minn = min(Z_gh.min(axis = 1))
  if minn < 0:
    minn = 0
  maxx = max(Z_gh.max(axis = 1))
  
  plt.xlim(min(X_gh), max(X_gh))
  plt.ylim(min(Y_gh), max(Y_gh))
  ax.set_zlim(minn, maxx)
  
  #ax.set_xlabel("X [nodes]")
  #ax.set_ylabel("Y [nodes]")
  #ax.set_zlabel("Z [nodes]")
  #ax.invert_zaxis() # FLIP Z AXIS
  ax.scatter(XX_gh, YY_gh, Z_gh)
  ax.plot_surface(XX_fs, YY_fs, Z_fs, color = 'r', alpha = .1) # 'dodgerblue'





# FOR ALL GHOSTS OF glvl LEVEL, PLOTT ALL AUXILIARY POINTS (I, F, A)
def show_all_ghosts(fname):
  from lib_fwi_fs import Read_GhostData_File_txt
  ghosts, intersects, ficts, auxs = Read_GhostData_File_txt(fname)
  
  for i in range(len(ghosts)):
    x, y = ghosts[i][ :2]
    plot_all_auxs(int(x), int(y), ghosts, intersects, ficts, auxs)
  #    plt.ylim(zmax, zmin)
  #    plt.gca().set_aspect('equal', adjustable='box')
    #plt.savefig(str(x).rjust(7, '0') + '-' + str(y) + '.png', format = 'png')
    #plt.close()
    #plt.show()
  
  for j, x in enumerate(X):
    i_list = pick_ghost_by_x(ghosts, x)

    for i in i_list:
      
      fig, ax = plt.subplots(1, 1, sharex = True)
      fig.set_size_inches(longer_side, square_side, forward = True) # IN 
      #plt.subplot(2, 1, 1)
      #plot_ghost_interp(-1, ghosts[i], intersects[i], ficts[i], auxs[i])
      #plt.subplot(2, 1, 2)
      Plot_FS(fs_file)
      plot_all_ghosts(ghost_coords_file)
      plot_auxs_2D(ghosts[i], intersects[i], ficts[i], auxs[i])
      xyz = str(ghosts[i][0]) + '-' + str(ghosts[i][1]) + '-' + str(ghosts[i][2])
      # STATIONARY WINDOW - BETTER ANIMATION
      #Zoom_To_Point_2D(x_mid, z_mid, len(X) * 2, True)
      # MOVING WINDOW - BETTER ZOOM => QC
      #Zoom_To_Point_2D(ghosts[i][0], ghosts[i][2], 8, True)
      plt.xlim(ghosts[i][0] - 5, ghosts[i][0] + 5) #FIXME ADAPT TO STORDER
      plt.ylim(ghosts[i][2] - 2, ghosts[i][2] + 7)
      plt.gca().invert_yaxis()
      plt.gca().set_aspect('equal', adjustable = 'box')
      plt.savefig(pattern + '_' + xyz + '.png', format = 'png', dpi = fig_resol, bbox_inches = 'tight')
      plt.close()

  print('lib_fwi_fs/show_middle_ghosts: END')  
  

  
  


  
  
# PLOTT ALL AUXILIARY POINTS (I, F, A) ASSOCIATED WITH A SINGLE GHOST NODE
def plot_all_auxs(x, y, z, ghosts, intersects, ficts, auxs):
  print('lib_fwi_fs/plot_all_auxs: START')  
  
  from lib_generic import cherrypick_x, cherrypick_xy
  if False: #FIXME #len(Y_fs) > 1: # 3D (HANDLE TICKS)
    fig, ax = Plot_FS_3D(X_fs, Y_fs, Z_fs)
    ghost, intersect, ficts, auxs = cherrypick_xy(x, y, ghosts, intersects, ficts, auxs)
    plot_auxs_3D(ax, ghost, intersect, ficts, auxs)
  
  else: # 2D
    
    n = 0
    index = 'err'
    for g in ghosts:
      if (g[0] == x) and (g[1] == y) and (g[2] == z):
        index = n
        break
      n += 1

    #ghost, intersect, ficts, auxs = cherrypick_x(x, ghosts, intersects, ficts, auxs)
    ghost = ghosts[index]
    intersect = intersects[index]
    ficts = ficts[index]
    auxs = auxs[index]
    ghost = ghost[ :3]
    plot_auxs_2D(ghost, intersect, ficts, auxs) # [0] <= COMPATIBILITY W/ 3D  
  
  print('lib_fwi_fs/plot_all_auxs: END')
  
# PLOTT ALL 3D AUXILIARY POINTS (I, F, A) ASSOCIATED WITH A SINGLE GHOST NODE
def plot_auxs_3D(ax, ghost, intersect, ficts, auxs):
  colors = iter(cm.rainbow(np.linspace(0, 1, len(ficts))))
  ax.scatter(ghost[0], ghost[1], ghost[2], c = 'b')
  ax.scatter(intersect[0], intersect[1], intersect[2], c = 'm')
  
  X, Y, Z = Split_Tuples(ficts)
  ax.scatter(X, Y, Z, c = 'k', s = psize)
  # PLOTT NORMAL
  ax.plot([ghost[0], X[-1]], [ghost[1], Y[-1]], [ghost[2], Z[-1]], c = 'grey', lw = 3, alpha = 0.5)
  ##
  
  for a in auxs:
    X, Y, Z = Split_Tuples(a)
    ax.scatter(X, Y, Z, c = next(colors))    

# PLOTT ALL 2D AUXILIARY POINTS (I, F, A) ASSOCIATED WITH A SINGLE GHOST NODE
def plot_auxs_2D(ghost, intersect, ficts, auxs):
  from lib_generic import Split_Tuples
  plt.scatter(ghost[0], ghost[2], c = 'k', s = psize)
  plt.scatter(intersect[0], intersect[2], c = 'w', s = psize)
  # PLOTT FICTS
  X, Y, Z = Split_Tuples(ficts)
  plt.scatter(X, Z, c = 'k', s = psize)
  # PLOTT NORMAL
  plt.plot([ghost[0], X[-1]], [ghost[2], Z[-1]], c = 'k')#, lw = nwidth)#, alpha = nopac)
  ##
  # PLOTT AUXS
  #colors = iter(cm.rainbow(np.linspace(0, 1, len(ficts)))) # IT SEEMS AUXS ARE IN WRONG ORDER 
  zmax = epsi
  for a in auxs:
    X, Y, Z = Split_Tuples(a)
    if max(Z) > zmax:
      zmax = max(Z)
    plt.scatter(X, Z, c = 'b', s = psize) #next(colors))  
  #plt.ylim(zmax + pad, ghost[2] - pad)
  




#
def plot_ghost_interp(wf_file, ghost, isect, ficts, auxs):
  print('lib_fwi_fs/plot_ghost_interp: START')
  
  from lib_fwi_generic import read_wavefield
  import numpy as np
  
  mode = 'auxs'
  
  if wf_file != -1:
    X, Z, wf = read_wavefield(wf_file)
  
  if mode == 'geometry':
    xf, zf = [], []
    for i in range(len(ficts)):
      xf.append(ficts[i][0])
      zf.append(ficts[i][2])
    
    xx = [ghost[0], isect[0], rpoint[0]]
    zz = [ghost[2], isect[2], rpoint[2]]
    
    plt.scatter(xf, zf)
    plt.scatter(xx, zz)
    
    #plt.gca().invert_yaxis()
  
  elif mode == 'auxs':
    from lib_fwi_fs import prepare_auxs_interp_geometry
    
    xx = prepare_auxs_interp_geometry(ficts, auxs)
    
    wf_a = list(np.ones(len(auxs[0]))) 
    #wf_a_random = [1, 10, 2, 1] 
     
    plt.subplots(len(ficts), 1, sharex = True)
    #fig.set_size_inches(longer_side, square_side, forward = True)
    
    zz = []
    for i in range(len(ficts)):
      wf_f = 0.
      if wf_file != -1:
        wf_a = []
      for j in range(len(auxs[i])):
        x_aux = int(auxs[i][j][0])
        y_aux = int(auxs[i][j][1])
        z_aux = int(auxs[i][j][2])
        
        if wf_file != -1:
          wf_aux = wf[x_aux][z_aux] #FIXME 2D AS YET
          wf_a.append(wf_aux)
        
        c_aux = auxs[i][j][3]
        wf_f += c_aux * wf_a[j]
      f_zz = wf_a #+ [wf_f]
      
      plt.subplot(len(ficts), 1, i + 1)
      plt.scatter(xx[i][:-1], f_zz, c = 'b')
      plt.scatter(0., wf_f, c = 'r')
      plt.xlim(-5, 5)
      #plt.ylim(100, 750)
      
    #plt.scatter(xx, zz)
      
  elif mode == 'ficts':
    from lib_fwi_fs import prepare_fict_interp_geometry
    
    
    xx = prepare_fict_interp_geometry(ghost, isect, ficts)
    
    #wf_f = list(np.zeros(len(ficts)))
    wf_f = [1, 2, -10, 4]
    
    wf_rpoint = 0.
    for i in range(len(ficts)):
      c_fict = ficts[i][3]
      wf_rpoint += c_fict * wf_f[i]
    
    zz = [-wf_rpoint, 0., wf_rpoint] + [wf_f]
    
    plt.scatter(xx, zz)
    
  else:
    eprint('liplot_fwi_fs/plot_ghost_interp: ERROR! UNKNOWN MODE!')
    quit()
  
  print('lib_fwi_fs/plot_ghost_interp: END')  



#
def plot_inext_nodes(inext, fs_file, ghost_file):
  #FIXME: READ IT FROM params_mod.f90
  in_flag = 0
  acc_flag = -1
  ext_flag = -666
  ##
  
  x = np.array([float(i[0]) for i in inext])
  z = np.array([float(i[2]) for i in inext])
  val = np.array([float(i[3]) for i in inext])
  
  inn, acc, ext = [], [], []
  
  for i in range(len(x)):
    if val[i] == in_flag:
      inn.append([x[i], z[i]])
    
    elif val[i] == acc_flag:
      acc.append([x[i], z[i]])
  
    elif val[i] == ext_flag:
      ext.append([x[i], z[i]])                      
  
    else:
      eprint('my_plotlib/plot_inext_nodes: Error! Unknown flag.\n')
      quit()
  
  grid = val.reshape((int(z.max()), int(x.max()))) # YES, THIS ORDER
  
  
  #fig, ax = Plot_FS(fs_file)
  
  plt.xlim(min(x), max(x))
  plt.ylim(max(z), min(z)) 
  
  #plot_all_ghosts(ghost_file)








#
# PRESUMABLY OBSOLETE
#


# SHOW GEOMETRY BEFORE AND AFTER ROTATION IN ONE PLOT
def plot_flat_n_rotated(project, rot_rec, rot_fs, rot_sou, rec, fs, sou, nz):
  fig, ax = plt.subplots(1, 1, sharex = True)
  fig.set_size_inches(longer_side, square_side, forward = True) # IN INCHES
    
  rot_rec_x, rot_rec_z = [r[0] for r in rot_rec], [r[1] for r in rot_rec]
  rot_fs_x, rot_fs_z = [r[0] for r in rot_fs], [r[1] for r in rot_fs]
  
  X, Y, Z = fs
  Z = [z[0] for z in Z]
  
  plt.plot(X, Z, c = 'grey')
  plt.plot(rot_fs_x, rot_fs_z, c = 'grey')
  plt.scatter([r[0] for r in rec], [r[1] for r in rec])
  plt.scatter(rot_rec_x, rot_rec_z, c = 'y')
  plt.scatter(sou[0], sou[1], s = 200, c = 'b', marker = '*')
  plt.scatter(rot_sou[0], rot_sou[1], s = 200, c = 'y', marker = '*')
  
  #pad = 40
  #plt.xlim([-pad, nx1 + pad])
  plt.ylim([-0.5, nz])
  #plt.xticks(np.arange(nx1))
  #plt.yticks(np.arange(nx3))
  #Display_Every_Nth_Tick(ax, 20)
  plt.gca().invert_yaxis()
  plt.gca().set_aspect('equal', adjustable='box')
  plt.grid()
  plt.savefig(project + '_rotated.png', dpi = 150)





def OLD_Plot_FS(Z, **kwargs):
  """
  Framework function to plot a gridded free surface. 
  
  Parameters
  ----------
  Z : array / list
    Gridded 3D data structure of (nx1*nx2*1) shape. 
  **kwargs : keyword arguments, optional
    Just passing it down.
      
  Returns
  -------
  0
  
  Notes
  -----
  Works both for 2D/3D.
  It has similar role and structure to Plot_Model.
  
  """    
  this_func = this_lib + 'Plot_FS: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  # COLLAPSE THE ARRAY 3D->2D (Z-TRACES OF FS HAVE ONLY 1 VALUE EACH)
  Z = [[j[0] for j in i] for i in Z]
  Z = np.array(Z)
  
  nx1 = Z.shape[0]
  nx2 = Z.shape[1]
   
  # 3D 
  if (nx1 > 1) and (nx2 > 1):
    Plot_FS_3D(Z, **kwargs)
  
  # 2D #FIXME: DO WE NEED IT?
  elif ((nx2 == 1) and (nx1 > 1)):
    x = list(range(1, nx1 + 1))
    y = [i[0] for i in Z]
    ax = plt.gca()
    plt.plot(x, y)
    plt.xlabel('x')
    plt.ylabel('z')    
    
  elif ((nx1 == 1) and (nx2 > 1)): 
    x = list(range(1, nx2 + 1))
    y = Z[0]
    plt.plot(x, y)
    plt.xlabel('y')
    plt.ylabel('z')
    
  else:
    eprint(this_func + 'Error! Wrong geometry.\n')
    eprint(this_func + 'nx1: ' + str(nx1) + '\n') 
    eprint(this_func + 'nx2: ' + str(nx2) + '\n')
    quit()
    
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def OLD_Plot_FS_3D(Z, **kwargs):
  """
  Plot a 3D free surface as a 3D surface plot / 2D projection (map) / 2D slice / ...
  
  Parameters
  ----------
  Z : array / list
    Gridded 3D data structure of (nx1*nx2) shape. 
  
  **kwargs : keyword arguments (optional)
      Current capabilities:
      plot_type : str
        At the moment possible types are: 'map', 'surf', 'slice'.
        Default: 'map'.
        
  Returns
  -------
  0
  
  Notes
  -----
  As the output of: 
  
  >>> X, Y = meshgrid(x, y)
  
  follows X.shape = Y.shape = (len(y), len(x)), the shape of Z
  must also be (len(y), len(x)), i.e. as if x was faster index:
  
  >>> for y in range(len(y)):
  >>>   for x in range(len(x)):
  >>>     Z[y][x] = ...
  
  This requires transposition with respect to the .vtr convention.
  
  """     
  this_func = this_lib + 'Plot_FS_3D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  # CHECK KEYWORD ARGUMENTSs 
  try:
    plot_type = kwargs['plot_type']
  except KeyError:
    plot_type = 'map'
    eprint(this_func + "Warning. plot_type not specified - assuming 'map'.\n")

  # PLOT
  if plot_type == 'surf':
    from lib_generic_PLOTT import Plot_Surface
    Plot_Surface(Z, **kwargs)
  
  elif plot_type == 'map':
    from lib_generic_PLOTT import Plot_Map
    Plot_Map(Z, **kwargs)
    
  elif plot_type == 'slice': 
    from lib_generic_PLOTT import Plot_Slice
    Plot_Slice(Z, **kwargs)
    
  elif plot_type == 'cube': # FIXME: NOT YET ADAPTED TO FS
    from lib_generic_PLOTT import Plot_Cube
    Plot_Cube(Z, **kwargs)
    
  else:
    raise ValueError('Unknown plot_type: ' + plot_type + '\n')
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------




