"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permision writing to k.chrapkiewicz17@imperial.ac.uk.

This library provides high-level plotting functions 
for FWI I/O objects such as models, wavefields, etc.

"""

## MODULES
import numpy as np
import matplotlib.pyplot as plt

from lib_generic_CONST import *#verbos_func
from lib_fwi_project_CONST import *
from lib_generic import *
#from lib_fwi_generic import Read_vtr, Save_vtr
from lib_io_fullwave import *
from lib_generic_PLOTT import Plot, Plot_File, Plot_Data

## CONSTS
this_lib = 'lib_fwi_generic_PLOTT.py/'

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
# JUXTAPOSITIONS
# -------------------------------------------------------------------------------


def Plot_Juxtaposition(Z_list, **kwargs):
  """
  FIXME: NEEDS TO BE MERGED WITH Plot_2D
  IN SOME CLEVER WAY.
  
  Parameters
  ----------
  Z_list : list of arrays
    Currently max 2 arrays.
  
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  It can be used for models and other arrays too.
  NOTE: Therefore merge with Plot_2D would be great. 
  
  """
  
  this_func = this_lib + 'Plot_Juxtaposition: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  from lib_math_signal import Data_Modify
  
  if len(Z_list) > 2:
    raise NotImplementedError('Cant juxtapose more than 2 files\n')
  
  juxt = Kwarg('juxt', 'alternate_cmaps', kwargs)
  cmap = Kwarg('cmap', 'seismic', kwargs)
  cbar = Kwarg('cbar', 0, kwargs)
  kwargs['alpha'] = Kwarg('alpha', 1, kwargs)
  
  if juxt == 'alternate_cmaps':
    kwargs['cmap'] = cmap
    # N SYNTH, N OBS, N SYNTH, ...
    Plot_Data(Z_list[0], interleave=Z_list[1], plt_func='imsh', yflip=1, **kwargs) #, fig_size=[15,8]
  
  elif juxt == 'alternate__wiggles':
    # N SYNTH, N OBS, N SYNTH, ...
    pass

  elif juxt == 'wigg_beside_wigg':
    # FOR EACH CHANNEL 2 TRACES (SYNTH AND OBS) CLOSE TO EACH OTHER
    pass 

  elif juxt == 'wigg_on_wigg':
    # RED/BLUE INTERTWINED TRACES (SAME AS WAVELET PLOTS)
    pass 
  
  elif juxt == 'wigg_on_cmap':
    from lib_generic_PLOTT import Plot_Wiggles
    kwargs['cmap'] = cmap
    kwargs['cbar'] = cbar #FIXME: DOESN'T SHOW UP
    kwargs['alpha'] = 1
    Plot_Data(Z_list[1], plt_func='imsh', **kwargs) # 0.2
    kwargs['yflip'] = True #NOTE: AD-HOC
    Plot_Wiggles(0.5*(Z_list[0]/np.max(Z_list[0])), gap=1, **kwargs)
    
  
  elif juxt == 'cmap_on_cmap': 
    # FIXME: IT PLOTS 2 FIGURES
    kwargs['cmap'] = 'Greys'
    kwargs['cbar'] = 0
    Plot_Data(Z_list[0], plt_func='imsh', alpha=1, **kwargs) #, fig_size=[15,8]
    kwargs['cmap'] = cmap
    kwargs['cbar'] = cbar
    Plot_Data(Z_list[1], plt_func='imsh', alpha=1, **kwargs) # 0.2
  
  elif juxt == 'diff':
    #NOTE: BEST FOR MODELS (EVEN THOUGH IN FULLWAVE YOU CAN PLOT & DUMP UPDATES)
    kwargs['cmap'] = cmap  
    g1 = Z_list[0]
    g2 = Z_list[1]
    
    # NORMALIZE
    kwargs['normalize'] = 'rms'
    g1 = Data_Modify(g1, **kwargs) 
    g2 = Data_Modify(g2, **kwargs)     
    # EQUALIZE Z-DIM # NOTE: WE REQUIRE/ASSUME THIS
    #if len(g1[0][0]) < len(g2[0][0]):
    #  g2 = np.array([[i[:len(g1[0][0])] for i in j] for j in g2])
    #elif len(g1[0][0]) > len(g2[0][0]):
    #  g1 = np.array([[i[:len(g2[0][0])] for i in j] for j in g1])
    Plot_Data(g1-g2, plt_func='imsh', yflip=1, **kwargs)
      
  else: 
    raise ValueError('Unknown value of juxt: ' + juxt)

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
# COMPOSITE PLOT INCL. WAVEFIELD
# -------------------------------------------------------------------------------


def Plot_All(proj, **kwargs):
  """
  Plot model + free surface + 
  + sources and receivers + wavefield.
  
  Parameters
  ----------
  sr : tuple 
    (x, y, z) of a single source / receiver
  
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  NOTE wavefield requires extended (etop etc.) domain.
  
  Table showing what to extend for each configuration
  wavefield and model types ('non-model' stands for 
  Gradient, Precon etc. which are already extended)
  
         |  model | non-model |  
  ----------------------------|
  fw on  | mod,sr |   sr      |
  -------|--------------------|
  fw off |  none  |   sr      | 
  -----------------------------
  
  """
  
  this_func = this_lib + 'Plot_All: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_io_fullwave import Read_vtr
  from lib_generic_PLOTT import Plot_Slices_XYZ, Plot_Slice
  from lib_fwi_generic import Set_Extra_Nodes 
  
  kwargs['cmap'] = Kwarg('cmap', 'inferno', kwargs)
  units = Kwarg('units', 'nodes', kwargs)
  slice_coord = Kwarg('slice_coord', None, kwargs) # THIS ALLOWS TO PLOT A SINGLE SLICE
  kwargs['xyz'] = Kwarg('xyz', np.array([0,0,0]), kwargs)
  mod_type = Kwarg('mod_type', 'model', kwargs)
  sr = Kwarg('sr', None, kwargs)
  fw = Kwarg('fw', None, kwargs)
  mod = Kwarg('mod', None, kwargs)
  fs = Kwarg('fs', None, kwargs)
  fname = proj.inp.path + proj.name + '-Runfile.key'
  etop, eleft, efront, ebot, eright, eback = Set_Extra_Nodes(fname, **kwargs)  
  
  # IMPLEMENTATION OF THE TABLE (SEE Notes ABOVE)  
  extend_sr = True
  extend_mod = False
  
  if fw and (mod_type == 'model'):
    extend_mod = True
    
  if extend_mod or (mod_type != 'model'):
    kwargs['xyz'] += np.array([eleft,
                               efront,
                               etop])
  
  if not fw and (mod_type == 'model'):
    extend_sr = False
  
  vols = []
  surfs = []
  scatts = []
  
  
  # -------------------------------------------------------------------------
  # READ SOURCE/RECEIVER 
  # -------------------------------------------------------------------------
  if sr:
    if extend_sr:
      sr[0] += eleft
      sr[1] += efront
      sr[2] += etop
    scatt = [sr] #NOTE: NOT-PROJECTED COORDS WILL BE ACCURATE REAL NUMBERS
    scatts.append(scatt)
    sr_int = [int(i) for i in sr] # INT FOR SLICING
    kwargs['xyz'] = sr_int
  else:
    scatt = None  
  
  
  # -------------------------------------------------------------------------
  # READ MODEL (IT CAN ALSO BE A GRADIENT ETC.)
  # -------------------------------------------------------------------------  
  if type(mod) == type(np.array([])) or isinstance(mod, str):
    #print this_func, 'mod', mod 
    #if 'RawGrad' in mod or 'RawPrec' in mod:
    mod = Read_Array(mod, **kwargs)
    
    #print this_func, 'mod, min, max', np.min(mod), np.max(mod)
    
    if extend_mod:
      mod_ext = np.ones((mod.shape[0] + eleft + eright, 
                         mod.shape[1] + efront + eback, 
                         mod.shape[2] + etop + ebot)) 
      # NOTE: COLOR RANGE PRESERVED AND MARGINS STILL 
      # DISTINCTIVE (UNLESS HORIZONTAL SLICE THROUGH WATER)    
      mod_ext *= np.min(mod)
      
      
      # NOTE: WE CAN'T SIMPLY USE eleft:-eright IF eright=0, ETC.
      mod_ext[eleft:mod_ext.shape[0]-eright, 
              efront:mod_ext.shape[1]-eback, 
              etop:mod_ext.shape[2]-ebot] = mod
      mod = np.array(mod_ext)
    vols.append(mod)
  else:
    if verbos > 0:
      print((this_func + 'mod=None => cmap=seismic\n'))
    kwargs['cmap'] = 'seismic'
  
  # -------------------------------------------------------------------------
  # READ FS
  # -------------------------------------------------------------------------  
  if fs:
    fs = Read_Array(fs, **kwargs)
    if units == 'm':
      fs = fs * proj.inp.fs.dxs[2] # CONVERT TO METRES
      surfs.append(fs)
    
    
  # -------------------------------------------------------------------------
  # READ WAVEFIELD (FORWARD/BACKWARD/...)
  # -------------------------------------------------------------------------   
  if fw:
    fw = Read_Array(fw, **kwargs)
    # NOTE: WAVEFIELD IS ALREADY EXTENDED!
    vols.append(fw)
  
  
  # -------------------------------------------------------------------------
  # PLOT EVERYTHING SLICED
  # -------------------------------------------------------------------------   
  if units == 'm':
    kwargs['box'] = proj.box # NEEDED FOR CONVERSION TO METRES VIA extent   
  
  plot_args = dict(kwargs)
  plot_args['vols'] = vols
  plot_args['surfs'] = surfs  
  plot_args['scatts'] = scatts
  plot_args['marker'] = Kwarg('marker', "*", kwargs)
  plot_args['marker_size'] = Kwarg('marker_size', 250, kwargs)
  plot_args['marker_color'] = Kwarg('marker_color', "white", kwargs)
  plot_args['marker_edge'] = Kwarg('marker_edge', "k", kwargs)  
  
  if slice_coord: # SINGLE SLICE
    if 'xyz' in plot_args: 
      if slice_coord == 'x':
        plot_args['coord_value'] = plot_args['xyz'][0]
      elif slice_coord == 'y':
        plot_args['coord_value'] = plot_args['xyz'][1]
      elif slice_coord == 'y':
        plot_args['coord_value'] = plot_args['xyz'][2]
    del plot_args['xyz']
  
  Plot_Slice(**plot_args) # THIS IS A FRAMEWORK FOR SINGLE/XYZ SLICES
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
# PLOT MODELS
# -------------------------------------------------------------------------------


def Plot_Model(Z, **kwargs):
  """
  Plot a 3D model.
  
  Parameters
  ----------
  Z : array / list
    Gridded 3D data structure of (nx1*nx2) shape. 
  **kwargs : keyword arguments, optional
      Current capabilities:
      plot_type : str
        Default: 'slice'.
      
  Returns
  -------
  0
  
  Notes
  -----
  
  """     
  this_func = this_lib + 'Plot_Model: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  # CHECK KEYWORD ARGUMENTS
  
  cmap = Kwarg('cmap', 'inferno', kwargs)
  kwargs['cmap'] = cmap
  #kwargs['yflip']
  
  try:
    plot_type = kwargs['plot_type']
  except KeyError:
    plot_type = 'slice'
    if verbos > 1:
      eprint(this_func + "Warning. plot_type not specified - assuming 'slice'.\n")
  
  # PLOT
  if plot_type == 'slice': 
    from lib_generic_PLOTT import Plot_Slice
    Plot_Slice(vols=[Z], **kwargs)
  
  # THIS IS CHOSEN IN Plot_Slice
  elif plot_type == 'slices_xyz':
    from lib_generic_PLOTT import Plot_Slices_XYZ
    Plot_Slices_XYZ(vols=[Z], **kwargs)
  
  elif plot_type == 'slice_anim':
    from lib_generic_PLOTT import Plot_Slice_Animation
    Plot_Slice_Animation(vols=[Z], **kwargs)
  
  elif plot_type == 'cube':
    from lib_generic_PLOTT import Plot_Cube
    Plot_Cube(Z, **kwargs)
  
  else:
    raise ValueError('Unknown plot_type: ' + plot_type)
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
# PLOT FREE SURFACE
# -------------------------------------------------------------------------------


def FS_Plot_Full(proj_name, **kwargs): #FIXME D-C?
  """
  Plot project's sources and 
  receivers.
  
  Parameters
  ----------
  proj_name : str
    Project name assumed to be provided 
    as a first argument of the FWI code
    being a prefix of all FWI input files.  
  **kwargs : keyword arguments, optional
      Just passing it down.

  Returns
  -------
  0
  
  Notes
  -----
  
  """    
  this_func = this_lib + 'FS_Plot_Full: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  proj_path = Kwarg('proj_path', './', kwargs)
  adjust = Kwarg('adjust', False, kwargs)
  #plot_type = Kwarg('plot_type', 'map', kwargs)
  fname = Kwarg('fname', proj_name+'-FreeSurf.vtr', kwargs)
  
  nx1_fs, nx2_fs, nx3_fs, fs_z = Read_vtr(fname, **kwargs)
  xcoord = Kwarg('xcoord', nx1_fs/2, kwargs)
  ycoord = Kwarg('xcoord', nx2_fs/2, kwargs)  
  
  fig, axes = plt.subplots(nrows=2, ncols=2, figsize=[15,8])
  
  ax = fig.add_subplot(2, 2, 1, projection='3d')
  Plot_Data(fs_z, data_type='FS', plot_type='surf', ax=ax)
  
  plt.subplot(2,2,2)
  Plot_Data(fs_z, data_type='FS', plot_type='map')
  #plt.gca().remove()
  
  plt.subplot(2,2,3)
  Plot_Data(fs_z, data_type='FS', plot_type='slice', slice_coord='y', coord_value=ycoord, cbar=0)
  
  plt.subplot(2,2,4)
  Plot_Data(fs_z, data_type='FS', plot_type='slice', slice_coord='x', coord_value=xcoord, cbar=0)
    
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_FS_Roughness(fs_z, **kwargs): #FIXME D-C?
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
  this_func = this_func = this_lib + 'Plot_FS_Roughness: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  from lib_generic_PLOTT import Plot_Map
  
  curv = Kwarg('curv', 'hess', kwargs)
  
  fs_z = np.array([[i[0] for i in j] for j in fs_z])
  
  if curv == 'hess': # FIXME: IT'S A GRADIENT NOT HESSIAN, BTW CURVATURE HAS MORE COMPLEX DEF.
    # HESSIAN
    grad = np.gradient(fs_z)
    shp = np.shape(np.array(grad))
    abses = np.zeros((shp[1], shp[2]))
    for x in range(shp[1]):
      for y in range(shp[2]):
        abs_grad = np.sqrt(grad[0][x][y]**2 + grad[1][x][y]**2)
        abses[x][y] = abs_grad
    
    Plot_Map(abses, cmap='viridis', **kwargs)
  
  elif curv == 'neigh':
    # MY MEASURE OF ROUGHNESS 
    from lib_fwi_fs import Nearest_Neighbours
    
    fs_shp = np.shape(fs_z) 
    print(fs_shp)
    max_diffs = np.zeros((fs_shp[0], fs_shp[1]))
    for x in range(fs_shp[0]):
      for y in range(fs_shp[1]):
        max_diff = 0
        for dx in [-1, 1]:
          for dy in [-1, 1]:
            try:
              diff = abs(fs_z[x][y] - fs_z[x+dx][y+dy])
            except IndexError:
              continue 
            if diff > max_diff:
              max_diff = diff
        max_diffs[x][y] = max_diff  
    
    plt.figure()
    Plot_Map(max_diffs, cmap='viridis', **kwargs)
  
  #print this_func + 'END'
  return 0  


# -------------------------------------------------------------------------------


def Plot_FS(Z, **kwargs):
  """
  Plot a 3D free surface
  
  Parameters
  ----------
  Z : array / list
    Gridded 3D data structure of (nx1*nx2) shape. 
  **kwargs : keyword arguments, optional
      Current capabilities:
      plot_type : str
        Default: 'slice'.
      
  Returns
  -------
  0
  
  Notes
  -----
  
  """     
  this_func = this_lib + 'Plot_FS: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  # CHECK KEYWORD ARGUMENTS
  
  cmap = Kwarg('cmap', 'viridis_r', kwargs)
  kwargs['cmap'] = cmap
  #kwargs['yflip']
  
  try:
    plot_type = kwargs['plot_type']
  except KeyError:
    plot_type = 'slice'
    if verbos > 0:
      eprint(this_func + "Warning. plot_type not specified - assuming 'slice'.\n")
  
  # PLOT
  if plot_type == 'slice': 
    from lib_generic_PLOTT import Plot_Slice
    Plot_Slice(surfs=[Z], **kwargs)
  
  elif plot_type == 'slices_xyz':
    from lib_generic_PLOTT import Plot_Slices_XYZ
    Plot_Slices_XYZ(surfs=[Z], **kwargs) # NOTE: surfs
  
  elif plot_type == 'map':
    from lib_generic_PLOTT import Plot_Map
    Z = np.array([[j[0] for j in i] for i in Z])
    Plot_Map(Z, **kwargs)

  elif plot_type == 'surf':
    from lib_generic_PLOTT import Plot_Surface
    Plot_Surface(Z, **kwargs)  
  
  else:
    raise ValueError('Unknown plot_type: ' + plot_type)
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
# PLOTT SOURCES & RECEIVERS
# -------------------------------------------------------------------------------


def Plot_SR(Ss, Rs, **kwargs):
  """
  Plot project's sources and 
  receivers.
  
  Parameters
  ----------
  proj_name : str
    Project name assumed to be provided 
    as a first argument of the FWI code
    being a prefix of all FWI input files.  
  **kwargs : keyword arguments, optional
      Just passing it down.

  Returns
  -------
  0
  
  Notes
  -----
  
  """    
  this_func = this_lib + 'Plot_SR: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_fwi_generic import SR_Read
  from lib_generic_PLOTT import Plot_Slices_XYZ
  
  proj_path = Kwarg('proj_path', './', kwargs)
  adjust = Kwarg('adjust', False, kwargs)
  plot_type = Kwarg('plot_type', 'proj', kwargs)
  separate = Kwarg('separate', False, kwargs)
  #reciprocity = Kwarg('reciprocity', False, kwargs)
  pairs = Kwarg('pairs', False, kwargs)
  #annot_R = Kwarg('annot_R', False, kwargs)
  
  try:
    fig_size = kwargs['fig_size']
    plt.figure(figsize=fig_size)  
  except KeyError:
    pass  
  
  c_R = Kwarg('c_R', 'r', kwargs)
  c_S = Kwarg('c_S', 'grey', kwargs)
  s_R = Kwarg('s_R', 50, kwargs)
  s_S = Kwarg('s_S', .5, kwargs)
  
  ##print this_func, 'adjust', adjust
  #dims, Ss, Rs = SR_Read(proj_name, adjust, **kwargs)
  ##if reciprocity: # NOTE: NOT HERE?
  #  #Ss, Rs = Rs, Ss
  
  # CONVERT TO LIST # FIXME: NECESSARY?
  Ss = [Ss[key] for key in Ss]
  Rs = [Rs[key] for key in Rs]
  
  #if pairs and not separate:
    #xy_list = Plot_SR_Pairs(Ss, Rs, **kwargs)  
    #kwargs['curves'] = xy_list
  
  Plot_Slices_XYZ(scatts=[Ss, Rs], **kwargs)
  #Plot_Points(Ss, plot_type=plot_type, c=c_S, s=s_S, label='input Ss',  **kwargs)  
  #if separate:
  #  plt.figure()
  #Plot_Points(Rs, plot_type=plot_type, c=c_R, s=s_R, label='input Rs', **kwargs) #alpha=.5)  
  #plt.legend()
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_SR_Pairs(Ss, Rs, **kwargs):
  """
  Plot line for each source-receiver pair.
  
  Parameters
  ----------
  Ss : list 
    List of sources XY.
  
  Rs : list 
    List of receivers XY.  
  
  **kwargs : keyword arguments, optional
      Just passing it down.

  Returns
  -------
  0
  
  Notes
  -----
  It is much more meaningful for trave-times.
  
  """    
  this_func = this_lib + 'Plot_SR_Pairs: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  max_offset = Kwarg('max_offset', big, kwargs)
  #try:
    #dx = kwargs['dx']
  #except KeyError:
    #eprint(this_func + 'Error. You need to specify dx for info.\n')
    #quit()
  
  
  plt.title('Interstation lines only for offset < ' + str(max_offset) + ' km.')
  
  xy_list = []
  i = 0
  for R in Rs:
    for S in Ss:
      #print 'R, S', R, S
      rx, ry = R[0], R[1]
      sx, sy = S[0], S[1]
      
      offset = np.sqrt((rx-sx)**2 + (ry-sy)**2) #* dx / 1000. # KM
      
      if offset <= max_offset:      
        x = [rx, sx]
        y = [ry, sy]
        xy_list.append([x, y])
        #plt.plot(x, y, '-', c='grey', alpha=.1)
        
      i += 1
  
  #print i
    
  #print this_func + 'END'
  return xy_list


# -------------------------------------------------------------------------------
# PLOTT SEISMIC TRACES (WAVELET & DATA)
# -------------------------------------------------------------------------------


def Plot_Wavelet_Full(fname_or_Z, dt, **kwargs):
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
  
  this_func = this_lib + 'Plot_Wavelet_Full: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  proj_path = Kwarg('proj_path', './', kwargs)
  xlim_time = Kwarg('xlim_time', None, kwargs)
  xlim_freq = Kwarg('xlim_freq', None, kwargs)
  figsize = Kwarg('figsize', [11, 5], kwargs) # DELIBERATELY NOT figsize_default
  title = Kwarg('title', None, kwargs)
  f_max = Kwarg('f_max', f_max_default, kwargs)
  xlim_freq = [0, 1.5 * f_max]
  
  #fname = proj_path + proj_name + '-RawSign.sgy'
  
  fig, axes = plt.subplots(nrows=1, ncols=2, figsize=figsize)
  if title:
    fig.suptitle(title)
    del kwargs['title']
  
  plt.subplot(1,2,1)
  Plot(fname_or_Z, data_type='wavelet', dt=dt, xlim=xlim_time, **kwargs)
  plt.subplot(1,2,2)
  Plot(fname_or_Z, data_type='wavelet', spectrum='ampl', dt=dt, xlim=xlim_freq, **kwargs)
  fig.tight_layout()

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Trace_Evolution(trace, **kwargs):
  """
  Plot a trace in a wiggly style of suxwigb
  as evolving during source-excitation.
  
  Parameters
  ----------
  trace : list / array 
    Time series A(t).
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  
  
  """
  
  this_func = this_lib + 'Plot_Trace_Evolution: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  snap_step = Kwarg('snap_step', 100, kwargs)
  # NOTE: THE PATH IS './' UNLESS IT'S A PART OF fname_core
  fname_core = Kwarg('fname_core', 'trace_evolution', kwargs) 
  
  nsamp = len(trace)
  
  for i in range(nsamp):
    if i % snap_step != 0:
      continue
    
    plt.figure()
    Plot_Trace(trace, mark_sample=i, **kwargs)
    snap_no = str(i).rjust(7, '0')
    plt.savefig(fname_core + '_' + snap_no + '.png')
    plt.close()
    
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Trace(trace, **kwargs):
  """
  Plot a trace in a wiggly style of suxwigb.
  
  Parameters
  ----------
  trace : list / array 
    Time series A(t).
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  Taking the spectrum should be separated from Plot_Trace 
  or called from external.  
  
  """
  
  this_func = this_lib + 'Plot_Trace: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_generic_PLOTT import Plot_2_Series
  from lib_math_num import Derivative
  
  ylabel = Kwarg('ylabel', 'amplitude [-]', kwargs)
  spectrum = Kwarg('spectrum', None, kwargs)
  mark_sample = Kwarg('mark_sample', -1, kwargs)
  
  deriv = Kwarg('deriv', 0, kwargs) # DERIVATIVE (0 - NOTHING, 1 - FIRST, ...)
  scale = Kwarg('scale', 1, kwargs)
  
  try:
    dt = kwargs['dt']
    xlabel = 'time [s]'
  except KeyError:
    eprint(this_func + 'Set default dt=1 => X-axis scaling may be wrong.\n')
    dt = 1
    xlabel = 'sample'
  
  # DIFFERENTIATE
  if deriv == 0:
    pass 
  elif deriv == 1:
    trace = Derivative(trace, dx=dt)
  else:
    raise ValueError('This derivative is not implemented yet.')
  
  # SCALE
  if scale != 1:
    trace = np.array(trace) * scale
  
  N = len(trace)
  T = dt*N
  
  df = 1. / T
  
  if verbos > 2:
    print(this_func, 'T=', T)
    print(this_func, 'df=', df)
    
  abcissas = np.array([df*n if n<N/2 else df*(n-N) for n in range(N)])
  # OR:
  #freqs = np.fft.fftfreq(data.size, time_step)
  # ABOVE BUT BUILT-IN
  #idx = np.argsort(freqs)
  #plt.plot(freqs[idx], ps[idx])
  
  if spectrum == 'ampl':
    trace = abs(np.fft.fft(trace))
    #print len(abcissas), len(trace)
    #abcissas = np.pad(abcissas, [1,1], 'constant')
    #trace = np.pad(trace, [1,1], 'constant')
    #print this_func, 'abc', abcissas
    #print this_func, 'trace', trace
    xlabel = 'frequency [Hz]'
    ylabel = 'amplitude [-]'
    
  elif spectrum == 'power':
    trace = abs(np.fft.fft(trace))**2
    xlabel = 'frequency [Hz]'
    ylabel = 'power [-]'
    
  elif spectrum == 'phase':
    trace = np.angle(np.fft.fft(trace))
    xlabel = 'frequency [Hz]'
    ylabel = 'phase [rad]'
  
  else: # TIME DOMAIN
    abcissas = np.arange(len(trace)) * dt
  
  x_axis = np.zeros(len(trace)) # X-axis arrow (y=0) 
  #print this_func, 'x_axis', x_axis
  #plt.scatter(abcissas, trace)
  #plt.plot(abcissas, '.')
  #plt.plot(trace, '-.')
  #plt.figure()
  Plot_2_Series(abcissas, x_axis, trace, c1='b', c2='r', xlabel=xlabel, ylabel=ylabel, **kwargs)
  
  # MARK CURRENT SAMPLE E.G. FOR WAVEFIELD EVOLUTIONARY SEQUENCE
  if mark_sample >= 0:
    plt.scatter([abcissas[mark_sample]], [trace[mark_sample]], s=50, c='k')#, facecolors='k', edgecolors='k')

  #print this_func + 'END'
  return abcissas, trace


# -------------------------------------------------------------------------------


def Plot_Picks(fname_vtr, **kwargs):
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
  
  this_func = this_lib + 'Plot_Picks: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  fname_asc = fname_vtr[:-len('.vtr')] + '_t1.ascii'
  c = Read_File(fname_asc)
  picks = [float(i[0]) / 2.5 for i in c]
  print(len(picks))
  marker_style = dict(color = 'lime', linestyle = 'none', marker = 'o', fillstyle = 'none',
                      markersize = 10)
  plt.gca().plot(picks, **marker_style)


  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
# PLOTT FILES OF DIFFERENT FORMATS
# -------------------------------------------------------------------------------


def Plot_sgy(fname, **kwargs):
  """
  Framework function to plot gridded, 
  3D data stored in a .sgy file
  through conversion to .vtr.
  
  Parameters
  ----------
  fname : file name 
    It must include extension.
    It can include path if needed. 
  **kwargs : keyword arguments, optional
      Just passing it down.

  Returns
  -------
  0
  
  Notes
  -----
  .vtr is a Fullwave3D's native binary format. 
  
  """    
  this_func = this_lib + 'Plot_sgy: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_io_fullwave import Read_sgy
  from lib_generic_PLOTT import Plot_Data
  
  Z = Read_sgy(fname, **kwargs)
  
  if len(Z) != 0:
    Plot_Data(Z, **kwargs)
  
  else: 
    fname_vtr = Strip(fname, **kwargs) + '.vtr'
    eprint(this_func + 'File in .vtr format: ' + fname_vtr + 
           ' not found. Skipping this figure. You can provide nx to automatically convert sgy 2 vtr\n')  
  
  #print this_func + 'END'
  return 0
  

# -------------------------------------------------------------------------------


def Plot_su(fname, **kwargs):
  """
  Framework function to plot gridded, 
  3D data stored in a .su file
  through conversion to .vtr.
  
  Parameters
  ----------
  fname : file name 
    It must include extension.
    It can include path if needed. 
  **kwargs : keyword arguments, optional
      Just passing it down.

  Returns
  -------
  0
  
  Notes
  -----
  .vtr is a Fullwave3D's native binary format. 
  
  """    
  this_func = this_lib + 'Plot_su: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  o = Bash('convert_su2sgy.sh ' + fname)
  fname_sgy = fname[:-len('.su')] + '.sgy'
  
  Plot_sgy(fname_sgy, **kwargs)
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_ttr(fname, **kwargs):
  """
  Through ttr->vtr
  
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
  
  this_func = this_lib + 'Plot_ttr: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_fwi_generic_PLOTT import Plot_vtr
  
  ext = fname.split(".")[-1]
  try:
    fname_core = kwargs['fname_core']
  except KeyError:
    fname_core = fname[ :-(len(ext)+1)]
  
  # FIXME: IT ASSUMES proj_name-Rest
  proj_name = Kwarg('proj_name', Split(fname_core, '-')[0], kwargs) 
  
  o = Bash('convert_ttr2vtr ' + proj_name)
  
  fname = fname[ :-len(ext)] + 'vtr'
  if verbos > 0:
    print(this_func, 'fname: ', fname)
  
  Plot_vtr(fname, **kwargs)

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_vtr(fname, **kwargs):
  """
  Framework function to plot gridded, 
  3D data stored in a .vtr file.
  
  Parameters
  ----------
  fname : file name 
    It must include extension.
    It can include path if needed. 
  **kwargs : keyword arguments, optional
      Just passing it down.

  Returns
  -------
  0
  
  Notes
  -----
  .vtr is a Fullwave3D's native binary format. 
  
  """    
  this_func = this_lib + 'Plot_vtr: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_fwi_generic import Read_vtr
  from lib_generic_PLOTT import Plot_Data
  
  
  file2 = Kwarg('file2', '', kwargs)

  # PROBABLY EVERY Plot_format NEEDS IT 
  # UNLESS EVERYTHING IS CONVERTED TO Plot_vtr
  fnames = [fname] + [file2]
  Z_list = []
  for fname in fnames:
    if fname != '':
      nx1, nx2, nx3, Z = Read_vtr(fname, **kwargs)
      Z_list.append(Z)
  
  if file2 != '':
    Plot_Juxtaposition(Z_list, **kwargs)
  
  else: 
    Plot_Data(Z_list[0], **kwargs)
    
  #print this_func + 'END'  
  return 0


# -------------------------------------------------------------------------------

