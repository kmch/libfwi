"""
Author: Kajetan Chrapkiewicz, 2019. 
All rights reserved. Ask for permision writing to k.chrapkiewicz17@imperial.ac.uk.

This library provides framework procedures for...
At the moment fully implemented are:

1. ...

"""

## MODULES
import numpy as np
import matplotlib.pyplot as plt

from lib_generic_CONST import *
from lib_generic import *
#from lib_fwi_generic import Read_vtr, Save_vtr
from lib_io_fullwave import *
from lib_fwi_project import Project_Filenames

from lib_generic_PLOTT import *

## CONSTS
this_lib = 'lib_fwi_project_PLOTT.py/'

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
# PLOT ALL PROJECT FILES
# -------------------------------------------------------------------------------


def Plot_Project(proj_name, **kwargs):
  """
  Plot throughput (input + output)
  of a Fullwave3D project for QC.
  
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
  this_func = this_lib + 'Plot_Project: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  Plot_Project_Input(proj_name, **kwargs)
  Plot_Project_Output(proj_name, **kwargs)
   
  #print this_func + 'END'  
  return 0


# -------------------------------------------------------------------------------
# PLOTT PROJECT'S INPUT
# -------------------------------------------------------------------------------


def Plot_Project_Input(proj_name, **kwargs):
  """
  Plot input of a Fullwave3D project for QC.
  
  Parameters
  ----------
  proj_name : str
    Project name assumed to be provided 
    as a first argument of the FWI code
    being a prefix of all FWI input files.  
  **kwargs : keyword arguments, optional
      Just passing it down.
      
  proj_path = Kwarg('proj_path', './', kwargs)  
  prefix = proj_path + proj_name
  
  plot_mod = Kwarg('plot_mod', False, kwargs)
  plot_sr = Kwarg('plot_sr', False, kwargs)
  plot_wavelet = Kwarg('plot_wavelet', False, kwargs)
  plot_data = Kwarg('plot_data', False, kwargs)
  plot_fs = Kwarg('plot_fs', False, kwargs)      

  Returns
  -------
  0
  
  Notes
  -----
  
  """    
  this_func = this_lib + 'Plot_Project_Input: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  # CHECK KEYWORD ARGUMENTS  
  proj_path = Kwarg('proj_path', './', kwargs)  
  prefix = proj_path + proj_name
  
  plot_mod = Kwarg('plot_mod', False, kwargs)
  plot_sr = Kwarg('plot_sr', False, kwargs)
  plot_sign = Kwarg('plot_sign', False, kwargs)
  plot_data = Kwarg('plot_data', False, kwargs)
  plot_fs = Kwarg('plot_fs', False, kwargs)
  
  #nx, ny = 2, 2
  #fig, ax = plt.subplots(nx, ny, figsize=[15, 10])
  #i = 1
  
  # THESE ARE SEPERATE FIGURES. OVERLAYING IS POSSIBLE 
  # INSIDE BACKEND FUNCTIONS
  if plot_mod:  
    #plt.subplot(nx,ny,i)
    #i += 1
    #try:
    Plot_Project_Input_Models(proj_name, **kwargs)
    #except:
    #eprint('Errors in Plot_Project_Input_Models')
  
  if plot_fs:
    plt.figure()
    Plot_Project_Input_FS(proj_name, **kwargs) 
  
  if plot_sr:
    plt.figure()
    Plot_Project_Input_SR(proj_name, **kwargs)
    
  if plot_sign:
    plt.figure()
    Plot_Project_Input_Wavelets(proj_name, **kwargs)
  
  if plot_data:
    plt.figure()
    Plot_Project_Input_Data(proj_name, **kwargs)

  #print this_func + 'END'  
  return 0


# -------------------------------------------------------------------------------


def Plot_Project_Input_Models(proj_name, **kwargs):
  """
  Plot input of a Fullwave3D project for QC.
  
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
  NOTE: it will need some care when superimposing 
  other stuff.
  
  """    
  this_func = this_lib + 'Plot_Project_Input_Models: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_fwi_project import Project_Suffices
  
  pars_gen = Kwarg('pars_gen', {}, kwargs)
  io = Kwarg('io', 'fw3d', pars_gen)
  
  pars_mod = Kwarg('pars_mod', {}, kwargs)
  layout = Kwarg('layout', '2+1', pars_mod) # SET DEFAULT TO '2+1'
  if 'layout' in kwargs:
    del kwargs['layout']
  
  suffices = Project_Suffices(io)
  
  # LOOP OVER ALL POSSIBLE MODEL-PARAMETERS
  for key in ['TrueVp']:#, 'TrueVs', 'epsi', 'delta', 'q']:
      fname = proj_name + suffices[key]
      Plot_Project_Input_Model(fname, proj_name, layout=layout, **kwargs)
      ##plt.gcf().suptitle('Model of ' + key)
      #Plot_Project_Input_Model(proj_name, key=key, **kwargs)

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Project_Input_Model(fname, proj_name, **kwargs): # NOTE: ONLY vtr
  """
  Plot input of a Fullwave3D project for QC.
  
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
  NOTE: it will need some care when superimposing 
  other stuff.
  
  """    
  this_func = this_lib + 'Plot_Project_Input_Model: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_fwi_generic_PLOTT import Plot_Model
  from lib_fwi_project import Project_Suffices
  
  pars_gen = Kwarg('pars_gen', {}, kwargs)
  io = Kwarg('io', 'fw3d', pars_gen)
  
  plot_type = Kwarg('plot_type', 'slices_xyz', kwargs)
  try:
    del kwargs['plot_type']
  except KeyError:
    pass
  
  # READ FS
  #proj_fnames = Project_Filenames(proj_name)
  fname_fs = proj_name + Project_Suffices(io)['FS']
  n1, n2, n3, surf = Read_vtr(fname_fs, **kwargs)
  
  # PLOT MODEL + FS + ...
  n1, n2, n3, model = Read_vtr(fname, **kwargs)
  Plot_Model(model, surfs=[surf], plot_type=plot_type, **kwargs)
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Project_Input_FS(proj_name, **kwargs):
  """
  NOTE: GONNA BE A BIG ONE.
  
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
  this_func = this_lib + 'Plot_Project_Input_FS: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  #from lib_fwi_fs import Read_GhostData_File_txt#, Index_Ghost_From_Coords
  #from lib_fwi_fs import Plot_Ghost_Data
  #from lib_fwi_generic import Plot_FS
  
  pars_FS = Kwarg('pars_FS', {}, kwargs)
  layout = Kwarg('layout', '2+1', pars_FS) # SET DEFAULT TO '2+1'
  if 'layout' in kwargs:
    del kwargs['layout']
  
  
  #fig, ax, i = Subplots(1, 3)
  
  # ORIGINAL FS
  #i = Subplot(fig, i)
  fname = proj_name + '-FreeSurf.vtr'
  Plot_File(fname, data_type='FS', plot_type='slices_xyz', layout=layout, **kwargs)
  
  #try: 
  #  #EXTEN. FS
  #  i = Subplot(fig, i)
  #  fname = proj_name + '-FreeSurf_exten.vtr'
  #  Plot_Gridded_Data_From_File(fname, data_type='FS', plot_type='map', **kwargs)
  #  
  #  # EXTEN. & INTERP. FS
  #  
  #  # FS ROUGHNESS
  #  nx1_fs, nx2_fs, nx3_fs, fs_z = Read_vtr(fname, **kwargs)
  #  i = Subplot(fig, i)
  #  Plot_FS_Roughness(fs_z, **kwargs) # DETERMINE FRACTAL DIMENSION   
  #
  #except:
  #  if verbos > 2:
  #    eprint('File ' + proj_name + '-FreeSurf_exten.vtr. Run fsprep.')
  #  pass
  
  # MODEL NOTE: WE DON'T HAVE THE EXTENDED VERSION THOUGH)
  #fname = prefix + '-TrueVp.vtr'  
  #nx1_m, nx2_m, nx3_m, model = Read_vtr(fname)

  # GHOST DATA
  #fname = prefix + '-GhostData.txt'
  #Plot_Ghost_Data(fname, fs_z=fs_z, **kwargs) # PROHIBITIVE FOR REAL 3D IN THE CURRENT SHAPE  
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Project_Input_SR(proj_name, **kwargs):
  """
  
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
  this_func = this_lib + 'Plot_Project_Input_SR: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_fwi_generic_PLOTT import Plot_SR
  
  pars_SR = Kwarg('pars_SR', {}, kwargs)
  xlim = Kwarg('xlim', [0, 100], pars_SR) # FIXME: READ nx1, ...
  ylim = Kwarg('ylim', [100, 0], pars_SR)
  #ext = Kwarg('ext', 'pgy', kwargs)
  layout = Kwarg('layout', '2+1', pars_SR) # SET DEFAULT TO '2+1'
  if 'layout' in kwargs:
    del kwargs['layout']
    
  Plot_SR(proj_name, ext='pgy', s_R=150, s_S=50, xlim=xlim, ylim=ylim, layout=layout, **kwargs)

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Project_Input_Wavelets(proj_name, **kwargs):
  """
  Plot input of a Fullwave3D project for QC.
  
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
  this_func = this_lib + 'Plot_Project_Input_Wavelets: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_fwi_generic import Runfile_Read
  from lib_fwi_project import Project_Suffices
    
  pars_gen = Kwarg('pars_gen', {}, kwargs)
  io = Kwarg('io', 'fw3d', pars_gen)  
  fname = proj_name + Project_Suffices(io)['RawSign']
  
  #runfile_name = proj_name + interfix_runfile + '.' + ext_runfile
  #runfile = Runfile_Read(runfile_name)
  
  #print runfile
  
  if True:
    Plot_Project_Input_Wavelet(fname, **kwargs)
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Project_Input_Wavelet(fname, **kwargs):
  """
  Plot input of a Fullwave3D project for QC.
  
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
  this_func = this_lib + 'Plot_Project_Input_Wavelet: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_fwi_generic_PLOTT import Plot_Wavelet_Full
  
  pars_gen = Kwarg('pars_gen', {}, kwargs)
  dt = Kwarg('dt', None, pars_gen)
  if not dt:
    raise KeyError('You need to provide dt')
  
  pars_src = Kwarg('pars_src', {}, kwargs)
  f_max = Kwarg('f_max', f_max_default, pars_src)
  
  # PLOT
  Plot_Wavelet_Full(fname, dt, f_max=f_max, **kwargs)

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Project_Input_Data(proj_name, **kwargs):
  """
  NOTE: GONNA BE A BIG ONE.
  
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
  this_func = this_lib + 'Plot_Project_Input_Data: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
# PLOT PROJECT'S OUTPUT
# -------------------------------------------------------------------------------


def Plot_Project_Output(proj_name, **kwargs):
  """
  Plot output of a Fullwave3D project.
  
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
  this_func = this_lib + 'Plot_Project_Output: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  # CHECK KEYWORD ARGUMENTS  
  data = Kwarg('data', True, kwargs)  
  fw = Kwarg('fw', True, kwargs)
  SR = Kwarg('SR', False, kwargs)

  # DATA
  if SR:
    plt.figure()
    Plot_Project_Input_SR(proj_name, **kwargs)  
  
  if data: 
    Plot_Project_Output_Data(proj_name, **kwargs)
  
  # FORWARD WAVEFIELD
  if fw:
    Plot_Project_Output_Wavefield(proj_name, **kwargs)

  #print this_func + 'END'  
  return 0


# -------------------------------------------------------------------------------


def Plot_Project_Output_Data(proj_name, **kwargs):
  """
  Plot synthetic data.
  
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
  this_func = this_lib + 'Plot_Project_Output_Data: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  proj_path = Kwarg('proj_path', './', kwargs)  
  fname = proj_path + proj_name + '-Observed-Time.ttr'
  
  trace1 = True # FIXME: TMP
  if trace1:
    Plot_Project_Output_Data_1Trace(proj_name, **kwargs)
  
  else:
    #Project_Output_Plot_Gathers #NOTE: CHECK IT!!!
    Plot_Gridded_Data_From_File(fname, data_type='model', cmap='seismic', **kwargs) # FIXME: TMP
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Project_Output_Data_1Gather(proj_name, **kwargs):
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
  this_func = this_lib + 'Plot_Project_Output_Data_1Trace: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Project_Output_Data_1Trace(proj_name, **kwargs):
  """
  FIXME: ONLY TEMPORARILY (INCORPORATE INTO Plot_Gather)
  
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
  this_func = this_lib + 'Plot_Project_Output_Data_1Trace: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  proj_path = Kwarg('proj_path', './', kwargs)  
  fname = proj_path + proj_name + '-Observed-Time.ttr'

  # SET UP A FIGURE
  nx, ny = 1, 2
  fig, ax = plt.subplots(nx, ny, figsize=[20, 8])
  i = 1
  # PLOT
  plt.subplot(nx,ny,i)
  i += 1
  plt.title('p')
  Plot(fname, data_type='wavelet', **kwargs) # FIXME
  plt.subplot(nx,ny,i)
  i += 1
  plt.title('dp/dt')
  Plot(fname, data_type='wavelet', deriv=1, **kwargs) 
  #Plot_Gridded_Data_From_File(fname, data_type='model', plot_type='slice', cmap="seismic")


  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Project_Output_Wavefield(proj_name, **kwargs):
  """
  Plot wavefield snapshots.
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  NOTE: shot_ids needs to be a list!
  
  """
  this_func = this_lib + 'Plot_Project_Output_Wavefield: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_fwi_generic import SR_Read_Shot_IDs
  
  
  
  # READ SHOTS
  shot_ids = SR_Read_Shot_IDs(proj_name, **kwargs)
  shot_ids = Kwarg('shot_ids', [shot_ids[0]], kwargs) # PLOTT FIRST SHOT 
  #shot_id = Kwarg('shot_id', '*', kwargs)
  
  # ITERATE OVER SHOTS
  for shot_id in shot_ids:
    Plot_Project_Output_Wavefield_1Shot(proj_name, shot_id, **kwargs)
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Plot_Project_Output_Wavefield_1Shot(proj_name, shot_id, **kwargs):
  """
  Plot wavefield snapshots.
  
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
  this_func = this_lib + 'Plot_Project_Output_Wavefield_1Shot: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_fwi_generic import Set_Extra_Nodes, SR_Read_Shot_IDs
  from lib_fwi_generic import Plot_Wavefield_Propagation  
  
  proj_path = Kwarg('proj_path', './', kwargs)
  prefix = proj_path + proj_name
    
  
  # READ FILE NAMES
  fnames = Get_Files(proj_path, proj_name + '-fw-*0' + shot_id + '-*vtr') #FIXME: THIS *0
  if len(fnames) == 0: # NOTE: TO PLOTT WHILE STILL RUNNING
    fnames = Get_Files(proj_path, 'fw-*0' + shot_id + '-*vtr') 
  
  # SET EXTRA NODES TO SHIFT PLOTS
  etop, eleft, efront, ebot, eright, eback = Set_Extra_Nodes(prefix + '-Runfile.key', **kwargs)
  if verbos > 4:
    print(this_func, 'eleft, efront, etop', eleft, efront, etop)
  
  # CALL THE ADVANCED PLOTTER
  Plot_Wavefield_Propagation(proj_name, fnames, shift_x=eleft, shift_y=efront, shift_z=etop, **kwargs) 
  
  #print this_func + 'END'
  return 0  


# -------------------------------------------------------------------------------

