"""
Author: Kajetan Chrapkiewicz, 2018. 
All rights reserved. Ask for permision writing to k.chrapkiewicz17@imperial.ac.uk.

This library provides procedures automating 
work with a single FWI project.

"""

## MODULES
import numpy as np
import matplotlib.pyplot as plt

## MY MODULES
from lib_generic import *
from lib_fwi_project_CONST import io_default, dx_default, z_sea_default
##


## CONSTS
this_lib = 'lib_fwi_project.py/'


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
# INPUT
# -------------------------------------------------------------------------------


def Project_Input(proj_name, **kwargs):
  """
  NOTE: REPLACE WITH AN OBJECT project.input.
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  These all should be objects (OOP).
  
  """
  this_func = this_lib + 'Project_Input: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  prep = Kwarg('prep', False, kwargs)
  check = Kwarg('check', False, kwargs)
  plot = Kwarg('plot', False, kwargs)
  
  if prep:
    Project_Input_Prepare(proj_name, **kwargs)
  
  if check:
    Project_Input_Check(proj_name, **kwargs)
  
  if plot:
    Project_Input_Plot(proj_name, **kwargs)

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Project_Input_Prepare_All(proj_name, **kwargs):
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
  Creating a new project usually makes sense
  only if it requires re-running Segyprep.
  
  """
  this_func = this_lib + 'Project_Input_Prepare_All: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from . import lib_fwi_project_CLASS
  from lib_fwi_generic import Segyprep_Input_Prepare, Segyprep_Run
  
  pars_gen = Kwarg('pars_gen', {}, kwargs)
  io = Kwarg('io', io_default, pars_gen) 
  dx = Kwarg('dx', dx_default, pars_gen) 
  pars_FS = Kwarg('pars_FS', {}, kwargs)
  z_sea = Kwarg('z_sea', z_sea_default, pars_FS)
  overwrite = Kwarg('overwrite', 0, pars_FS)
  
  # SEGYPREP FILES
  Segyprep_Input_Prepare(proj_name, **kwargs)
  Segyprep_Run(proj_name, z_sea, dx, **kwargs)
  
  # OTHER FILES
  pfiles = [
    lib_fwi_project_CLASS.FS_File(proj_name, io),
    lib_fwi_project_CLASS.Runfile(proj_name, io),
    lib_fwi_project_CLASS.PBS_File(proj_name, io)
  ] 
  
  for pfile in pfiles:
    pfile.Prepare(**kwargs) 
    
    
  # RUN FSPREP
  o = Bash('fsprep.sh ' + proj_name + ' ' + str(overwrite))
  if verbos > 5:
    print(o)

  #print this_func + 'END'
  return 0


# CHECK -------------------------------------------------------------------------


def Project_Input_Check(proj_name, dx, dt, vel, **kwargs):
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
  Creating a new project usually makes sense
  only if it requires re-running Segyprep.
  
  """
  this_func = this_lib + 'Project_Input_Check: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_math_num import Check_Stability, Check_Accuracy, Max_Resolution, Propagation_Dists
  
  proj_path = Kwarg('proj_path', './', kwargs)
  f_max = Kwarg('f_max', None, kwargs)
  
  v_min = Kwarg('v_min', None, kwargs)
  v_max = Kwarg('v_max', None, kwargs)
  
  if not v_min:
    v_min = vel 
  
  if not v_max:
    v_max = vel
  
  
  
  fname = proj_path + proj_name + '-NumChecks.log'
  f = open(fname, 'w')
  
  Check_Stability(f, dx, dt, v_max, 'low')
  if f_max:
    Check_Accuracy(f, dx, v_min, f_max, 'low')
    Max_Resolution(f, dx, v_min, f_max)
    f_min = 0.1 # Hz # FIXME
    
    #Propagation_Dists(f, dx, nx1, nx2, nx3, vel_min, vel_max, f_min, f_max)
  
  f.close()
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Project_Input_Check_Dims(fs_file, mod_file, **kwargs):
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
  Creating a new project usually makes sense
  only if it requires re-running Segyprep.
  
  """
  this_func = this_lib + 'Project_Input_Check_Dims: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_io_fullwave import Read_vtr
  
  f1, f2, f3, fs = Read_vtr(fs_file, **kwargs)
  m1, m2, m3, mo = Read_vtr(mod_file, **kwargs)
  
  if (f1 != m1):
    raise ValueError('No. of inline nodes different for ' + fs_file + ': ' + str(f1) + ', and ' + mod_file + ' ' + str(m1))
  
  if (f2 != m2):
    raise ValueError('No. of x-line nodes different for ' + fs_file + ': ' + str(f2) + ', and ' + mod_file + ' ' + str(m2))
  
  if verbos > 0:
    print(this_func, 'Check successful - dimensions of ' + mod_file + ' match ' + fs_file) 
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
# OUTPUT
# -------------------------------------------------------------------------------


def Project_Output(proj_name, **kwargs):
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
  These all should be objects (OOP).
  
  """
  this_func = this_lib + 'Project_Output: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  prep = Kwarg('prep', True, kwargs)
  check = Kwarg('check', True, kwargs)
  plot = Kwarg('plot', True, kwargs)  
  
  if prep:
    Project_Output_Prepare(proj_name, **kwargs)
  
  if check:
    Project_Output_Check(proj_name, **kwargs)
  
  if plot:
    Project_Output_Plot(proj_name, **kwargs)

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Project_Output_Prepare(proj_name, **kwargs):
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
  These all should be objects (OOP).
  
  """
  this_func = this_lib + 'Project_Output_Prepare: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Project_Output_Check(proj_name, **kwargs):
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
  These all should be objects (OOP).
  
  """
  this_func = this_lib + 'Project_Output_Check: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Project_Output_Plot(proj_name, **kwargs):
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
  These all should be objects (OOP).
  
  """
  this_func = this_lib + 'Project_Output_Plot: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  # SPATIAL PROPERTIES
  #fnames = Get_Files(proj_path, proj_name+'-*CP*TrueVp.sgy')
  #Plot_Evolution(proj_name, fnames=fnames, data_type='model', **kwargs)
  #Plot_Evolution(proj_name, fnames=fnames, data_type='dm', **kwargs)
  
  #fnames = Get_Files(proj_path, proj_name+'-*CP*RawGrad.sgy')
  #Plot_Evolution(proj_name, fnames=fnames, data_type='grad', **kwargs)
  
  #fnames = Get_Files(proj_path, proj_name+'-*CP*Grad.sgy')
  #Plot_Evolution(proj_name, fnames=fnames, data_type='grad', **kwargs)
  
  #fnames = Get_Files(proj_path, proj_name+'-*CP*Precond.sgy')
  #Plot_Evolution(proj_name, fnames=fnames, data_type='grad', **kwargs)
  


  # TEMPORAL PROPERTIES
  #Project_Output_Plot_Gathers(proj_name, **kwargs) # NOW ONLY FOR START MOD 
  # ONCE WE DAMP .ttr AT EVERY TIME-STEP, WE WILL PLOTT EVOLUTION OF IT TOO



  # ALL COMBINED INTO WAVEFIELD SNAPSHOTS
  #fnames = Get_Files(proj_path, proj_name+'-fw*StartVp.sgy', **kwargs) # START MOD
  #Plot_Propagation(proj_name, fnames=fnames, data_type='fw', plot_fs=True, plot_mod=True, plot_sr=True, plot_data=True, plot_wavelet=True)
  
  #fnames = Get_Files(proj_path, proj_name+'-fw*FinalVp.sgy', **kwargs) # FINAL MOD
  #Plot_Propagation(proj_name, fnames=fnames, data_type='fw', plot_fs=True, plot_mod=True, plot_sr=True, plot_data=True, plot_wavelet=True)
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Project_Output_Plot_Gathers(proj_name, **kwargs):
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
  These all should be objects (OOP).
  
  """
  this_func = this_lib + 'Project_Output_Plot_Gathers: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_fwi_generic import Read_vtr
  from lib_su import SU_Get_Traces_No
  from lib_math_signal import Gather_Modify_Normalize
  from lib_generic_PLOTT import Plot_Gridded_Data#_From_File
  
  
  # PLOTT OBS. AND SYNTH. GATHERS SEPARATELY
  gathers = []
  for suffix in ['-OutSeis.sgy', '-Synthetic.sgy']:
    fname = proj_name + suffix
    fname_vtr = fname[:-len('.sgy')] + '.vtr'
    try:
      nx1, nx2, nx3, g = Read_vtr(fname_vtr)
    except:
      ntraces = SU_Get_Traces_No(fname)
      o = Bash('sgy2vtr.sh ' + fname + ' ' + str(ntraces))
      nx1, nx2, nx3, g = Read_vtr(fname_vtr)
    #max_abs = max(abs(np.min(g)), abs(np.max(g)))
    g_normed = Gather_Modify_Normalize(g)
    #g_normed = g / float(max_abs)
    #Plot_Gridded_Data(g_normed, cmap='seismic', fig_size=[15,8], plt_func='imsh')
    gathers.append(g_normed)

  # JUXTAPOSE THE TWO
  #Plot_Gridded_Data(gathers[0], interleave=[gathers[1]], cmap='seismic', fig_size=[15,8], plt_func='imsh')

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
# GENERIC
# -------------------------------------------------------------------------------


def Project_Duplicate(pnew, pold, files2dupl, **kwargs):
  """
  Copy / link files from an old project 
  into a new one.
  
  Parameters
  ----------
  pnew : Proj
    New project.
  pold : Proj 
    Project to duplicate.
  
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Project_Duplicate: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  #print this_func, 'pnew', pnew.name
  #print this_func, 'pold', pold.name  
  
  #fnames = pnew.files2dupl # THESE ARE FILES OF OLD PROJECT
  
  for fname in files2dupl:
    print(fname)
    #fname = Path_Leave(fname, **kwargs)
    #print fname
    #core = Core(fname, pold.name, None, **kwargs)
    #nfname = pnew.name + core
    
 #   # (RE-)ADD PATHS
 #   fname = pold.path + fname
 #   nfname = pnew.path + nfname
 #   #print this_func, 'fname', fname
 #   #print this_func, 'nfname', nfname
 #   
 #   
 #   # DO THE JOB
 #   Duplicate(fname, nfname, **kwargs)
 # 
 # # FIXME
 # pbs_file = pnew.name + pnew.suffices['pbs']
 # o = Bash('sed -i -e s/' + pold.name + 
 #          '/' + pnew.name + '/g ' + pbs_file)
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
# META-DATA
# -------------------------------------------------------------------------------


def Project_Suffices(io, **kwargs):
  """
  Suffices of all 
  
  Parameters
  ----------
  io : str 
    API for IO: fw3d, segy, etc.
    
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  fnames : dict 
    Dictionary of fnames.
  
  Notes
  -----
  NOTE Needs to be completed
  - epsilon, delta
  - Vs
  - other damped 'wavefields' (pre-conditioner etc.)
  
  """
  this_func = this_lib + 'Project_Suffices: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  #pars_gen = Kwarg('pars_gen', {}, kwargs)
  #io = Kwarg('io', 'fw3d', kwargs)# pars_gen)
  
  suffices = {}
  if io == 'fw3d':
    suffices['TrueVp'] = '-TrueVp.vtr'
    suffices['TrueVs'] = '-TrueVs.vtr'
    suffices['TrueDelta'] = '-TrueDelta.vtr' # D-C
    suffices['TrueEpsil'] = '-TrueEpsilon.vtr' # D-C
    suffices['TrueQp'] = '-TrueQp.vtr' # D-C
    suffices['TrueQs'] = '-TrueQs.vtr' # D-C      
    suffices['StartVp'] = '-StartVp.vtr'
    suffices['StartVs'] = '-StartVs.vtr'
    suffices['StartDelta'] = '-StartDelta.vtr' # D-C
    suffices['StartEpsil'] = '-StartEpsilon.vtr' # D-C
    suffices['StartQp'] = '-StartQp.vtr' # D-C 
    suffices['StartQs'] = '-StartQs.vtr' # D-C     
    suffices['S'] = '-PointSources.pgy'
    suffices['R'] = '-PointReceivers.pgy'
    suffices['Sign'] = '-SourceSig-Time.ttr'
    suffices['SignIdx'] = 'dummy' # NOTE
    suffices['Obser'] = '-Observed-Time.ttr'
    suffices['ObserIdx'] = '-Observed-0000.ttr'
    suffices['Synth'] = suffices['Obs'] 
    suffices['SynthIdx'] = suffices['ObsIdx']     
    suffices['RawSeis.sgy'] = 'dummy' # NOTE   
    suffices['CP_Vp_files'] = '-CP*-Vp.vtr' # D-C
    suffices['Templ'] = 'dummy'
    suffices['TemplIdx'] = 'dummy'
    suffices['OutSeis'] = 'dummy'
    
    
  elif io == 'sgy':
    suffices['TrueVp'] = '-TrueVp.sgy'
    suffices['TrueVs'] = '-TrueVs.sgy'
    suffices['TrueDelta'] = '-TrueDelta.sgy' # D-C
    suffices['TrueEpsil'] = '-TrueEpsilon.sgy' # D-C
    suffices['TrueQp'] = '-TrueQp.sgy' # D-C
    suffices['TrueQs'] = '-TrueQs.sgy' # D-C    
    suffices['StartVp'] = '-StartVp.sgy'
    suffices['StartVs'] = '-StartVs.sgy'    
    suffices['StartDelta'] = '-StartDelta.sgy' # D-C
    suffices['StartEpsil'] = '-StartEpsilon.sgy' # D-C
    suffices['StartQp'] = '-StartQp.sgy' #NOTE: THIS IS INVERTED FOR ONLY BY NUNO'S qPSO, I THINK
    suffices['StartQs'] = '-StartQs.sgy' #NOTE: THIS IS INVERTED FOR ONLY BY NUNO'S qPSO, I THINK    
    suffices['S'] = '-Sources.geo'
    suffices['R'] = '-Receivers.geo'
    suffices['Sign'] = '-Signature.sgy'
    suffices['SignIdx'] = '-Signature.idx'
    suffices['Obser'] = '-Observed.sgy'
    suffices['ObserIdx'] = '-Observed.idx'
    suffices['Synth'] = '-Synthetic.sgy'
    suffices['SynthIdx'] = '-Synthetic.idx' 
    suffices['RawSeis.sgy'] = '-RawSeis.sgy' 
    suffices['CP_Vp_files'] = '-CP*-Vp.sgy'    
    suffices['Templ'] = '-Template.sgy'
    suffices['TemplIdx'] = '-Template.idx'
    suffices['OutSeis'] = '-OutSeis.sgy'
   
 
  
  else:
    raise ValueError('Unknown IO API: ' + io)

  # GENERIC
  suffices['SynOut'] = '-SynOut.log'
  suffices['SynErr'] = '-SynErr.log'
  suffices['InvOut'] = '-InvOut.log'
  suffices['InvErr'] = '-InvErr.log'

  suffices['RawSign'] = '-RawSign.sgy'
  suffices['RawSeis.txt'] = '-RawSeis.txt' # SIC (!)  
  suffices['SegyPrep'] = '-SegyPrep.key'
  
  suffices['ModPrep'] = '-ModPrep.key'
  

  suffices['FS'] = '-FreeSurf.vtr'
  
  suffices['runfile'] = '-Runfile.key'
  suffices['runfile_synth'] = '-Runfile_synth.key'
  suffices['runfile_inv'] = '-Runfile_inv.key'
  suffices['runfile_skelet'] = '-Skeleton.key'
  suffices['pbs'] = '-Run.pbs'
  
  
  suffices['StartRawGrad'] = '-StartRawGrad.vtr' # NOTE: THIS IS CHEATING 
  suffices['StartRawPrec'] = '-StartRawPrec.vtr' # NOTE: THIS IS CHEATING 
  
  suffices['CP_Runfiles'] = '-CP*-Runfile.key'  
  #suffices['CP_RawGrad_files'] = '-CP*-RawGrad.vtr'    
  suffices['fw_files'] = '-fw-*.vtr'
  suffices['bw_files'] = '-bw-*.vtr'
  
  # EXTERNAL (META-) DATA
  suffices['data'] = '-Observed_*.sgy'
  suffices['PickP'] = '-tlPick_*P.dat' # POSSIBLE TO cat IT BY A DEFAULT CLASS METHOD
  
  # ANOMALIES
  suffices['Checker'] = '-Checker.vtr' # NOTE
  
  
  #print this_func + 'END'
  return suffices


# -------------------------------------------------------------------------------


def Project_Filenames(proj_name, io, mode, **kwargs):
  """
  Get the file names of the input/output 
  of the project using given I/O API.
  
  Parameters
  ----------
  io : str 
    API for IO: fw3d, segy, etc.
    
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  fnames : dict 
    Dictionary of fnames.
  
  Notes
  -----
  FIXME: ADD epsilon, delta etc.
  
  """
  this_func = this_lib + 'Project_Filenames: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  in_out = Kwarg('in_out', 'in', kwargs)
  if in_out != 'in' and in_out != 'out':
    raise ValueError('Wrong in_out: ' + in_out)

  all_suffices = Project_Suffices(io, **kwargs)
  suffices = {}  
  
  keys = []
  
  if mode == 'sp':
    if in_out == 'in':
      keys += ['SegyPrep',
              'RawSign',
              'TrueVp',
              'TrueVs',
              'RawSeis.txt']
    elif in_out == 'out':
      keys += ['S',
              'R',
              'Sign',
              'SignIdx',
              'OutSeis',
              'Templ', 
              'TemplIdx']
      
  if mode == 'synth' or mode == 'synthetic' or mode == 'both':
    if in_out == 'in':
      keys += ['TrueVp', 
              'TrueVs', 
              'S', 
              'R', 
              'Sign', 
              'SignIdx', 
              'Templ', 
              'TemplIdx', 
              'OutSeis', 
              'FS',
              'runfile', 
              'runfile_synth', 
              'runfile_skelet',
              'pbs'] 

    elif in_out == 'out':
      keys += ['Synth', 
              'SynthIdx',
              'SynOut',
              'SynErr']
    
  if mode == 'inv' or mode == 'inversion' or mode == 'tomography' or mode == 'both': 
    if in_out == 'in':
      keys += ['TrueVp', 
              'TrueVs', 
              'StartVp',
              'StartVs',
              'S', 
              'R', 
              'Sign', 
              'SignIdx', 
              'Templ', 
              'TemplIdx', 
              'FS',
              'runfile', 
              'runfile_inv', 
              'runfile_skelet',
              'pbs'] 

    elif in_out == 'out':
      keys += ['Synth', 
              'SynthIdx',
              'InvOut',
              'InvErr']

  #else:
    #raise ValueError('Unknown mode: ' + mode)

  
  fnames = {}
  for key in keys:
    fnames[key] = proj_name + all_suffices[key]
  
  #print this_func + 'END'
  return fnames


# -------------------------------------------------------------------------------


def Read_FWI_Project_Meta_Data(meta_path, **kwargs):
  """
  
  
  Parameters
  ----------
  
  
  
  Returns
  -------
  
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Read_FWI_Project_Meta_Data: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
   
  gx_list, gy_list, sx_list, sy_list, obs_list, shot_lines = Read_OBS_Data(meta_path, **kwargs)


  #print this_func + 'END'
  return gx_list, gy_list, sx_list, sy_list, obs_list, shot_lines


# -------------------------------------------------------------------------------


def Read_OBS_Data(meta_path, **kwargs):
  """
  
  
  Parameters
  ----------
  meta_path : str
    Path to meta-data files.
  stations_IDs_file : str
    
  
  
  Returns
  -------
  
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Read_OBS_Data: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_class_fwi import OBS, Shot_line
  
  try:
    stations_IDs_file = kwargs['stations_IDs_file']
  except KeyError:
    stations_IDs_file = "stations_IDs.txt"

  try:
    shot_info_file = kwargs['shot_info_file']
  except KeyError:
    shot_info_file = "Santorini_ShotLineTablewShotCorrections_8-19-2016.csv"

  # GET STATION IDS
  c = Read_File(meta_path + stations_IDs_file)
  station_ids = [i[0] for i in c]
  
  # FILL-IN THE DICTIONARY
  obs_list = {} # SHORT NAME FOR RAPID TYPING
  for i, station_id in enumerate(station_ids):   
    sensor = station_id[0]
    station_no = station_id[1: ] # JUST NUMBERS TO TAKE ADVANTAGE OF ARITHMETICS
    obs_list[station_no] = OBS(ID=station_no, sensor=sensor)
  
  
  c = Read_File(meta_path + 'gx.txt')
  gx_list = [float(i[0]) for i in c]
  
  c = Read_File(meta_path + 'gy.txt')
  gy_list = [float(i[0]) for i in c]

  print(this_func, 'No. of stations in current dataset: ', len(gx_list))


  c = Read_File(meta_path + 'sx.txt')
  sx_list = [float(i[0]) for i in c]
      
  c = Read_File(meta_path + 'sy.txt')
  sy_list = [float(i[0]) for i in c]  
  
  for i, station_id in enumerate(station_ids): #FIXME: IT USES station_ids BECAUSE gx_list ARE IN THIS ORDER   
    station_no = station_id[1: ]
    obs_list[station_no].x = gx_list[i]
    obs_list[station_no].y = gy_list[i]
  
  
  # SHOT INFO 
  shot_lines = {}

  content = Read_File(meta_path + shot_info_file) # READ .csv FILE (NOT GENERIC)
  for line in content:
    if len(line) < 2:
        continue 
    l = line[0].split(',')
    if len(l) < 5:
        continue 
    line_no = l[0]
    if line_no == '':
        continue
        
    shot_start = l[2]
    shot_end = l[3]
    
    # CREATE AN OBJECT INSTANCE
    shot_line = Shot_line(line_no, shot_start, shot_end)
    
    # ADD IT TO THE DICTIONARY
    shot_lines[line_no] = shot_line

  print(this_func, 'Total no. of shot lines len(shot_lines): ', len(shot_lines))
  
  #print this_func + 'END'
  return gx_list, gy_list, sx_list, sy_list, obs_list, shot_lines


# -------------------------------------------------------------------------------

