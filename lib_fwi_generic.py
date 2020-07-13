"""
Author: Kajetan Chrapkiewicz, 2018. 
All rights reserved. Ask for permision writing to K.Chrapkiewicz17@imperial.ac.uk.

This library provides framework procedures for...
At the moment fully implemented is:

1. ...

"""

## MODULES
import numpy as np
from lib_generic_CONST import *
from lib_fwi_project_CONST import *
from lib_generic import *
from lib_io_fullwave import Read_pgy, Read_geo, Save_pgy, Read_vtr
#from lib_fwi_project import Proj_File, Model_File
from lib_fwi_project_CONST import io_default, dx_default, z_sea_default

## CONSTS
this_lib = 'lib_fwi_generic.py/'

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
# DATA EXTRACTION
# -------------------------------------------------------------------------------


def Slice(X, **kwargs):
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
  this_func = this_lib + 'Slice: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  slice_type = Kwarg('slice_type', 'data', kwargs)
  
  if slice_type == 'data':
    X, labels = Slice_Data(X, **kwargs)
  
  elif slice_type == 'points':
    X, labels = Slice_Points(X, **kwargs)

  else:
    raise ValueError('Wrong slice_type.')
  
  #print this_func + 'END'
  return X, labels


# -------------------------------------------------------------------------------


def Slice_Data(ZZZ, **kwargs):
  """
  Only for .vtr-like data format!
  That's the only format we use throughout.
  Everything else should be converted to it!
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  FS will be collapsed into a curve.
  
  """
  this_func = this_lib + 'Slice_Data: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  slice_coord = Kwarg('slice_coord', slice_coord_default, kwargs) # Y IS MOST OFTEN COLLAPSED
  data_type = Kwarg('data_type', 'vol', kwargs)
  
  #if coord_value == 'all':
  #  i1 = 0
  #  i2 = -1
  #else:
  #  i1 = coord_value
  #  i2 = coord_value + 1
  
  if slice_coord == 'x':
    coord_value = Kwarg('coord_value', len(ZZZ)/2, kwargs)
    ZZ = ZZZ[coord_value]
    title = 'Slice for in-line node no. ' + str(coord_value)
    xlabel = 'cross-line node'
    ylabel = 'depth node'

  elif slice_coord == 'y':
    coord_value = Kwarg('coord_value', len(ZZZ[0])/2, kwargs)
    ZZ = [i[coord_value] for i in ZZZ]
    title = 'Slice for cross-line node no. ' + str(coord_value)
    xlabel = 'in-line node'
    ylabel = 'depth node'
  
  elif slice_coord == 'z':
    coord_value = Kwarg('coord_value', len(ZZZ[0][0])/2, kwargs)
    ZZ = [[j[coord_value] for j in i] for i in ZZZ]
    title = 'Slice for depth node no. ' + str(coord_value)
    xlabel = 'in-line node'
    ylabel = 'cross-line node'
    
  else:
    raise ValueError('Error. Wrong slice coord: ' + slice_coord + '\n')
  
  if data_type == 'FS': # GET RID OF THE DUMMY Z-DIMENSION    
    ZZ = [i[0] for i in ZZ]

  ZZ = np.array(ZZ)    

  labels = [title, xlabel, ylabel]
  
  #print this_func + 'END'
  return ZZ, labels


# -------------------------------------------------------------------------------


def Slice_Points(points, **kwargs):
  """
  Project points onto 1 of the planes.
  => coord_value is actually dummy!
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  This seems to be quite generic. 
  
  """
  this_func = this_lib + 'Slice_Points: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  slice_coord = Kwarg('slice_coord', 'y', kwargs)
  #coord_value = Kwarg('coord_value', len(ZZZ[0])/2, kwargs) # NOTE: IT'S ACTUALLY DUMMY  
  
  if slice_coord == 'x':
    X1, X2 = [[i[1] for i in points], [i[2] for i in points]]
    title = 'Projection onto YZ plane'
    xlabel = 'cross-line node'
    ylabel = 'depth node'    
    
  elif slice_coord == 'y':
    X1, X2 = [[i[0] for i in points], [i[2] for i in points]]
    title = 'Projection onto XZ plane'
    xlabel = 'in-line node'
    ylabel = 'depth node'    

  elif slice_coord == 'z':
    X1, X2 = [[i[0] for i in points], [i[1] for i in points]]
    title = 'Projection onto XY plane'
    xlabel = 'in-line node'
    ylabel = 'cross-line node'    

  else:
    raise ValueError('Error. Wrong slice coord: ' + slice_coord + '\n')
  
  points2d = list(zip(X1, X2))
  
  labels = [title, xlabel, ylabel]
  
  #print this_func + 'END'
  return points2d, labels


# -------------------------------------------------------------------------------


def SR_Read_Shot_IDs(proj_name, **kwargs):
  """
  Read IDs of all shots. Useful when we 
  want to plot the input/output for 
  the given shot(s).
  
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
  this_func = this_lib + 'SR_Read_Shot_IDsc: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  
  
  dims, sources, receivers = SR_Read(proj_name, adjust=False, **kwargs)
  
  if verbos > 4:
    print(this_func, sources)
  
  ids = []
  for key in sources:
    ids.append(key)

  #print this_func + 'END'
  return ids


# -------------------------------------------------------------------------------


def Wavefield_Extract_Shot_ID(fname, **kwargs):
  """
  Extract the shot ID from the filename 
  of the wavefield snapshot.
  
  Parameters
  ----------
  fname : str 
    File name of the wavefield snapshot. 
    It can contain a path. 
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  Returns
  -------
  shot_id : str 
    Shot ID.
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Wavefield_Extract_Shot_ID: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  shot_id = Path_Leave(fname)
  shot_id = Split(shot_id, 'csref')[1]
  shot_id = Split(shot_id, '-')[0]
  shot_id = shot_id.lstrip('0') # STRIP OFF LEADING ZEROS
  
  print(this_func, 'shot_id', shot_id)
  
  #print this_func + 'END'
  return shot_id


# -------------------------------------------------------------------------------


def Set_Extra_Nodes(runfile_name, **kwargs): #NOTE: IMPORTANT
  """
  Adjust sources & receivers coordinates 
  to conform to extended grid 
  (by extra nodes, see the runfile).
  
  Parameters
  ----------
  runfile_name : str
    Name of the runfile. 
    It should include  extension. 
    It can include path if needed.
  
  Returns
  -------
  extratop : float 
    No. of extra nodes added to the top edge.
  extraleft  : float 
    No. of extra nodes added to the left edge.
  extrafront  : float 
    No. of extra nodes added to the front edge.
    
  Notes
  -----
  NOTE: it must be consistent with what Fullwave3D does!
  It should also be synced with io_mod.f90/Set_Aux_Variables.
  
  """
  this_func = this_lib + 'Set_Extra_Nodes: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')       
  
  
  
  runfile = Runfile_Read(runfile_name, **kwargs)
  
  nx2 = int(runfile['NX2']) 
  
  try:
    # NOTE: WE NEED INT FOR ARRAY INDICES
    btop       = int(runfile['btop'])  
    bbot       = int(runfile['bbot'])  
    bleft      = int(runfile['bleft'])  
    bright     = int(runfile['bright'])  
    bfront     = int(runfile['bfront'])  
    bback      = int(runfile['bback'])  
    extratop   = int(runfile['extratop'])  
    extrabot   = int(runfile['extrabot'])  
    extraleft  = int(runfile['extraleft'])  
    extraright = int(runfile['extraright'])  
    extrafront = int(runfile['extrafront'])  
    extraback  = int(runfile['extraback']) 
  except KeyError:
    eprint(this_func + 'Error. Description of boundaries in the Runfile is not complete.\n')
    eprint(this_func + '(all 6 b* and 6 e* parameters are needed, even for 2D) \n')
    quit()
  
  if (btop == 0) and (extratop < 2):
    extratop += 2

  if (bbot == 0) and (extrabot < 2):
    extrabot += 2

  if (bleft == 0) and (extraleft < 2):
    extraleft += 2

  if (bright == 0) and (extraright < 2):
    extraright += 2
  
  # ONLY FOR 3D:
  if nx2 > 1:
    if (bfront == 0) and (extrafront < 2):
      extrafront += 2
    
    if (bback == 0) and (extraback < 2):
      extraback += 2    
  
  if verbos > 1:
    print(this_func, 'extratop, extraleft, extrafront', extratop, extraleft, extrafront)  
    print(this_func, 'extrabot, extraright, extraback', extrabot, extraright, extraback)  
  
    
  #print this_func + 'END'
  return extratop, extraleft, extrafront, extrabot, extraright, extraback


# -------------------------------------------------------------------------------
# SEGYPREP INPUT
# -------------------------------------------------------------------------------



# SEGYPREP RUN(FILE) ------------------------------------------------------------


def Segyprep_Create(proj, **kwargs):
  """
  Prepare a runfile of SegyPrep 
  arguments.
  
  Parameters
  ----------
  proj_name : str
    Project name assumed to be provided 
    as a first argument of the FWI code
    being a prefix of all FWI input files

  Returns
  -------
  segyprep : dict 
    Dictionary of SegyPrep parameters.
  
  Notes
  -----
  Add cross-line of receivers?
  
  """
  this_func = this_lib + 'Segyprep_Create: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
    
  from lib_generic import Save_Dict
  
  proj_path = proj.path + '/inp/'
  fname = proj_path + proj.name + '-SegyPrep.key'

  dims = proj.dims 
  nx1, nx2, nx3 = dims 
  dx = proj.dx 
  dt = proj.dt
  dt_ms = dt * 1000 # ms
  ns = proj.ns 
  ttime = proj.ttime 
  ttime_ms = ttime * 1000 # ms
  
  Ss = Kwarg('Ss', [], kwargs) 
  Rs = Kwarg('Rs', [], kwargs)
  fs_z_max = Kwarg('fs_z_max', nx3 / 2, kwargs)
  
  geometry = Kwarg('geometry', 'sgy', kwargs)
  reciprocity = Kwarg('reciprocity', 0, kwargs)
  
  segyprep = {}
  
  problem = Kwarg('problem', None, kwargs)
  if problem:
    segyprep['problem'] = problem
  else:
    segyprep['problem'] = proj.problem  
  
  io = Kwarg('io', None, kwargs)
  if io:
    segyprep['io api'] = io
  else:
    segyprep['io api'] = proj.io
  
    
  segyprep['nx1'] = nx1
  segyprep['nx2'] = nx2
  segyprep['nx3'] = nx3
  segyprep['dx'] = dx
  segyprep['ttime'] = ttime_ms # YES _ms!
  segyprep['dtms'] = dt_ms
  segyprep['reciprocity'] = reciprocity
  segyprep['geometry'] = geometry
  
  
  if segyprep['geometry'] == 'sgy':
    segyprep['geometry'] = 'segy'
    
    segyprep['x origin'] = proj.box[0] # FIXME: SHOULD BE box.x1 ETC.
    segyprep['x shift'] = 0
    segyprep['y origin'] = proj.box[2]
    segyprep['y shift'] = 0
    segyprep['z type'] = 'elevation' 
    
    #try:
    #  segyprep['x origin'] = kwargs['x1']
    #  segyprep['x shift'] = 0
    #  segyprep['y origin'] = kwargs['y1']
    #  segyprep['y shift'] = 0
    #  segyprep['z type'] = 'elevation' 
    #  
    #except KeyError:
    #  raise IOError('For sgy-geometry you have to specify x1, y1 (bottom-left corner of model box) in metres.')
  
  elif (len(Ss) == 0) and (len(Rs) == 0):
    segyprep['geometry'] = 'regular'
    
    segyprep['rec depth'] = max(int(nx3 / 2), fs_z_max + 10) * dx # m
    segyprep['sou depth'] = int(fs_z_max + 5) * dx # m
    
    segyprep['rec dx'] = 1 * dx # m
    rec_xpad = 4
    segyprep['rec nx'] = nx1 - 2 * rec_xpad
    segyprep['rec x origin'] = rec_xpad * dx # m  
    
    segyprep['sou dx'] = 0 * dx # m
    segyprep['sou nx'] = 1
    segyprep['sou x origin'] = int(nx1 / 2) * dx # m 
    
    segyprep['rec dy'] = 0 * dx # m
    segyprep['rec ny'] = 0
    segyprep['rec y origin'] = int(nx2 / 2) * dx # m
    
    segyprep['sou dy'] = 0 * dx # m
    segyprep['sou ny'] = 0
    segyprep['sou y origin'] = int(nx2 / 2) * dx # m
  
  elif (len(Ss) == 1) and (len(Rs) == 1):
    segyprep['geometry'] = 'regular' # NOTE: OTHERWISE SegyPrep WILL SET GEOMETRY TO segy!
    
    if len(Ss) != 1 or len(Rs) != 1:
      eprint(this_func + 'Error. SR lists should be 1-long.\n')
      quit()
    
    S = Ss[0]
    R = Rs[0]
                 
    segyprep['sou x origin'] = (S[0] - 1) * dx # NOTE: DOUBLE-CHECK THIS '-1' (DEPTH DOESN'T NEED IT)
    segyprep['sou y origin'] = (S[1] - 1) * dx
    segyprep['sou depth'] = S[2] * dx
    segyprep['sou nx'] = 1
    segyprep['sou ny'] = 0
      
    segyprep['rec x origin'] = (R[0] - 1) * dx # NOTE: DOUBLE-CHECK THIS '-1' (DEPTH DOESN'T NEED IT)
    segyprep['rec y origin'] = (R[1] - 1) * dx
    segyprep['rec depth'] = R[2] * dx
    segyprep['rec nx'] = 1
    segyprep['rec ny'] = 0
    
  else:
    raise ValueError('Lists of SR longer than > 1 require Prepare_SR AFTER running SegyPrep. Not yet implemented')
    
    #segyprep['geometry'] = 'regular' # NOTE: OTHERWISE SegyPrep WILL SET GEOMETRY TO segy!
    
    #S = Ss[0]
    #R = Rs[0]
    #             
    #segyprep['sou x origin'] = (S[0] - 1) * dx 
    #segyprep['sou y origin'] = (S[1] - 1) * dx
    #segyprep['sou depth'] = S[2] * dx
    #segyprep['sou nx'] = len(Ss) # NOTE: THIS IS CRUCIAL FOR CORRECT Observed-0000.vtr FILE
    #segyprep['sou ny'] = 0       ##
    #  
    #segyprep['rec x origin'] = (R[0] - 1) * dx
    #segyprep['rec y origin'] = (R[1] - 1) * dx
    #segyprep['rec depth'] = R[2] * dx
    #segyprep['rec nx'] = len(Rs) # NOTE: THIS IS CRUCIAL FOR CORRECT Observed-0000.vtr FILE
    #segyprep['rec ny'] = 0       ##
  
  segyprep['fixed array'] = 'yes' # CHANNEL NUMBERS <-> UNIQUE PHYSICAL LOCATIONS
  segyprep['unique'] = 'yes'
  segyprep['FFID'] = 'yes' # CHECK IT!
  
  outseis = Kwarg('outseis', 1, kwargs)
  if outseis:
    segyprep['outseis'] = 'yes'
  else:
    segyprep['outseis'] = 'no'
    
  segyprep['text'] = 'yes'
  segyprep['debug'] = 'yes'
  segyprep['retain'] = 'yes'
  
  
  Save_Dict(fname, segyprep)  
  
  if verbos > 4:
    command = 'cat ' + fname
    o, e = Bash2(command)
    print(this_func, 'Output of command ' + command + ': ')
    print(o, e)
  
  #print this_func + 'END' 
  return segyprep


# -------------------------------------------------------------------------------


#def Segyprep_Create_From_Template


# -------------------------------------------------------------------------------


def Segyprep_Run_n_Shift_Z(proj_name, z_sea, dx, **kwargs):
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
  It assumes it's called in the right directory
  (can't use cd).
  
  FIXME: CATCH WARNINGS OUTPUT BY SEGYPREP
  
  """
  this_func = this_lib + 'Segyprep_Run_n_Shift_Z: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_generic import Bash
  
  proj_path = Kwarg('proj_path', './', kwargs)
  #overwrite = Kwarg('overwrite', 'no', kwargs)
 
  o = Bash('fwi_run_segyprep.sh ' + proj_name, **kwargs)
  #o, e = Bash2('printf "' + overwrite + r'\n"' + " | segyprep_v3.16 " + proj_name) # NOTE  r'\n"' PRINTS \n NOT NEWLINE
  if verbos > 5:
    print(o)
    o, e = Bash2('kate ' + proj_name + '-SegyPrep.log')
  
  try:
    # SHIFT VERTICALLY
    SR_Modify_Shift_Z_geo(proj_name, z_sea, dx, proj_path=proj_path)  
    
    # SHOW THE CONTENTS
    for suffix in ['-Sources.geo', '-Receivers.geo']:
      fname = proj_name + suffix
      if verbos > 3:
        o = Bash('cat ' + fname)
        print(this_func, 'Contents of', fname, '\n', o) 
  except:
    eprint(this_func + 'Could not shift .geo files. Probably due to fw3d  format.\n')
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
# FULLWAVE INPUT (EXCLUDING EXTERNAL LIBS)
# -------------------------------------------------------------------------------


# SOURCES & RECEIVERS -----------------------------------------------------------


def SR_Prepare(proj_name, **kwargs):
  """
  Prepare sources and receivers
  input files.
  
  
  Parameters
  ----------
  proj_name : str
    Project name assumed to be provided 
    as a first argument of the FWI code
    being a prefix of all FWI input files  
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  kwargs : dict 
    Updated of Ss and Rs keywords.
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'SR_Prepare: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  pars_gen = Kwarg('pars_gen', {}, kwargs)
  io = Kwarg('io', 'fw3d', pars_gen)
  dims = Kwarg('dims', dims_default, pars_gen)
  nx1, nx2, nx3 = dims
  pars_SR = Kwarg('pars_SR', {}, kwargs)
  Ss = Kwarg('Ss', None, pars_SR)
  Rs = Kwarg('Rs', None, pars_SR)
  reciprocity = Kwarg('reciprocity', False, pars_SR)
  
  if io == 'sgy':
    # IT WILL BE TAKEN FROM sgy HEADERS
    pass    
  
  else:
    if not (Ss and Rs):
      mid = [nx1 // 2, nx2 // 2, nx3 // 2]
      off = max(5, mid[0] // 5)
    
      # CASE 1
      s = [mid[0] - off, mid[1], mid[2] - off] 
      # CASE 2 
      r = [mid[0], mid[1], mid[2] + off] 
      
      Ss = [s]
      Rs = [r]
    
    if reciprocity:
      Ss, Rs = Rs, Ss
    
    SR_Create(proj_name, Ss=Ss, Rs=Rs, **kwargs)
    
    for Xs, key in zip([Ss, Rs], ['Ss', 'Rs']):
      if not key in pars_SR:
        pars_SR[key] = Xs
    kwargs['pars_SR'] = pars_SR
    
  #print this_func + 'END'
  return kwargs


# -------------------------------------------------------------------------------


def SR_Create(proj_name, **kwargs):
  """
  Create sources and receivers
  input files.
  
  
  Parameters
  ----------
  proj_name : str
    Project name assumed to be provided 
    as a first argument of the FWI code
    being a prefix of all FWI input files  
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'SR_Create: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  # READ INPUT
  Ss = Kwarg('Ss', None, kwargs)
  Rs = Kwarg('Rs', None, kwargs)
  dims = Kwarg('dims', None, kwargs)
  proj_path = Kwarg('proj_path', './', kwargs)
  fmt = Kwarg('fmt', 'pgy', kwargs) 
  
  if not (Ss and Rs):
    eprint(this_func + 'Error. You must provide both SR lists\n')
    quit()
    
  if fmt == 'pgy':
    SR_Create_pgy(proj_name, **kwargs)
  
  elif fmt == 'sofi':
    from . import lib_io_sofi
    lib_io_sofi.SR_Create(proj_name, **kwargs)
    
  else:
    raise ValueError('Wrong value of fmt: ' + fmt)
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def SR_Create_pgy(proj_name, **kwargs):
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
  this_func = this_lib + 'SR_Create_pgy: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  pars_gen = Kwarg('pars_gen', {}, kwargs)
  dims = Kwarg('dims', dims_default, pars_gen)
  Ss = Kwarg('Ss', None, kwargs)
  Rs = Kwarg('Rs', None, kwargs)
  
  proj_path = Kwarg('proj_path', './', kwargs)
  
  if not dims:
    try:
      n1, n2, n3, model = Read_vtr(proj_path + proj_name + '-TrueVp.vtr', **kwargs)
      #print this_func, n1, n2, n3
      dims = [n1, n2, n3]
    except:
      eprint(this_func + 'Error. Dims are not provided and cannot be read from TrueVp file either.\n')
      quit()
  
  for Xs, suffix in zip([Ss, Rs], ['-PointSources.pgy', '-PointReceivers.pgy']):
    fname = proj_path + proj_name + suffix
    if 'dims' in kwargs:
      del kwargs['dims']
    Save_pgy(fname, Xs, dims, **kwargs)
  
  #print this_func + 'END'
  return 0


# MODIFY ------------------------------------------------------------------------


def SR_Modify_Adapt_To_Extras(proj_name, sources, receivers, **kwargs):
  """
  FIXME: WORKS ONLY FOR .geo
  
  Adjust sources & receivers coordinates 
  to conform to extended grid 
  (by extra nodes, see the runfile).
  
  Parameters
  ----------
  proj_name : str
    Project name assumed to be provided 
    as a first argument of the FWI code
    being a prefix of all FWI input files
  sources : list 
    List of points [x, y, z]
  receivers : list 
    List of points [x, y, z]
  
  Returns
  -------
  sources : list 
    List of points [x, y, z]
  receivers : list 
    List of points [x, y, z]
  
  Notes
  -----
  It will try to read .pgy files in the first place.
  If failed, it will try to read .geo files.
  
  """
  this_func = this_lib + 'SR_Modify_Adapt_To_Extras: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')     
  
  io = Kwarg('io', 'geo', kwargs)
  
  if io == 'fw3d':
    raise 'IO '
  
  eprint(this_func + 'Adjusting sources & receivers.\n')
  
  etop, eleft, efront, ebot, eright, eback = Set_Extra_Nodes(proj_name + '-Runfile.key')
  
  for dic in [sources, receivers]:
    for key in dic:
      x, y, z = dic[key]
      dic[key] = [x+eleft, y+efront, z+etop]
  
  #xs, ys, zs = Split_Tuples(sources)
  #xr, yr, zr = Split_Tuples(receivers)
  
  #xs = [i+eleft for i in xs]
  #xr = [i+eleft for i in xr]
  
  #ys = [i+efront for i in ys]
  #yr = [i+efront for i in yr]
  
  #zs = [i+etop for i in zs]
  #zr = [i+etop for i in zr] 
  
  #sources = zip(xs, ys, zs)
  #receivers = zip(xr, yr, zr)
  
  #print this_func + 'END'  
  return sources, receivers


# -------------------------------------------------------------------------------


#def SR_Modify_Shift_Z(proj_name, shift_z, dx, **kwargs):


# -------------------------------------------------------------------------------


#def SR_Modify_Shift_Z_pgy(proj_name, shift_z, dx, **kwargs):


# -------------------------------------------------------------------------------


def SR_Modify_Shift_Z_geo(proj_name, shift_z, dx, **kwargs):
  """
  This necessary in 'segy' geometry because SR positions 
  from SGY headers refer to sea surface which is not 
  the position of the free surface. And that's what is 
  referred to by fullwave.
  
  Parameters
  ----------
  proj_name : str
    Project name assumed to be provided 
    as a first argument of the FWI code
    being a prefix of all FWI input files
  shift_z : float 
    No. of nodes to add to original z-coordinate of S/R.
  dx : float 
    Size of the grid cell (m).
  
  **kwargs 
  
  Returns
  -------
  0
  
  sources : list 
    List of points [x, y, z]
  receivers : list 
    List of points [x, y, z]
  
  Notes
  -----
  NOTE: IT NEEDS SEGYPREP RE-RUN EVERYTIME (.geo ARE MODIFIED)
  
  """
  this_func = this_lib + 'SR_Modify_Shift_Z: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  proj_path = Kwarg('proj_path', './', kwargs)
  
  for suffix in ['-Sources', '-Receivers']:
    ext = '.geo'
    fname = proj_path + proj_name + suffix + ext 
    interfix = '_not_corrected'
    fname_nc = proj_path + proj_name + suffix + interfix + ext
    command = 'cp ' + fname + ' ' + fname_nc
    o = Bash(command)
  
    c = Read_File(fname_nc)
    header = c[0]
    header_str = ''
    for word in header:
      header_str += word + ' ' # FIXME: BETTER FORMATTING
    data = c[1: ]
    f = open(fname, 'w')
    f.write(header_str + '\n')
  
    for line in data:
      source_id = line[0]
      x = line[1]
      y = line[2]
      z = line[3]    
      z = str(shift_z * dx + float(z))
      f.write(source_id + ' ' + x  + ' ' + y  + ' ' + z  + '\n')
    f.close()
  
  #print this_func + 'END'
  return 0


# CHECK --------------------------------------------------------------------------


def SR_Check_Model_Bounds(proj_name, adjust=True, **kwargs):
  """
  Check sources & receivers positions.
  
  Parameters
  ----------
  proj_name : str
    Project name assumed to be provided 
    as a first argument of the FWI code
    being a prefix of all FWI input files
  adjust : bool 
    If true adjust SR coordinates to conform 
    to extended grid (extra nodes, see the runfile)
  **kwargs : keyword arguments (optional)
    Current capabilities:
    proj_path : str 
      Path to the project.
      Default: './'
    fs_file : str 
      File storing gridded FS z-coordinates.
      Default: 
      proj_path + proj_name + -FreeSurf.vtr (if adjust=False)
      proj_path + proj_name + -FreeSurf_exten.vtr (if adjust=True)
    margin_z : float
      
      Default: 0.
    
    
  Returns
  -------
  0
  
  Notes
  -----
  
  """
  this_func = this_lib + 'SR_Check_FS_Position: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')    
  
  try:
    proj_path = kwargs['proj_path']
  except KeyError:
    proj_path = "./"  
  
  try:
    fs_file = kwargs['fs_file']
  except KeyError:
    fs_file = proj_path + proj_name
    if adjust:
      fs_file += '-FreeSurf_exten.vtr'
    else:
      fs_file += '-FreeSurf.vtr'
  
  try:
    margin_z = kwargs['margin_z']
  except KeyError:
    margin_z = 0 
  
  dims, sources, receivers = SR_Read(proj_name, adjust, proj_path=proj_path)
  nx1_fs, nx2_fs, nx3_fs, fs_z = Read_vtr(fs_file)
  
  err = False
  for ptype, points in zip(['Source', 'Receiver'], [sources, receivers]):
    for i, p in enumerate(points):
      pass
      #px, py, pz = points[p]
      #fs_z_p = fs_z[int(px)-1][int(py)-1][0] # TO COMPLY WITH PYTHON INDEXING
    
      #if pz <= fs_z_p + margin_z:
        #if pz < fs_z_p + margin_z:
          #eprint(this_func + '' + ptype + ' no. ' + str(i + 1) + ' above FS + margin ' + str(margin_z) + '.\n')
        #elif pz == fs_z_p + margin_z:
          #eprint(this_func + '' + ptype + ' no. ' + str(i + 1) + ' exactly on FS + margin ' + str(margin_z) + '.\n')
        
        #eprint(this_func + ptype + ': ' + str(p) + '\n')
        #eprint(this_func + 'Z-coord of the free surface at this point: ' + str(fs_z_p) + '\n')
        #err = True
        
  if err:
    raise ValueError('Positions of sources/receivers are not correct.')
  
  print(this_func, 'Check successful - all sources and receivers are below the free surface.')
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def SR_Check_FS_Position(proj_name, adjust=True, **kwargs):
  """
  Check sources & receivers positions.
  
  Parameters
  ----------
  proj_name : str
    Project name assumed to be provided 
    as a first argument of the FWI code
    being a prefix of all FWI input files
  adjust : bool 
    If true adjust SR coordinates to conform 
    to extended grid (extra nodes, see the runfile)
  **kwargs : keyword arguments (optional)
    Current capabilities:
    proj_path : str 
      Path to the project.
      Default: './'
    fs_file : str 
      File storing gridded FS z-coordinates.
      Default: 
      proj_path + proj_name + -FreeSurf.vtr (if adjust=False)
      proj_path + proj_name + -FreeSurf_exten.vtr (if adjust=True)
    margin_z : float
      
      Default: 0.
    
    
  Returns
  -------
  0
  
  Notes
  -----
  In nodes.
  
  
  """
  this_func = this_lib + 'SR_Check_FS_Position: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')    
  
  proj_path = Kwarg('proj_path', './', kwargs)
  margin_z = Kwarg('margin_z', 0, kwargs)  
  try:
    fs_file = kwargs['fs_file']
  except KeyError:
    fs_file = proj_path + proj_name
    if adjust:
      fs_file += '-FreeSurf_exten.vtr' # NOTE
    else:
      fs_file += '-FreeSurf.vtr'
  
  
  dims, sources, receivers = SR_Read(proj_name, adjust, proj_path=proj_path)
  nx1_fs, nx2_fs, nx3_fs, fs_z = Read_vtr(fs_file)
  
  err = False
  for ptype, points in zip(['Source', 'Receiver'], [sources, receivers]):
    for i, p in enumerate(points):
      px, py, pz = points[p]
      fs_z_p = fs_z[int(px)-1][int(py)-1][0] # TO COMPLY WITH PYTHON INDEXING
    
      if pz <= fs_z_p + margin_z:
        if pz < fs_z_p + margin_z:
          eprint(this_func + '' + ptype + ' no. ' + str(i + 1) + ' above (FS + margin ' + str(margin_z) + ').\n')
        elif pz == fs_z_p + margin_z:
          eprint(this_func + '' + ptype + ' no. ' + str(i + 1) + ' exactly on (FS + margin ' + str(margin_z) + ').\n')
        
        eprint(this_func + ptype + ' ID: ' + str(p) + '\n')
        eprint(this_func + 'Z-coord of the ' + ptype + ' ' + str(pz) + '\n')
        eprint(this_func + 'Z-coord of the free surface at this point: ' + str(fs_z_p) + '\n')
        
        err = True
  
  # THIS IS DONE AFTER THE FULL LOOP TO DISPLAY ALL FAULTY  CASES
  if err:
    raise ValueError('Positions of sources/receivers are not correct.')
  
  print(this_func, 'Check successful - all sources and receivers are below the free surface.')  
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def SR_Check_Sparsity(proj_name, dx, **kwargs):
  """
  
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  spacings : list 
    Of srcs and rcvrs respectively.
  
  Notes
  -----
  Be mindful of reciprocity.
  
  """
  this_func = this_lib + 'SR_Check_Sparsity: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from numpy.linalg import norm
  
  dims, srcs, rcvrs = SR_Read(proj_name, dx=dx, **kwargs)
  
  spacings = []
  
  for SRs, name in zip([srcs, rcvrs], ['sources', 'receivers']):
  
    min_dists = []
    for key1 in SRs:
      xy1 = np.array(SRs[key1][:2])
      min_dist = 1e5
      for key2 in SRs:
        if key1 == key2:
          continue
        xy2 = np.array(SRs[key2][:2])
        dist = norm(xy1 - xy2)
        if dist < min_dist:
          min_dist = dist
      min_dists.append(min_dist)
      #print 'nearest OBS for OBS ', key1, ' is: ', min_dist * 50 / 1000., ' km away'
      
    spacing_km =  (sum(min_dists) / len(min_dists)) * dx / 1000.
    print('average ', name, ' spacing is ',  spacing_km, ' km')
    spacings.append(spacing_km)

  #print this_func + 'END'
  return spacings


# -------------------------------------------------------------------------------


def SR_Read(proj_name, adjust=True, **kwargs):
  """
  Read sources & receivers position 
  from Fullwave's geometry files.
  
  Parameters
  ----------
  proj_name : str
    Project name assumed to be provided 
    as a first argument of the FWI code
    being a prefix of all FWI input files
  adjust : bool 
    If true adjust SR coordinates to conform 
    to extended grid (extra nodes, see the runfile)
  **kwargs : keyword arguments (optional)
    Current capabilities:
    proj_path : str 
      Path to the project.
      Default: './'
    dx : float 
      Size of the grid cell.
    
  Returns
  -------
  nx1, nx2, nx3 : int
    Dimension of the grid as in the geom.
    files.
  sources : list 
    List of points [x, y, z]
  receivers : list 
    List of points [x, y, z]
  
  Notes
  -----
  It will try to read .pgy files in the first place.
  If failed, it will try to read .geo files.
  
  """
  this_func = this_lib + 'SR_Read: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')    
  
  proj_path = Kwarg('proj_path', './', kwargs)
  prefix = proj_path + proj_name
  io = Kwarg('io', 'sgy', kwargs)
  
  ## CHECK EXISTENCE
  #if Exists(prefix + '-PointSources.pgy'):
  #  ext = 'pgy'
  #elif Exists(prefix + '-Sources.geo'):
  #  ext = 'geo'
  #else:
  #  raise IOError('Neither .pgy nor .geo source file has been found.')
  
  
  # READ
  if io == 'fw3d':
    nx1, nx2, nx3, sources = Read_pgy(prefix + '-PointSources.pgy', **kwargs) 
    nx1_r, nx2_r, nx3_r, receivers = Read_pgy(prefix + '-PointReceivers.pgy', **kwargs)
    # REDUNDANT CHECK?
    if nx1 == nx1_r and nx2 == nx2_r and nx3 == nx3_r:
      dims = [nx1, nx2, nx3]
    else:
      raise ValueError('Dims in SR files differ.')
    
  elif io == 'sgy':
    # REDUNDANT CHECK?
    try:
      dx = kwargs['dx']
    except KeyError:
      raise IOError('You need to provide dx to read the .geo files.')
    
    del kwargs['dx'] # SUPRISING ERROR: KEYWORD CAN'T ANY OF REGULAR ARGUMENTS!
    dims, sources = Read_geo(prefix + '-Sources.geo', dx, **kwargs) 
    dims, receivers = Read_geo(prefix + '-Receivers.geo', dx, **kwargs)   
  
  else:
    raise ValueError('Unknown io: ' + io)
  
  # CATCH INCONSISTENCY
  #if nx1 != nx1_r or nx2 != nx2_r or nx3 != nx3_r:
  #  eprint(this_func + 'Error. Inconsistent dimensions between S. and R. file \n')
  #  quit()
    
  if adjust:
    sources, receivers = SR_Modify_Adapt_To_Extras(prefix, sources, receivers)
  
  return dims, sources, receivers


# -------------------------------------------------------------------------------


def SR_Read_All_PROTEUS(**kwargs):
  """
  Read true positions in PROTEUS-experiment
  frame of local coordinates of all 
  OBSes and shots.
  
  Parameters
  ----------
    
  Returns
  -------
  
  Notes
  -----
  
  """
  this_func = this_lib + 'SR_Read_All_PROTEUS: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')    
  
  path = '/home/kmc3817/heavy_PhD/meta_data/'
  
  dx = 1 # WILL PRESERVE METRES (INSTEAD OF CONVERTING INTO NODES)
  if 'dx' in kwargs:
    del kwargs['dx']
  dims, sources = Read_geo(path + 'all-Sources.true', dx, **kwargs) 
  dims, receivers = Read_geo(path + 'all-Receivers.true', dx, **kwargs)   
  
  return sources, receivers



# RUNFILE -----------------------------------------------------------------------


def Runfile_Create(proj, **kwargs):
  """
  Prepare a dictionary of Runfile arguments.
  
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
  # FIXME: merge with make_segyprep
  
  """
  this_func = this_lib + 'Runfile_Create: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  #from lib_generic import Save_Dict

  create_from = Kwarg('create_from', 'template', kwargs)
  
  if create_from == 'scratch':
    Runfile_Create_From_Scratch(proj, **kwargs)
    
  elif create_from == 'template':
    Runfile_Create_From_Template(proj, **kwargs)
    
  else:
    raise ValueError('create_from: ' + create_from)
  
  #print this_func + 'END' 
  return 0


# -------------------------------------------------------------------------------


def Runfile_Create_From_Scratch(proj, **kwargs):
  """
  FIXME: NOT READY YET
  
  Parameters
  ----------
  
  Returns
  -------
  0
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Runfile_Create_From_Scratch: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
   
  #ball = 1 # FIXME
  #eall = 10
  
  
  #print "path", proj_path
  #print "proj_name", proj_name
  #print "file_name", file_name
  #print "ball", ball
  #print "eall", eall
  
  #runfile = {}
  #file_name = proj_name + "-Runfile.key"
  
  #runfile["problem"] = "synthetic" # "synthetic", "tomography", "both
  #runfile["equation"] = "acoustic"
  #runfile["anisotropy"] = "none"
  #runfile["domain"] = "time"
  #runfile["kernel"] = "low"
  #runfile["io"] = "fw3d" 
  
  #runfile["etime"] = 2
  
  #runfile["nx1"] = 200
  #runfile["nx2"] = 1
  #runfile["nx3"] = 200
  #runfile["dx"] = 50
  
  #runfile["ncomp"] = 1
  #runfile["nshots"] = 1
  #runfile["nrecs"] = 101
  #runfile["maxrc"] = 101
  
  #runfile["btop"] = 0
  #runfile["bbot"] = ball
  #runfile["bleft"] = ball
  #runfile["bright"] = ball
  #runfile["bfront"] = ball
  #runfile["bback"] = ball
  
  #runfile["extratop"] = 0
  #runfile["extrabot"] = eall
  #runfile["extraleft"] = eall
  #runfile["extraright"] = eall
  #runfile["extrafront"] = eall
  #runfile["extraback"] = eall 
  
  # f = open(proj_path + file_name, "w")
  
  #f.write
  
  
  
  # for i in runfile:
  #     f.write(i + " : " + str(runfile[i]) + "\n")
  # f.close()
  
  # !cd "{proj_path}" && cat {file_name}  

  #Save_Dict(file_name, runfile, **kwargs)

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Runfile_Create_From_Template(proj, **kwargs):
  """
  Prepare the project runfile 
  by copying the template (my format) and 
  updating it with the skeleton generated
  by SegyPrep.
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
  
  Returns
  -------
  0
  
  Notes
  -----
  It assumes Skeleton exists!
  
  """  
  this_func = this_lib + 'Runfile_Create_From_Template: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')    
  
  proj_name = proj.name
  proj_path = proj.path + '/inp/'
  z_sea =  proj.z_sea  

  ibfs = Kwarg('ibfs', 0, kwargs)
  topabs = Kwarg('topabs', False, kwargs)
  b_abs = Kwarg('b_abs', 40, kwargs)
  e_abs = Kwarg('e_abs', 50, kwargs)


  # GET ALL THE FILES
  machine = Kwarg('machine', proj.machine, kwargs)
  paths = Paths_Prepare(machine)
  template = paths['templates'] + paths['template_runfile']    
  skeleton = proj_path + proj_name + '-Skeleton.key'
  runfile = proj_path + proj_name + '-Runfile.key'
  #runfile_synth = proj_path + proj_name + '-Runfile_synth.key'
  #runfile_inv = proj_path + proj_name + '-Runfile_inv.key'
  
  # COPY TEMPLATE INTO A NEW RUNFILE
  o, e = Bash2('cp ' + template + ' ' + runfile)
  #print template
  
  # NOTE: READ PROJECT-SPECIFIC VALUES FROM SKELETON GENERATED BY SEGYPREP
  sgprep_skelet = Runfile_Read(skeleton) # NOTE IS IT NEEDED ANYMORE?
  
  #print sgprep_skelet #['problem']
  
  # UPDATE THE TEMPLATE RUNFILE BASED ON THE SKELETON 
  #FIXME pROBLEM
  Runfile_Update(runfile, problem=proj.problem, io=sgprep_skelet['IO'], etime=sgprep_skelet['ttime'],
                 bbot=str(b_abs), bleft=str(b_abs), bright=str(b_abs), bfront=str(b_abs), bback=str(b_abs), 
                 extrabot=str(e_abs), extraleft=str(e_abs), extraright=str(e_abs), extrafront=str(e_abs), extraback=str(e_abs), 
                 ibfs=str(ibfs), seaLevel=str(z_sea),  
                 NX1=sgprep_skelet['NX1'], NX2=sgprep_skelet['NX2'], NX3=sgprep_skelet['NX3'], 
                 NRECS=sgprep_skelet['NRECS'], MAXRC=sgprep_skelet['MAXRC'], 
                 NCOMP=sgprep_skelet['NCOMP'], NSHOTS=sgprep_skelet['NSHOT'])
  
  if topabs:
    Runfile_Update(runfile, btop=str(b_abs), extratop=str(e_abs))
  
  #!cd "{proj_path}" && cat {proj_name}'-Runfile.key'
  
  Runfile_Update(runfile, **kwargs) # FIXME: DOUBLE-CHECK
  
  # SPLIT INTO SYNTH/INV FLOWS
  #o = Bash('cp ' + runfile + ' ' + runfile_synth)
  #o = Bash('cp ' + runfile + ' ' + runfile_inv)
  
  #o, e = Bash2('cat ' + runfile_synth)
  #print 'oo', o
  #print runfile_inv
  #Runfile_Update(runfile_inv, problem='tomography')
  #o, e = Bash2('cat ' + runfile_synth)
  #print 'oo00', o  
  
  #o = Bash('sed -i -e "s/synthetic/tomography/g" ' + runfile_inv) # THIS DOESN'T WORK
  
  if verbos > 1:
    o, e = Bash2('cat ' + runfile)
    print(this_func, 'Contents of ' + runfile + ':\n', o)
    #o = Bash('cat ' + runfile_inv)
    #print this_func, 'Contents of ' + runfile_inv + ':\n', o
  
  #print this_func + 'END'
  return 0


# ------------------------------------------------------------------------------


def Runfile_Read(file_name, **kwargs):
  """
  READ Runfile.key AND RETURN RECORDS OF FORMAT: 
  '1-word-key : value' 
  AS A DICTIONARY 

  Parameters
  ----------
  file_name : str 
    Full name of the file including 
    a path if needed.  
  **kwargs : keyword arguments (optional)
  
  Returns
  -------
  0
  
  Notes
  -----
  Reads only 1-word parameters!!!
  
  It should work for SegyPrep.key as well!
  
  """  
  this_func = this_lib + 'Runfile_Read: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')    
  
  
  
  if verbos > 3:
    print(this_func, 'Runfile to read: ', file_name)
  
  content = Read_File(file_name)
  params = {}
  for record in content:
    # GET ONLY RECORDS OF CERTAIN FORMAT:
    if (len(record) == 3) and (record[1] == ':'):
      key = record[0]
      value = record[2]
      params[key] = value
      
    # SPECIAL FORMAT
    if (record[0] == 'MAX') and (record[1] == 'TIME'):
      key = 'ttime'
      value = str(float(record[3]) / 1000.) # CONVERT ms TO s
      params[key] = value
  
  #print this_func + 'END'  
  return params


# -------------------------------------------------------------------------------


def Runfile_Update(file_name, **kwargs): 
  """
  CHANGE RECORDS OF RUNFILE, THE REST REMAINS EXACTLY THE SAME
   
  Parameters
  ----------
  file_name : str 
    Full name of the file including 
    a path if needed.  
  **kwargs : keyword arguments (optional)
  
  Returns
  -------
  0
  
  Notes
  -----
  
  """  
  this_func = this_lib + 'Runfile_Update: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  content = Read_File_Not_Split(file_name)
  
  f = open(file_name, 'w')
  
  for line in content:
    nline = line
    line = line.split(None)
    if len(line) == 0:
      continue
    
    for key in kwargs:
      value = kwargs[key]
      if key == line[0]:
        if value == 'disable':
          nline = '!' + key  + '\n' # FIXME: ADD FULL LINE FOR AESTHETICS
        elif value == 'enable':
          raise NotImplementedError('value=enable')
        else:
          nline = key + ' : ' + str(value) + '\n'
        break 
  
    f.write(nline)
  f.close()
  
  #print this_func + 'END'
  return 0  


# -------------------------------------------------------------------------------


def Runfile_Set_Iterations(file_name, blocks, **kwargs): 
  """
  Set section E of the runfile.
  
  Parameters
  ----------
  file_name : str 
    Full name of the file including 
    a path if needed.  
  blocks : list 
    List of Iteration_Block instances.
  
  **kwargs : keyword arguments (optional)
  
  Returns
  -------
  0
  
  Notes
  -----
  
  """  
  this_func = this_lib + 'Runfile_Set_Iterations: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  # READ THE EXISITING RUNFILE
  content = Read_File_Not_Split(file_name)
  
  f = open(file_name, 'w')
  # REWRITE ALL SECTION BEFORE SECTION E (ITERATION BLOCKS)
  for line in content:
    nline = line
    line = line.split(None)
    if len(line) > 0 and line[0] == 'NBLOCK':
      break
    else:
      f.write(nline) # NOTE
  
  # ADD CUSTOM SECTION E
  f.write('NBLOCK : ' + str(len(blocks)) + '\n')
  
  for b in blocks:
    f.write('Frequency : ' + str(b.freq) + '\n')
    f.write('Iterations : ' + str(b.niters) + '\n')
    if b.minoff:
      f.write('minoff : ' + str(b.minoff) + '\n')
    if b.maxoff:
      f.write('maxoff : ' + str(b.maxoff) + '\n')      
    
  f.close()
  
  #print this_func + 'END'
  return 0  


# PBS SCRIPT --------------------------------------------------------------------


def PBS_File_Create(proj, **kwargs):
  """
  Create a PBS script file for 
  submitting a job to a queueing system. 
  
  Parameters
  ----------
  proj_name : str
    Project name assumed to be provided 
    as a first argument of the FWI code
    being a prefix of all FWI input files
    
  **kwargs : keyword arguments (optional)
    Current capabilities:
    path : str 
      Path to the project.
      Default: './'
    problem : str
      Fullwave3D problem
      - 'synth' for synthetic
      - 'inv' for inversion/tomography
      - 'both' for both
      Default: 'synth'.
  
  Returns
  -------
  0
  
  Notes
  -----
  NOTE:(source: Imperial MPI jobs) 
  as soon as the walltime is reached the content
  of $TMPDIR is DELETED.
  
  pbsexec before mpiexec gives 20 min grace period,
  i.e.it terminates the code run with mpiexec 20 min 
  before end of the requested wall time. This allows 
  the rest of the .pbs script to be executed. 
  
  pbsdsh2 is necessary for throughputting on ALL nodes 
  of a multi-node job.
  
  Jo doesn't use it at all (cd to dir). But copying was
  taught in HPC training.
  
  If * follows a prefix, it must be proceeded by \, e.g.:
  $work_dir/prefix\*.vtr
  
  """
  this_func = this_lib + 'PBS_File_Create: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  proj_name = proj.name
  path = proj.inp.path
  problem = proj.problem
  if problem == 'syn' or problem == 'synth':
    problem = 'synthetic'
  if problem == 'inv' or problem == 'inversion':
    problem = 'tomography'

  fname = path + proj_name + '-Run.pbs'
  f = open(fname, 'w')
  f = PBS_File_Write_Header(f, proj_name, proj.env.var, **kwargs) # NOTE!!!
  if problem == 'synthetic' or problem == 'both':
    f.write('\n') 
    f.write('#\n') 
    f.write('# SYNTHETIC-RUN FLOW\n') 
    f.write('#\n') 
    f.write('\n') 
    
    f.write('# DELETE OBSERVED DATA, OTHERWISE SYNTH WOULD FAIL\n')
    f.write('rm $this_dir/${project_name}-*.log\n')
    f.write('rm $this_dir/${project_name}-*fw*vtr\n')
    f.write('rm $this_dir/${project_name}-Observed-Time.tt?\n')
    f.write('rm $this_dir/${project_name}-Synthetic.*\n')
    f.write('rm $this_dir/${project_name}-Observed.*\n')

    f.write('\n')
    f.write('# COPY ALL INPUT TO WHERE THE COMPUTATION WILL BE PERFORMED\n')
    f.write('pbsdsh2 cp $this_dir/../inp/${project_name}-* $work_dir\n')
    f.write('\n')
    f.write('# RUN THE CODE\n')
    f.write('pbsexec mpiexec ${code_path} ${project_name} ${CPNUM} ${NTHREAD} ')
    f.write('1> ${project_name}-SynOut.log ' )
    f.write('2> ${project_name}-SynErr.log \n' )
    f.write('stat=$?\n')
    f.write('echo "Exit status: "$stat\n') 
    f.write('\n')
    f.write('# COPY RESULTS (BOTH SEGY AND FW3D IF PRESENT)\n')
    f.write('pbsdsh2 cp $work_dir/${project_name}-*.log $this_dir/../out\n')
    f.write('pbsdsh2 cp $work_dir/${project_name}-Observed-Time.ttr $this_dir/../out\n')
    f.write('pbsdsh2 cp $work_dir/${project_name}-Observed-Time.ttm $this_dir/../out\n')
    f.write('pbsdsh2 cp $work_dir/${project_name}-Synthetic.sgy $this_dir/../out\n')
    f.write('pbsdsh2 cp $work_dir/${project_name}-Synthetic.idx $this_dir/../out\n')
    f.write('pbsdsh2 cp $work_dir/*fw\*.vtr $this_dir/../out # ALL (SNAPSHOTS!)\n') # prefix\* IS NEEDED
    f.write('cd $this_dir/../out/\n')
    f.write('for f in *fw* \n')
    f.write('do \n')
    f.write('  mv $f  ${project_name}-$f\n') # (ANY VERSION OF) RENAME DOESN'T WORK ON CX1
    f.write('done \n')
    f.write('cd $this_dir\n')
  if problem == 'tomography' or problem == 'both':
    f.write('\n') 
    f.write('#\n') 
    f.write('# INVERSION-RUN FLOW\n') 
    f.write('#\n') 
    f.write('\n')     
    #f.write('cd $this_dir # MAKE SURE WE ARE IN THE RIGHT DIRECTORY\n')
    #f.write('\n') 
    f.write('unset SLAVES_WAVEFIELDSVTR # OTHERWISE TOO MANY SNAPSHOTS \n')
    f.write('\n') 
    #f.write('# CHOOSE THE APPROPRIATE RUNFILE FOR INVERSION \n')
    #f.write('cp $this_dir/${project_name}-Runfile_inv.key $this_dir/${project_name}-Runfile.key \n')
    #f.write('\n') 
    if problem == 'both': # CHECKERBOARD
      f.write('# ---------------- CHECKERBOARD BIT (problem="both") ---------------- \n')
      f.write('\n') 
      f.write('cp $this_dir/../out/${project_name}-Observed-Time.tt? $this_dir/ \n')
      f.write('cp $this_dir/../out/${project_name}-Synthetic.* $this_dir/ \n')
      f.write('cp $this_dir/${project_name}-Synthetic.sgy $this_dir/${project_name}-Observed.sgy \n')
      f.write('# JUST IN CASE Synthetic.idx IS NOT PRESENT: \n')
      f.write('cp $this_dir/${project_name}-Template.idx $this_dir/${project_name}-Observed.idx \n')
      f.write('# IF Synthetic.idx IS PRESENT: \n')
      f.write('cp $this_dir/${project_name}-Synthetic.idx $this_dir/${project_name}-Observed.idx \n')
      f.write('\n') 
      f.write('# ---------------- END OF CHECKERBOARD BIT ---------------------------\n')
      f.write('\n') 
    f.write('# COPY ALL INPUT TO WHERE THE COMPUTATION WILL BE PERFORMED \n')
    f.write('pbsdsh2 cp $this_dir/../inp/${project_name}-* $work_dir \n')
    f.write('\n') 
    f.write('# RUN THE CODE                  \n')
    f.write('pbsexec mpiexec ${code_path} ${project_name} ${CPNUM} ${NTHREAD} ')
    f.write('1> ${project_name}-InvOut.log ' )
    f.write('2> ${project_name}-InvErr.log \n' )
    f.write('\n') 
    f.write('# COPY RESULTS  \n')
    f.write('pbsdsh2 cp $work_dir/*CP\*vtr $this_dir/../out        \n')
    f.write('pbsdsh2 cp $work_dir/*CP\*sgy $this_dir/../out        \n')
    f.write('pbsdsh2 cp $work_dir/*CP\*key $this_dir/../out        \n')
    f.write('pbsdsh2 cp $work_dir/*SCHEDULER* $this_dir/../out        \n')
    f.write('cp $work_dir/*SLAVES* $this_dir/../out        \n')
    f.write('pbsdsh2 cp $work_dir/*LastCheckpoint.txt $this_dir/../out  \n')
  f.write('pbsdsh2 cp $work_dir/*Out.log $this_dir/../out  # COPY ALL INTERNAL LOG FILES \n')
  f.write('pbsdsh2 cp $work_dir/*Err.log $this_dir/../out  # \n')
  f.write('\n')  
  f.write('\n')
  f.write('end=`date +%s` \n')
  f.write('runtime=$((end-start)) # in sec \n')
  f.write("echo 'RUNTIME OF THE SCRIPT: '$runtime' s' \n")
  f.close()                                                                                                   
  
  if verbos > 1:
    o, e = Bash2('cat ' + fname)
    print(this_func, 'Contents of ' + fname + ':\n', o)
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def PBS_File_Write_Header(f, proj_name, env_vars, **kwargs):
  """
  Create a header (resources, modules etc.) 
  PBS script file for submitting a job to a queueing system. 
  
  Parameters
  ----------
  f : file
    Already open file.
  proj_name : str
    Project name assumed to be provided 
    as a first argument of the FWI code
    being a prefix of all FWI input files
  **kwargs : keyword arguments (optional)
    
  Returns
  -------
  f : file 
    File with filled header.
  
  Notes
  -----
  Nomenclature:
  process = MPI rank = rank = task = job (SET AS mpiprocs AND FULLWAVE'S -K WHEN RUN LOCALLY)
  thread = ompthread  = virtual core (SET BOTH AS ompthreads FULLWAVE'S numthreads) 
  cpu = logic_core = logical core (SET AS ncpus)
  
  Unlike local runs,
  "it's not necessary to add "-n" or any other flag to specify the number of ranks."  
  (http://www.imperial.ac.uk/admin-services/ict/self-service/research-support/
  rcs/computing/high-throughput-computing/configuring-mpi-jobs/)
  
  "default values of  mpiprocs==ncpus and ompthreads==1" 
  (-||-)
  
  """
  this_func = this_lib + 'PBS_File_Write_Header: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  code = Kwarg('code', 'rev690', kwargs)
  code_path = Kwarg('code_path', 
                    "/rds/general/user/kmc3817/home/PhD/fullwave3D/" + 
                    code + "/bin/fullwave3D.exe", kwargs)
  cpnum = Kwarg('cpnum', -1, kwargs) # CHECKPOINT NO. TO RESTART FROM (ARG. OF fullwave3D.exe)
  job_name = Kwarg('job_name', proj_name, kwargs)

  #NOTE: NITTY-GRITTY
  time, nnodes, nprocesses, nthreads, ncpus, mem_gb = PBS_File_Set_Resources(**kwargs)
  
  # SHEBANG
  f.write('#!/bin/bash\n')
  f.write('\n')

  # INITIAL COMMENTS
  f.write('##\n')
  f.write("# This script is supposed to be run in a 'project/out/' directory\n")
  f.write("# and a directory 'project/inp/' must exist too.\n")
  f.write('#\n')
  #f.write("# A variable: 'project_name' is meant to be passed as a script argument:\n")
  #f.write('#  $> qsub -v project_name=<project_name> this_script.pbs\n')
  #f.write('#\n')
  f.write('# For more explanation, see lib_fwi_generic.py/PBS_File_Write_Header.\n')
  f.write('#\n')
  f.write('##\n')
  f.write('\n')
  f.write('#PBS -N ' + job_name + '\n')
  f.write('#PBS -o ' + proj_name + '-Out.log\n')
  f.write('#PBS -e ' + proj_name + '-Err.log\n')
  f.write('#PBS -l walltime=' + time + '\n')
  f.write('#PBS -l select=' + str(nnodes) +
          ':mpiprocs=' + str(nprocesses) + 
          ':ompthreads=' + str(nthreads) + 
          ':ncpus=' + str(ncpus) + 
          ':mem=' + mem_gb + '\n')    
  f.write('#PBS -l place=scatter:excl\n') 
  f.write('\n')  
  f.write('# PROJECT INFO\n')
  f.write('project_name=' + proj_name + '\n')
  f.write('\n')
  f.write('# PATHS\n')
  f.write('code_path=' + code_path + '\n')
  f.write("echo 'path to the code-executable: '${code_path}\n")
  f.write('this_dir=$PBS_O_WORKDIR\n')
  f.write('work_dir=$TMPDIR\n') 
  f.write('\n')
  f.write("# DISABLE PINNING OF THE PROCESSES (MAKE ALL NODE CORES AVAILABLE TO ALL PROCESSES)\n")
  f.write('unset NCPUS\n')
  f.write('export I_MPI_PIN=no # (DOES NOT WORK IF RUNNING ON A WHOLE NODE) \n')
  f.write('\n')
  f.write("# ARGUMENTS OF fullwave3D.exe\n")
  f.write('CPNUM=' + str(cpnum) + '\n')
  f.write('NTHREAD=' + str(nthreads) + ' # EQUAL TO ompthreads\n')
  f.write('\n')
  f.write("# FULLWAVE'S ENVIRONMENTAL VARIABLES\n")
  f = PBS_File_Write_Env(f, env_vars, **kwargs) #NOTE
  f.write('# LOAD MODULES\n')
  f.write('module load mpi\n')
  f.write('module load intel-suite\n')
  f.write('\n')
  #f.write('ulimit -s unlimited\n') # FIXME: it was probably me fixing some memory runtime errors
  #f.write('\n')
  f.write('start=`date +%s`\n')
  f.write('\n')
  
  #print this_func + 'END'
  return f 


# -------------------------------------------------------------------------------


def PBS_File_Set_Resources(**kwargs):
  """

  """
  this_func = this_lib + 'PBS_File_Set_Resources: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
    
  q = Kwarg('q', 'multi', kwargs) # QUEUE
  #nshots = Kwarg('nshots', 1, kwargs) 
  time = PBS_File_Set_Time(**kwargs)
  nnodes = Kwarg('nnodes', 3, kwargs) # NO. OF NODES (3 FOR multi QUEUE)
  nprocesses = Kwarg('nprocesses', 2, kwargs) # NO. OF PROCESSES PER NODE
  if nprocesses < 2:
    eprint(this_func + 'Error. At least 2 processes (1 sched + 1 slave) needed per node\n')
    #quit()  
  ncpus = Kwarg('ncpus', 2, kwargs) 
  mem_gb = Kwarg('mem_gb', 62, kwargs) # DEFAULT IS FOR GENERAL QUEUE
  
  # NOTE: IT WILL OVERWRITE SOME OF THE ABOVE VALUES
  err_nnodes = False
  if q == 'mrwarn':
    node_type = Kwarg('node_type', 'mrwarn_both', kwargs)
    mem_gb = 128
    if node_type == 'mrwarn_old' or node_type == 'mrwarn_both':
      ncpus =  40
    elif node_type == 'mrwarn_new':
      ncpus =  48
    eprint(this_func + 'Automatically setting ncpus to ' + str(ncpus) + ' and mem_gb to ' + str(mem_gb) + '\n')
    
  elif q == 'general':
    ncpus = Kwarg('ncpus', 32, kwargs)
    mem_gb = Kwarg('mem_gb', 62, kwargs) # OR 124 
    #eprint(this_func + 'Automatically setting ncpus and mem_gb\n')
    if nnodes > 16:
      err_nnodes = True
  
  elif q == 'single':
    ncpus = 48
    mem_gb = 124 
    eprint(this_func + 'Automatically setting ncpus and mem_gb\n')
    if nnodes > 1:
      err_nnodes = True
 
  elif q == 'multi':
    ncpus = 12
    mem_gb = 46
    eprint(this_func + 'Automatically setting ncpus and mem_gb\n')
    if nnodes > 16 or nnodes < 3:
      err_nnodes = True
  
  elif q == 'debug':
    if ncpus > 8:
      eprint(this_func + 'Error. Wrong ncpus.\n')
      quit()        
    mem_gb = 96
    if nnodes > 1:
      err_nnodes = True
  
  elif q == 'long':
    if ncpus > 8:
      eprint(this_func + 'Error. Wrong ncpus.\n')
      quit()        
    mem_gb = 96
    if nnodes > 1:
      err_nnodes = True
  
  else:
    eprint(this_func + 'Error. Queue type not implemented.\n')
    quit()
  
  if err_nnodes:
    raise ValueError('Wrong no. of nodes: ' + str(nnodes) + ' for the queue ' + q)
  
  
  #NOTE
  if 'nprocesses' in kwargs:
    del kwargs['nprocesses']
  if 'ncpus' in kwargs:
    del kwargs['ncpus']
  nthreads = PBS_File_Set_Nthreads(ncpus, nprocesses, **kwargs)
  mem_gb = str(mem_gb) + 'gb'  
  
  return time, nnodes, nprocesses, nthreads, ncpus, mem_gb


# -------------------------------------------------------------------------------


def PBS_File_Set_Time(**kwargs):
  """
  Calculate time request for a given queue type.
  
  Parameters
  ----------
  q : str 
    Queue.
  
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'PBS_File_Set_Time: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  # MUST BE HH:00:00 OTHERWISE RAISES ERROR, AT LEAST FOR THE general QUEUE
  hours = Kwarg('hours', 1, kwargs) 
  err = False
  try:
    time = kwargs['time']
    err = True
  except:
    pass
  
  if err:
    eprint(this_func + 'Error. Parameter <time> is obsolete. Use <hours> instead\n')
    quit()  
  
  time = str(hours) 
  
  if hours < 10:
    time = '0' + time
  
  time += ':00:00'
  
  
  q = Kwarg('q', 'multi', kwargs) # QUEUE
  if q == 'debug':
    time = '00:10:00'

  #print this_func + 'END'
  return time


# -------------------------------------------------------------------------------


def PBS_File_Set_Nthreads(ncpus, nprocesses, **kwargs):
  """
  COMPUTE NO. OF THREADS (VIRTUAL CORES) PER NODE 
  
  """
  this_func = this_lib + 'PBS_File_Set_Nthreads: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  try:
    nthreads = kwargs['nthreads']
    eprint(this_func + 'nthreads should rather be computed than defined.\n')
  except KeyError:
    rest = ncpus % nprocesses
    if rest == 0:
      nthreads = int(ncpus / nprocesses) # FLOOR
    else:
      raise ValueError('mpiprocs * ompthreads should equal ncpus. ' + 
                       'Otherwise qsub error. ncpus: ' + str(ncpus) + '\n')
      # ACTUALLY:
      #"PBS -lselect=N:ncpus=Y:mem=Z:mpiprocs=P:ompthreads=W
      #You should ensure that PxW <= Y or the job may be aborted by PBS"
      #(http://www.imperial.ac.uk/admin-services/ict/self-service/research-support/rcs/computing/high-throughput-computing/configuring-mpi-jobs/)
  
  return nthreads


# -------------------------------------------------------------------------------


def PBS_File_Write_Env(f, env_vars, **kwargs):
  """
  
  """
  max_ram = Kwarg('max_ram', 40, kwargs)


  # BELOW IS A FULL FOR REV690
 
  #---- List of environment variables that influence behaviour of scheduler process...
  # SCHEDULER_SHOWLEVEL=[integer]:
  #   show messages at or above given level (default=2)
  # SCHEDULER_CHECKVERBOSEFILE="yes" or "on" or [path]:
  #   Switch on verbose logging for scheduler if it finds file named
  #   "fullwave3d-verbose-scheduler" exists (off once file has gone).
  # SCHEDULER_NTHREADS=[integer]:
  #   number of threads to use for scheduler process
  #   (default determined automatically)
  # SCHEDULER_CONVERGENCE=[real]:  (frequency domain)
  #   convergence level scheduler tells solvers on slaves to end (default=1E-4)
  # SCHEDULER_WRITEUPDATE="yes" or "on":
  #   switch on dumping of update VTR(s) at each iteration
  # SCHEDULER_GLOBPRECSTAB=[real]
  #   Stabilisation factor for global preconditior (default=0.25)
  # SCHEDULER_SPREADSRCS="no" or "off":
  #   switch off non-integer source locations
  # SCHEDULER_SPREADRCVRS="no" or "off":
  #   switch off non-integer receiver locations
  # SCHEDULER_DUMPMODELPROPS="yes" or "on":
  #   dump (non-encoded) properties (after they have been extended)
  # SCHEDULER_DUMPRAWMODELPROPS="yes" or "on":
  #   dump (non-encoded) properties right after reading (not extended)
  # SCHEDULER_DUMPGRAD="yes" or "on":
  #   dump tweaked gradient(s) at each iteration
  if env_vars['dump_grad']:
    f.write('export SCHEDULER_DUMPGRAD="yes"\n')
  # SCHEDULER_DUMPPREC="yes" or "on":
  #   dump tweaked preconditioner(s) at each iteration
  if env_vars['dump_prec']:
    f.write('export SCHEDULER_DUMPPREC="yes"\n')
  # SCHEDULER_DUMPRAWGRAD="yes" or "on":
  #   dump raw gradient(s) at each iteration  
  if env_vars['dump_rawgrad']:
    f.write('export SCHEDULER_DUMPRAWGRAD="yes"\n')
  # SCHEDULER_DUMPRAWPREC="yes" or "on":
  #   dump raw preconditioner(s) at each iteration  
  if env_vars['dump_rawprec']:
    f.write('export SCHEDULER_DUMPRAWPREC="yes"\n')
  # SCHEDULER_DUMPTTIPARAMS="yes" or "on":
  #   dump anisotropy params (i.e. encoded properties), extended
  # SCHEDULER_SLAVESHAVEDATA="yes" or "on":  (time domain only)
  #   slaves read observed data file (must be visible to slaves)
  f.write('export SCHEDULER_SLAVESHAVEDATA="yes"\n')
  # SCHEDULER_SLAVESHAVEMODEL="yes" or "on":  (time domain only)
  #   slaves read model data files (must be visible to slaves)
  f.write('export SCHEDULER_SLAVESHAVEMODEL="yes"\n')    
  # SCHEDULER_SLAVESHAVEUPDATE="yes" or "on":
  #   slaves read global update files (must be visible to slaves)
  # SCHEDULER_SLAVESTIMEOUT=[integer]:
  #   Override default timeout for scheduler to mark slave as dead
  # SCHEDULER_SLAVESMAXRECOVER=[integer]:  (default 5; <1 means no reuse)
  #   Maximum number of times a dead slave can recover before not reused
  # SCHEDULER_ERRSTREAM="stdout" or "stderr" or "both"
  #   Choose where scheduler error messages will go (default stdout)
  # SCHEDULER_SHOWTIMESTAMP="yes" or "on"
  #   Include time in each line of output
  if env_vars['sched_timestamp']:
    f.write('export SCHEDULER_SHOWTIMESTAMP="yes" \n')  
  # SCHEDULER_SHOWDATESTAMP="yes" or "on"
  #   Include date in each line of output
  # SCHEDULER_SHOWHOSTNAME="yes" or "on"
  #   Include hostname in each line of output
  # SCHEDULER_FORCESLAVESGETMODEL="yes" or "on":
  #   force slaves to get model props at the start of every shot
  # SCHEDULER_FORCESLAVESGETDATA="yes" or "on":
  #   force slaves to get trace data at the start of every shot
  # 
  #---- List of environment variables that influence behaviour of slave processes...
  # SLAVES_SHOWLEVEL=[integer]:
  #   show messages at or above given level (default=2)
  # SLAVES_CHECKVERBOSEFILE="yes" or "on" or [path]:
  #   Switch on verbose logging for a slave if it finds file named
  #   "fullwave3d-verbose-slave-[SLAVEID]" exists (off once gone).
  # SLAVES_TIMEOUT=[integer]:
  #   seconds to wait for responses (default=1000)
  # SLAVES_SOURCEMULT=[real]:  (time domain)
  #   multiplication factor for source wavelet input data
  #   (may help to scale some things in certain cases)
  # SLAVES_FORCESOURCESMOOTH="no"/"off"/"yes"/"on":  (time domain)
  #   force smoothing of source wavelets (default from functional API)
  # SLAVES_WAVEFIELDSVTR=[string]:  (frequency domain)
  #   prefix name for dumped VTRs of (complex) wavefields;
  #   full name takes form:  prefix-csrefNNNNN-iterNNNNNfwdN.vtr
  if 'dump_fw' in env_vars:
    snap_step = env_vars['dump_fw']
    f.write('export SLAVES_WAVEFIELDSVTR=' + str(snap_step) + '\n')  
  # SLAVES_WAVEFIELDSVTR=[real]:  (time domain)
  #   time (seconds) at which to dump wavefield VTRs;
  #   full name takes form:  TIME-csrefNNNNN-iterNNNNNfwdN.vtr
  # SLAVES_DUMPGRAD=[string]:  (time domain)
  #   prefix name for VTR dumps of gradient contributions;
  #   full name takes form:  prefix-csrefNNNNN-iterNNNNN.vtr
  # SLAVES_DUMPPREC=[string]:  (time domain)
  #   prefix name for VTR dumps of spatial preconditioner contributions;
  #   full name takes form:  prefix-csrefNNNNN-iterNNNNN.vtr
  # SLAVES_DUMPDAT=[string]:  (time domain)
  #   prefix name for source/observed/modelled TTR dumps;
  #   full name takes form:  prefix-csrefNNNNN-iterNNNNNfwdN.ttr
  if env_vars['dump_dat']:
    f.write("export SLAVES_DUMPDAT=${project_name}-SLAVES_DUMPDAT\n")  
  # SLAVES_DUMPLOWPASS=[string]:  (time domain)
  #   prefix name for TTR dumps of traces before & after low-pass filter;
  #   full name takes form:  prefix-csrefNNNNN-iterNNNNNfwdN.ttr
  # SLAVES_DUMPAGC=[string]:  (time domain)
  #   prefix name for TTR dumps of traces before & after AGC filter;
  #   full name takes form:  prefix-csrefNNNNN-iterNNNNNfwdN.ttr
  # SLAVES_DUMPCOMPARE=[string]:  (time domain)
  #   prefix name for TTR dumps of traces just before functional compares them;
  #   full name takes form:  prefix-csrefNNNNN-iterNNNNNfwdN.ttr
  if env_vars['dump_compare']:
    f.write("export SLAVES_DUMPCOMPARE=${project_name}-SLAVES_DUMPCOMPARE\n")    
  # SLAVES_DUMPRESIDS=[string]:  (time domain)
  #   prefix name for TTR dumps of residuals;
  #   full name takes form:  prefix-csrefNNNNN-iterNNNNNfwdN.ttr
  # SLAVES_DUMPADJOINT=[string]:  (time domain)
  #   prefix name for TTR dumps of adjoint source;
  #   full name takes form:  prefix-csrefNNNNN-iterNNNNN.ttr
  # SLAVES_DUMPDENTAB=[string]: (frequency domain only)
  #   prefix name for two-column ascii dumps of density table;
  #   full name takes form:  prefix-csrefNNNNN-iterNNNNNfwdN.txt
  # SLAVES_DUMPDENSITY=[string]: (time domain only)
  #   prefix name for vtr dumps of density model built from table;
  #   full name takes form:  prefix-csrefNNNNN-iterNNNNNfwdN.vtr
  # SLAVES_DUMPMODELPROPS=[string]: (time domain only)
  #   prefix name for vtr dumps of model properties;
  #   full name takes form:  prefix-csrefNNNNN-iterNNNNNfwdN.vtr
  # SLAVES_DUMPCSREFS=[list/ranges]
  #   Restrict dumps to those CSRefs in the list, which is a
  #   comma-separated set of values/ranges, e.g. "5,15-20,35-45,60"
  # SLAVES_GARDNERPOWER=[real]
  #   Value of power for Gardner's law (default=0.25)
  # SLAVES_GARDNERFACTOR=[real]
  #   Value of multiplying factor for Gardner's law (default=310.0)
  # SLAVES_SRCPRECSTAB=[real]
  #   Stabilisation factor for local source-side preconditioner (default=10)
  # SLAVES_RCVRPRECSTAB=[real]
  #   Stabilisation factor for local receiver-side preconditioner (default=10)
  # SLAVES_LOCALSTORE=[directory-path]  (time domain & tomography only)
  #   Local directory to store forward wavefields ready for backprop.
  f.write('export SLAVES_LOCALSTORE=$work_dir \n')  
  # SLAVES_LOCALSTORETRIGGER=[real]     (time domain & tomography only)
  #   Size of array to use in RAM (GBytes) before switching to disk storage,
  #   if enabled by SLAVES_LOCALSTORE above (default=1, i.e. one Gigabyte)
  f.write('export SLAVES_LOCALSTORETRIGGER=' + str(max_ram) + ' # GB\n')
  # SLAVES_ERRSTREAM="stdout" or "stderr" or "both" (default stdout)
  #   Choose where slave error messages will go
  f.write('export SLAVES_ERRSTREAM=stderr \n')
  # SLAVES_FORCEGRIDDECI=[integer]  (range is 1 to 4)
  #   Over-ride the grid decimation factor used to store forward wavefields
  # SLAVES_FORCETIMEDECI=[integer]  (range is 1 to 40)
  #   Over-ride the time decimation factor used to store forward wavefields
  # SLAVES_SHOWTIMESTAMP="yes" or "on"
  #   Include time in each line of output
  if env_vars['slave_timestamp']:
    f.write('export SLAVES_SHOWTIMESTAMP="yes" \n')
  # SLAVES_SHOWDATESTAMP="yes" or "on"
  #   Include date in each line of output
  # SLAVES_SHOWHOSTNAME="yes" or "on"
  #   Include hostname in each line of output
 
  f.write('\n')   
  
  return f

  
# -------------------------------------------------------------------------------


def PBS_File_Calculate_Resources(nshots, **kwargs): #NOTE: NOT USED
  """
  Calculate the optimal configuration
  of resources to request for a PBS job.
  
  Parameters
  ----------
  nshots : int 
    Max. no. of shots per single iteration.
  **kwargs : keyword arguments (optional)
    Current capabilities:

  Returns
  -------
  nnodes : int 
    No. of nodes.
  nprocesses : int 
    No. of processes.
  
  Notes
  -----
  
  """
  this_func = this_lib + 'PBS_File_Calculate_Resources: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  
  try:
    single_node = kwargs['single_node']
  except KeyError:
    single_node = False    
  
  nscheduler = 1 # CONST.
  
  # NOMENCLATURE USED BY MW
  nslaves = nshots
  nprocesses = nslaves + nscheduler
  
  #if single_node:
  for nnodes in range(1, nnodes_max+1):
    print(nnodes)
    


  #print this_func + 'END'
  return nnodes, nprocesses


# -------------------------------------------------------------------------------

