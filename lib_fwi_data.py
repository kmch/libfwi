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
from lib_io_fullwave import *
##


## CONSTS
this_lib = 'lib_fwi_data.py/'

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
# PICKS
# -------------------------------------------------------------------------------


def Picks_Prepare_For_Data(sgy_files, picks_path, **kwargs):
  """
  Prepare raw picks from tlPicker to facilitate its
  merging with data files and plotting.
  
  Parameters
  ----------
  sgy_files : list 
    List of .sgy files (with paths).
  picks_path : str 
    Path where all .dat files
    with picks reside.
  shot_lines_file : str
    Name (with path) of file containing:
    line_id : shot_id_min shot_id_max
  
  Returns
  -------
  0
  
  Notes
  -----
  Pick times are in ms. To overlay them on gathers 
  we need to convert them to samples 
  (or convert gathers to ms).
  
  Ascii files will be saved in the same directory 
  as .sgy files since their names are just .sgy names 
  with different suffix.
  
  """
  this_func = this_lib + 'Picks_Prepare_For_Data: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_math_signal import Signal_Modify_Simply
  from lib_io_fullwave_CONST import gap_flag
  
  verbos = Kwarg('verbos', 0, kwargs)
  t_shift = Kwarg('t_shift', -50, kwargs) # 50 ms UP (TOWARDS 0)
  t_window = Kwarg('t_window', 1500, kwargs)
  
  name_format = Kwarg('name_format', 'Synthetic', kwargs)
  prefix = Kwarg('prefix', 'tlPick_', kwargs)
  interfix = Kwarg('interfix', '_P_line', kwargs)
  suffix = Kwarg('suffix', '.dat', kwargs)
  save_path = Kwarg('save_path', './', kwargs)
  split_picks = Kwarg('split_picks', True, kwargs)
  shot_lines_file = Kwarg('picked_lines_file', "/home/kmc3817/heavy_PhD/meta_data/shot_lines_ranges.txt", kwargs)
  
  # LOAD DICTIONARY        
  #shot_ranges = Read_Shot_Lines_Ranges(shot_lines_file)
  
  picked = 0
  notpicked = 0
  skipped = 0
  
  if split_picks:
    Picks_Split_Lines(proj, picks_path, shot_lines_file, **kwargs)
  
  
  # FOR EACH .sgy FILE
  for sgy_file in sgy_files:
    fname = Path_Leave(sgy_file)
  
    
    station_no, line_no = Data_Filename_Extract_OBS_n_Line(fname, name_format)
    #if verbos > 2:
      #print this_func, 'fname: ', fname
      #print this_func, 'station_no, line_no', station_no, line_no
    
    fname_picks = picks_path + prefix + station_no + interfix + line_no + suffix

    # READ THE PICKS IF AVAILABLE; CYCLE IF NOT
    print(this_func, 'fname_picks', fname_picks)
    picks = Picks_Read_Single_Line_File(fname_picks)
    if len(picks) == 0:
      notpicked += 1
      continue
    else:
      picked += 1
    
    picks = Picks_Assign_To_Channels(picks, sgy_file)
    ch_list = [i[0] for i in picks] # CHANNELS, ACTUALLY NOT USED ANY MORE!
    t_list = [i[1] for i in picks]  # PICKED TIMES
    
    # SKIP TOO EARLY PICKS      
    if len(t_list) == 0 or (max(t_list) * 1000) <= abs(t_shift):
      skipped += 1
      #eprint(this_func + 'Skipping the file because all picks <= 0 after shift\n')
      continue
    
    # PREPARE TWO LISTS FOR OF TOP (ACTUAL PICKS) AND BOTTOM PICKS (TOP + t_window)
    # 1. MOVE t_shift EARLIER TO CAPTURE THE FIRST SWING (ALLOW FOR PICKING ERROR) 
    # 2. CONVERT SEC->MSEC
    # 3. ASSIGN SPECIAL VALUES FOR GAPS FLAGGED WITH A GAP-FLAG
    t1_list = Signal_Modify_Simply(t_list, scale=1000, shift=t_shift, flagged=[gap_flag, 0])  
    t2_list = Signal_Modify_Simply(t_list, scale=1000, shift=t_shift+t_window, flagged=[gap_flag, 0]) 
    
    
    # SAVE BOTH TOP (ACTUAL PICK) AND BOTTOM PICK (TOP + t_window) TO RESPECTIVE FILES
    for i, line in enumerate([t1_list, t2_list]):
      fname = sgy_file[ :-len('.sgy')] + '_t' + str(i+1) + '.ascii' #'_shift' + str(abs(t_shift) / 1000.) + '_twin' + str(t_window / 1000.) + '.ascii'
      f = open(fname, 'w')
      for value in line:
        f.write(str(value) + '\n')
      f.close()
    
  print(this_func, 'No. of picked lines', picked)
  print(this_func, 'No. of skipped lines', skipped)
  print(this_func, 'No. of not picked lines', notpicked)
    
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Picks_Read_Single_Line_File(fname, **kwargs):
  """
  Read picks from a file of a fixed format.
  
  Parameters
  ----------
  fname : str 
    File name.
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  picks : list 
    List of the form (channel_id, picked_time).
    Type of channel_id is preserved (str) but 
    picked_time is a number.
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Picks_Read_Single_Line_File: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  picks_format = Kwarg('picks_format', 'tlPick', kwargs)
  
  # SPECIFY STRUCTURE OF A GIVEN FILE-FORMAT
  if picks_format == 'tlPick':
    X_read = lambda c : [i[0] for i in c]
    T_read = lambda c : [float(i[1]) for i in c]
  
  else:
    eprint(this_func + 'Error. Unknown format of the pick-file\n')
    quit()
  
  # READ IF THE FILE EXISTS
  try:
    c = Read_File(fname, verbos=False)
    X = X_read(c)
    T = T_read(c)
    picks = list(zip(X, T))
  except:
    picks = []
  
  #print this_func + 'END'
  return picks


# -------------------------------------------------------------------------------


def Picks_Assign_To_Channels(picks, file_sgy, **kwargs):
  """
  For each channel of the data-file assign:
  - a pick  - if present in the picks
  - an error-flag - if not.
  
  Parameters
  ----------
  picks : list 
    List of the form (channel_id, picked_time).
    Type of channel_id is preserved (str) but 
    picked_time is a number.
  sgy_file : str 
    Data file for which the picks 
    are prepared.
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  Returns
  -------
  channels_picks : list 
    Updated list of picks.
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Picks_Assign_To_Channels: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_io_fullwave_CONST import key_shot_id, gap_flag
  
  key_chnl = Kwarg('key_chnl', key_shot_id, kwargs)
  
  suffix = '_' + key_chnl + '.txt' # FORMAT USED BY sugethw.sh
  file_channels = file_sgy[ :-len('.sgy')] + suffix
  
  try:
    c = Read_File(file_channels, verbos=0)
    channels = [i[0] for i in c] 
  except:
    o = Bash('su_gethw.sh fldr '+file_sgy)
    c = Read_File(file_channels)
    channels = [i[0] for i in c] 
    
  channels_picks = [[x, gap_flag] for x in channels]
  
  # ASSIGN PICKS FOR THESE CHANNELS WHERE THEY EXIST, THE REST WILL BE FLAGGED
  for i, ch in enumerate(channels):
    for pick in picks:
      x, t = pick
      if x == ch: # COMPARING 2 STRINGS
        channels_picks[i] = [x, t]
  
  #print this_func + 'END'
  return channels_picks


# -------------------------------------------------------------------------------
# VARIA
# -------------------------------------------------------------------------------


def Data_Filename_Extract_OBS_n_Line(fname, name_format, **kwargs):
  """
  Extract OBS and shot-line IDs from a filename 
  of a data-file.
  
  Parameters
  ----------
  fname : str 
    Name of the file (without path).
  name_format : str 
    One of the precisely formats
    of naming data-files used by UO, Jo
    and myself.
  
  Returns
  -------
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Data_Filename_Extract_OBS_n_Line: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  fname = Path_Leave(fname) # JUST IN CASE
  
  if name_format == 'UO':
    char_pos = len('MGL1521_?')
    station_no = fname[char_pos : char_pos + 3]
    char_pos = len('MGL1521_????_?_')
    line_no = fname[char_pos : -len('.sgy')]
  
  elif name_format == 'Synthetic':
    try:
      proj_name = kwargs['proj_name']
    except KeyError:
      eprint(this_func + 'Error. Need proj_name for this name_format\n')
      quit()

  elif name_format == 'jm': # Jo's PHASE COMPARISONS
    char_pos = 1
    station_no = fname[char_pos : char_pos + 3]
    char_pos = char_pos + 7
    line_no = fname[char_pos : -len('.sgy')]
  
  else:
    eprint(this_func + 'Error. Unknown name_format\n')
    quit()
  
  #print this_func, 'station_no, line_no', station_no, line_no
  
  #print this_func + 'END'
  return station_no, line_no


# -------------------------------------------------------------------------------


def Prepare_SU_Input_List(py_list, **kwargs):
  """
  Convert standard Python list into 
  a string read by a command line => Seismic Unix.
  
  Parameters
  ----------
  py_list : list 
    List of values.
  
  **kwargs : keyword arguments, optional
  
  Returns
  -------
  su_list : str 
    String with commas.
  
  Notes
  -----
  Format of the string for list [val_1, val_2, ..., val_n]
  'val_1,val_2,val_3,val_n'
  
  """
  this_func = this_lib + 'Prepare_SU_Input_List: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  su_list = ''
  for l in py_list:
    su_list += str(l) + ','
  
  # DELETE THE LAST COMMA
  su_list = su_list[ :-1]
  
  #print this_func + 'END'
  return su_list 


# -------------------------------------------------------------------------------


def Read_Stations_Picked_Lines(file_name, **kwargs):
  """
  Read IDs of picked lines for each station. 
  
  Parameters
  ----------
  file_name : string
    File name. It should include extension. 
    It can include path if needed.
  **kwargs : keyword arguments, optional
  
  Returns
  -------
  picked_lines : dict
    Dictionary where key is a station no.
    and the value is the list of picked lines IDs.
  
  Notes
  -----
  File structure:
  
  station1_no : pl1 pl2 ...
  station2_no : pl1 pl2 ...
  ...
  
  """
  this_func = this_lib + 'Read_Stations_Picked_Lines: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  #fname = Get_Files(meta_path, picked_lines_file)[0]
  c = Read_File(file_name)
  
  picked_lines = {}
  for l in c:
    sid = l[0]
    plines = l[2: ]
    picked_lines[sid] = plines
  
  #print this_func + 'END'
  return picked_lines


# -------------------------------------------------------------------------------





# RUBBISH 

def Picks_Line_Full_Info(line_picks, line_no, shot_ranges, **kwargs):
  """
  
  
  Parameters
  ----------
  
  
  
  Returns
  -------
  
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Picks_Line_Full_Info: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  shot_min, shot_max = shot_ranges[line_no]
  #all_line_shots = range(shot_min, shot_max + 1)
  full_line_X = list(range(shot_min, shot_max + 1))
  nshots = len(full_line_X)
  full_line_T = np.ones(nshots) * err_value
  
  for i, x in enumerate(full_line_X):
    for pick in line_picks:
      if int(pick[0]) == int(x):
        full_line_T[i] = float(pick[1])

  #print this_func + 'END'
  return full_line_X, full_line_T