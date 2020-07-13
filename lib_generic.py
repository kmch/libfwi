"""
Author: Kajetan Chrapkiewicz, 2018. 
All rights reserved. Ask for permision writing to K.Chrapkiewicz17@imperial.ac.uk.

This library provides framework procedures for...
At the moment fully implemented is:

1. ...

"""

## MODULES 
import numpy as np
from lib_generic_CONST import verbos_func


## CONSTS
this_lib = 'lib_generic.py/'


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
# SYSTEM
# -------------------------------------------------------------------------------


def Duplicate(source, destination, **kwargs):
  """
  Copy / move / link.
  
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
  
  this_func = this_lib + 'Duplicate: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  dupl_cmd = Kwarg('dupl_cmd', 'cp', kwargs)
  
  if verbos > 0:
    print((this_func + 'Duplicating ' + source + 
          ' to ' + destination + ' using ' + dupl_cmd + 
          ' command.'))
  
  command = dupl_cmd + ' ' + source + ' ' + destination
  o, e  = Bash2(command, **kwargs)

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Cp(source, destination='./', **kwargs):
  """
  Copy.
  
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
  
  this_func = this_lib + 'Cp: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  o = Bash('cp ' + source + ' ' + destination)
  print(this_func, o)

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Bash(command, **kwargs):
  """
  Call Bash command/script.
  This is actually a wrapper for subprocess/check_output.
  
  
  Parameters
  ----------
  command : str 
    Command as a 1 string.
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  oupt : str 
    Output.
  
  Notes
  -----
  
  
  """
  
  this_func = this_lib + 'Bash: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from subprocess import check_output
  
  path = Kwarg('path', './', kwargs)
  
  
  i = [_f for _f in command.split(' ') if _f]
  #print this_func, 'input of the callee: ', inpt
  
  try: 
    o = check_output(i)#, cwd=path, shell=True) # cwd - PATH WHERE IT WILL BE EXECUTED
  except:
    eprint(this_func + 'Execution of the command: ' + command + ' failed - check it manually.\n')
    o = 'No output'
    pass
  #print o
  
  # FIXME: NOT SURE WHAT shell=True IS DOING
  
  #except OSError:
  #  eprint(this_func + 'Problems with executing the command. Outputting empty array...\n')
  #  oupt = []
    
  #print this_func + 'END'
  return o


# -------------------------------------------------------------------------------


def Bash2(command, **kwargs):
  """
  Call Bash command/script in an improved way.
  
  Parameters
  ----------
  command : str 
    Command as a 1 string.
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  out : str 
    Standard output.
  err : str
    Standard error.
  
  Notes
  -----
  It respects wildcards.
  
  """
  
  this_func = this_lib + 'Bash2: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from subprocess import Popen, PIPE
  
  path = Kwarg('path', './', kwargs)
  
  command = 'cd ' + path + ' && ' + command # NEW
  
  if verbos > 5:
    print(this_func, 'Running ', command)
  
  proc = Popen(command, shell=True,
               stdout=PIPE, stderr=PIPE)
  
  out, err = proc.communicate()
  
  if err != '' and verbos > 0:
    eprint(this_func + 'Non-empty stderr of the command ' + command + ': \n')
    eprint(err)
  
  #print this_func + 'END'
  return out, err


# -------------------------------------------------------------------------------


def Exists(fname, **kwargs):
  """
  Check if the file exists.
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  Returns
  -------
  exists : bool 
    True if the file exists.
  
  Notes
  -----
  It is a tiny wrapper.
  
  """
  
  this_func = this_lib + 'Exists: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  #from os.path import isfile
  import os
  
  exists = False 
  
  if os.path.exists(fname):
    exists = True  
  
  else:
    exists = False
    if verbos > 0:
      eprint(this_func + 'File: ' + fname + ' does not exist.\n')

  #print this_func + 'END'
  return exists


# -------------------------------------------------------------------------------


def Paths_Prepare(machine, **kwargs):
  """
  Prepare dict. of paths 
  for a given machine.
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  paths : dict
    Dictionary containing all the useful 
    paths e.g. for jupy notebooks.
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Paths_Prepare: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  if machine == "hp":
    user = "kajetan"
  elif machine == "kmc":
    user = "kmc3817"
  elif machine == "toshiba":
    user = "chkajetan"
  else:
    print("Error. Unknown machine")
  
  paths = {}
  
  # GENERAL
  paths['home'] = "/home/" + user + "/"
  paths['ext_drive'] = "/media/" + user + "/TOSHIBA EXT/"
  paths['phd_local'] = paths['home'] + 'heavy_PhD/' #ext_drive + "BACKUP_PhD_local_07-05-2018/"  
  paths['phd_light'] = paths['home'] + 'Dropbox/light_PhD/' #NOTE THIS IS ACTUALLY FOR laptops ONLY
  paths['projects'] = paths['phd_local'] + "PROJECTS/"
  paths['proteus'] = paths['home'] + "Dropbox/PhD/team_PROTEUS/" # CAN"T BE ~/Dropbox
  paths['figs1'] = paths['phd_local'] + "figures1_draft/"
  paths['figs2'] = paths['phd_local'] + "figures2_better/"
  paths['figs3'] = paths['phd_local'] + "figures3_best/"
  
  # META DATA
  paths['meta'] = paths['phd_local'] + "meta_data/"  
  paths['stations_names_file'] = "stations_file_names.txt"
  paths['stations_IDs_file'] = "stations_IDs.txt"
  paths['picks_info_file'] = "Santorini_Station_Picks.csv"
  
  # DATA ETC.
  paths['obs_data'] = "/media/kmc3817/DATADRIVE1/heavy_PhD/DATA/Santorini_2015/seismic/OBS/segy_local_coords/"
  paths['land_data'] = "/media/kmc3817/DATADRIVE1/heavy_PhD/DATA/Santorini_2015/seismic/land/Santorini/local_coords_headers_by_MP/"
  paths['picks'] =  paths['phd_local'] + "picks/first_arriv/UO/time_corrected/"
  paths['mod'] = paths['phd_local'] + "start_mods/"
  paths['source_wavelet'] = paths['phd_local'] + "source_wavelets/"
  paths['fs'] = paths['projects'] + "topo/"  
  paths['shot_info_file'] = "shot_info_16-08-19.csv"
  paths['source_wavelet_file'] = "mp_source_3d_7.5Hz_17-10-17.sgy"
  paths['old_mod_file'] = "Model_170710_171413_it5_corr.sgy"
  paths['new_mod_file'] = "Ben_whole_model_24-04-18.sgy"
  paths['mod_file'] = paths['new_mod_file'] # FIXME: TMP!
  paths['fs_file_grd'] = "santorini_merged_local_50.grd" 
  
  # TEMPLATES
  paths['templates'] = paths['home'] + 'Dropbox/light_PhD/fwi_templates/'
  paths['template_runfile'] = 'template-Runfile.key'

  #print this_func + 'END'
  return paths


# -------------------------------------------------------------------------------
# I/O
# -------------------------------------------------------------------------------


def Read_ext(fname, **kwargs): # READ ANY EXTENSION
  """
  Read any extension.
  
  Parameters
  ----------
  
  Returns
  -------
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Read_Ext: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_io_fullwave import Read_vtr
   
  ext = Ext(fname)
  if ext == 'vtr':
    from lib_io_fullwave import Read_vtr
    n1, n2, n3, Z = Read_vtr(fname, **kwargs)
  
  elif ext == 'sgy':
    from lib_io_fullwave import Read_sgy
    Z = Read_sgy(fname, **kwargs)
    
  else:
    raise IOError('Unknown extension: ' + ext)
    
  #print this_func + 'END'   
  return Z


# -------------------------------------------------------------------------------


def Read_Array(fname_or_Z, **kwargs): # , function_read,
  """
  Generic wrapper for reading function returning 
  array if the argument is already one.
  
  Parameters
  ----------
  
  Returns
  -------
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Read_Array: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_io_fullwave import Read_vtr
   
  if verbos > 1:
    print(this_func, 'fname_or_Z', fname_or_Z)
    
  if isinstance(fname_or_Z, str):
    Z = Read_ext(fname_or_Z, **kwargs)  
  
  elif type(fname_or_Z) == type(np.array([])):
    Z = fname_or_Z
  
  else:
    raise TypeError('Unknown type of fname_or_Z: ' + str(type(fname_or_Z)))

  #print this_func + 'END'   
  return Z


# -------------------------------------------------------------------------------


def eprint(*args, **kwargs):
  """
  Print to stderr.
  
  """
  import sys
  sys.stderr.write(*args, **kwargs)


# -------------------------------------------------------------------------------


def Get_Files(path, pattern, **kwargs):
  """
  Get list of files located in path 
  defined by the REGEX pattern.
  
  Parameters
  ----------
  path : str 
    Full path to the directory supposed to
    contain the files.
  pattern : str 
    REGEX pattern defining names of the files;
    it can contain wildcards, e.g.
    pattern = 'project-*vtr'
  
  Returns
  -------
  files_list : list 
    Sorted alphabetically list of files names each of which includes
    the full path (!).  
  
  Notes
  -----
  Used in Jupyter notebooks.
  
  """
  this_func = this_lib + 'Get_Files: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from os import listdir as ls
  from fnmatch import filter as expand
  
  files_list = sorted([path + '/' + f for f in expand(ls(path), pattern)])
  
  if len(files_list) == 0:
    if verbos > 0:
      eprint(this_func + 'No files matching ' + 
             pattern + ' in ' + path + '\n')  
  
  if verbos > 6:
    print('pattern', pattern)
    print(this_func, 'path: ', path, expand(ls(path), pattern))
    print(this_func, 'files_list: ', files_list)
  
  #print this_func + 'END'   
  return files_list


# -------------------------------------------------------------------------------


def Read_File(file_name, **kwargs):
  # READ A FILE NON-ZERO LINES, SPLIT THEM WITH None AND RETURN CONTENT 
  #print 'lib_generic/Read_File: START'
  
  verbos = Kwarg('verbos', 1, kwargs)
  
  try:
    f = open(file_name,'r')    
  except IOError:
    raise IOError('lib_generic/Read_File: Error. File: ' + file_name + ' not found.')
    
  content = []
  for line in f:
    line = line.split(None)
    if len(line) != 0:
      content.append(line)
  f.close()
  
  #print 'lib_generic/Read_File: END'  
  return content


# -------------------------------------------------------------------------------


def Save_Dict(file_name, dictionary, **kwargs):
  """
  Save a dictionary data-structure to a file.
  
  Parameters
  ----------
  file_name : str 
    Full name of the file including 
    a path if needed.
  dictionary : dict 
    Dictionary to save.
  
  Returns
  -------
  0
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Save_Dict: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START') 
  
  
  
  if verbos > 2:
    print(this_func, 'Saving: ', file_name)
  
  f = open(file_name, 'w')
  for i in dictionary:
    f.write(i + ' : ' + str(dictionary[i]) + '\n')
  f.close()
  
  #print this_func + 'END' 
  return 0


# -------------------------------------------------------------------------------


def Save_bin(file_name_core, nx, ny, Z, **kwargs):
  """
  Save as a .bin file.
    
  Parameters
  ----------
  file_name_core : str 
    Name of the file without extension 
    including a path if needed.
  nx : int 
  
  ny : int 
  
  Z : array
  
  
  Returns
  -------
  0
  
  Notes
  -----
  
  """  
  this_func = this_lib + 'Save_bin: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  from scipy.io import FortranFile
  
  f = FortranFile(file_name_core + '.bin', 'w')
  for i in range(nx):
    for j in range(ny):
      trace = Z[i][j]
      f.write_record(np.array(trace, dtype = np.float32))  
  f.close()
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Rewrite_File(file_name, new_file_name, **kwargs): 
  """
  Rewrite content of a file into 
  a new one without changing anything.
  It is equivalent to copying but it 
  does not invoke any console command.
  
  Parameters
  ----------
  file_name : str 
    Full name of the file to rewrite from 
    including a path if needed.
  new_file_name : str 
    Full name of the file to write into
    including a path if needed.  
  
  Returns
  -------
  0
  
  Notes
  -----
  It avoids calling system-dependent
  function from 'sys' library etc.
  
  It probably doesn't work for binaries.
  
  """  
  this_func = this_lib + 'Rewrite_File: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  content = Read_File_Not_Split(file_name)
  fn = open(new_file_name, 'w')
  
  with open(file_name, 'r') as f:
    for line in f:
      fn.write(line)
  
  fn.close()
  
  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Read_File_Not_Split(file_name, **kwargs):
  """
  Read a file exactly as it is.
  
  Parameters
  ----------
  file_name : str 
    Full name of the file including 
    a path if needed.
  
  Returns
  -------
  content : list
    Content of a file as 
    a list of lines.
  
  Notes
  -----
  Most basic file-reading routine.
  
  FIXME: Think maybe of a complex Read_File(file_name, **kwargs) 
  function...
  
  """  
  this_func = this_lib + 'Read_File_Not_Split: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  f = open(file_name, 'r')
  
  content = []
  for line in f:
    content.append(line)
  
  f.close()
  
  #print this_func + 'END'
  return content


# -------------------------------------------------------------------------------
# PARSERS
# -------------------------------------------------------------------------------


def Ext_Change(fname, ext_new, **kwargs):
  """
  Change the name of file 
  so that it has different extension.
  
  Parameters
  ----------
  fname : str 
    Assumed format: name.extension
  ext_new : str 
    Without full stop.
    
  Returns
  -------
  fname : str
    New fname
  
  
  Notes
  -----
  Work out this returning scheme. # FIXME
  
  """
  this_func = this_lib + 'Path_Leave: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START') 
  
  ext_old = Ext(fname)
  fname = fname[ :-len(ext_old)] + ext_new
  
  #print this_func + 'END'  
  return fname


# -------------------------------------------------------------------------------


def Path_Leave(path, **kwargs):
  """
  Extract a file's name from a whole 
  path. 
  
  Parameters
  ----------
  path : str 
    Path including the file name.
    
  Returns
  -------
  file_name : str
    Extracted file name.
  
  
  Notes
  -----
  Work out this returning scheme. # FIXME
  
  """
  this_func = this_lib + 'Path_Leave: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START') 
  
  import ntpath
  
  head, tail = ntpath.split(path)
  
  #print this_func + 'END'  
  return tail or ntpath.basename(head)


# -------------------------------------------------------------------------------


def Path(path, **kwargs):
  """
  path from full file name 
  
  Parameters
  ----------
  path : str 
    Path including the file name.
  
  Notes
  -----
  Work out this returning scheme. # FIXME
  
  """
  this_func = this_lib + 'Path_Leave: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START') 
  
  import ntpath
  
  head, tail = ntpath.split(path)
  
  head += '/'
  
  #print this_func + 'END'  
  return head


# -------------------------------------------------------------------------------


def Ext(file_name, **kwargs):
  """
  Extract the extension from the 
  file name.
  
  Parameters
  ----------
  file_name : str 
    It can include path. It shouldn't include 
    more than one full stop that exactly separates
    the extension from the proceeding part of the 
    name.
  
  Returns
  -------
  ext : str 
    Extension (without the full stop).
  
  Notes
  -----
  
  """
  
  this_func = this_lib + 'Ext: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  split = [_f for _f in file_name.split('.') if _f]
  
  if len(split) >= 2:
    ext = split[-1]
  
    if len(split) > 2:
      eprint(this_func + 'File name contains more full stops. Only the last bit will be treated as extension.')
  
  else:
    eprint(this_func + 'Error. Bad splitting. Check if a full stop separates extension.\n')
    quit()

  #print this_func + 'END'
  return ext



# -------------------------------------------------------------------------------


def Strip(fname, **kwargs):
  """
  Strip the extension off the file name.
  
  Parameters
  ----------
  fname : str 
    String of the form:
    <string_without_dots>.<extension>
      
  Returns
  -------
  nfname: str 
    Stripped string.
  
  Notes
  -----
  It assumes the string:
    <string_without_dots>.<extension>
    
  """
  this_func = this_lib + 'Strip: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  # THIS ALONE DOESN'T WORK FOR PATHS ./fname.ext
  #nfname = filter(None, fname.split('.'))[0] 

  ext = Ext(fname)
  suffix = '.' + ext
  
  if 'suffix' in kwargs:
    del kwargs['suffix']
  nfname = Core(fname, ' ', suffix, **kwargs)
  
  #print this_func + 'END'
  return nfname


# -------------------------------------------------------------------------------


def Del(kwarg_name, kwargs, warning=False):
  """
  Delete kwarg from the dict.
  
  Parameters
  ----------
  kwarg_name : str
    Name of the argument.
  kwargs : dict 
    Dictionary of all keyword 
    arguments.
  
  Returns
  -------
  kwargs : dict
  
  Notes
  -----
  Name of the function is as short as possible
  because it's frequently typed.
  
  We can't not include kwargs as an argument.
  
  """  
  this_func = this_lib + 'Del: '
  #print this_func + 'START'

  try:
    del kwargs[kwarg_name]
  except KeyError:
    if warning:
      eprint('Value of ' + kwarg_name + ' was not specified.\n')
  
  
  #print this_func + 'END'
  return kwargs


# -------------------------------------------------------------------------------


def Kwarg(kwarg_name, default_value, kwargs, warning=False):
  """
  Read a value of a keyword argument
  assigning a default if it's missing 
  from the dict **kwargs.
  
  Parameters
  ----------
  kwarg_name : str
    Name of the argument.
  default_value : any 
    The value which will be assigned 
    if kwarg_name is not found 
    in kwargs.
  kwargs : dict 
    Dictionary of all keyword 
    arguments.
  
  Returns
  -------
  value : any 
    Value assigned to kwarg_name.
  
  Notes
  -----
  Name of the function is as short as possible
  because it's frequently typed.
  
  We can't not include kwargs as an argument.
  
  """  
  this_func = this_lib + 'Kwarg: '
  #print this_func + 'START'

  try:
    value = kwargs[kwarg_name]
  except KeyError:
    value = default_value 
    if warning:
      eprint('Value of ' + kwarg_name + ' not specified. Assuming ' 
             + default_value + '...\n')
  
  #print this_func + 'END'
  return value


# -------------------------------------------------------------------------------

 
def Core(string, prefix, suffix, **kwargs):
  """
  Get a core of the string i.e. strip 
  both its prefix and suffix.
  
  Parameters
  ----------
  string : str 
    String to strip.
  prefix : str 
    Prefix contained in the string.
  suffix : str 
    Suffix contained in the string.
  
  Returns
  -------
  core : str 
    Stripped string.
  
  Notes
  -----
  Prefix/suffix must by space (' ', not '') if we dont want any
  
  FIXME DOUBLE-CHECK ESPECIALLY FOR MULTIPLE SEPARATORS
  
  CAUTION DOES NOT WORK FOR STRINGS (E.G. PATHS) CONTAINING WHITE SPACES 
  
  """
  this_func = this_lib + 'core: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  beginning = list(filter(None, string.split(suffix)))[0]
  core = list(filter(None, beginning.split(prefix)))[0]
  
  #print this_func + 'END'
  return core 


# -------------------------------------------------------------------------------


def Split(string, separator, **kwargs):
  """
  Simple wrapper around filter function
  to split a string into chunks
  separated by a separator-substring.
  Chunks do not contain the separator.
  
  Parameters
  ----------
  string : str 
    String to strip.
  separator : str 
    Substring separating consecutive
    chunks.
  
  Returns
  -------
  chunks : list 
    List of chunks.
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Split: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  chunks = [_f for _f in string.split(separator) if _f]
  
  #print this_func + 'END'
  return chunks


# -------------------------------------------------------------------------------


def Is_Int(r, precision, **kwargs):
  """
  Check if the float-format number
  is a whole number.
  
  Parameters
  ----------
  r : float 
    Real number to check.
  presision : float 
    Precision of the check i.e. the 
    smallest difference between a
    non-integer and its nearest integer.
  
  Returns
  -------
  answer : bool 
    True if r is a whole number 
    within precision.
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Is_Int: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START') 
  
  if (abs(r - int(r)) < precision) or (abs(int(r) + 1 - r) < precision):
    answer = True
  else:
    answer = False
  
  #print this_func + 'END'
  return answer


# -------------------------------------------------------------------------------


def Nint(r, **kwargs):
  """
  Return nearest int of the float.
  
  Parameters
  ----------
  r : float 
    Real number to check.
  
  Returns
  -------
  nint : int 
    Nearest integer.
  
  Notes
  -----
  x.5 will be round down to x.
  
  """
  this_func = this_lib + 'Nint: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START') 
  
  lint = int(r) # LOWER
  uint = int(r) + 1 # UPPER 
  
  ldiff = abs(r - lint)
  udiff = abs(r - uint)
  
  if ldiff < udiff:
    nint = lint
  
  elif udiff < ldiff:
    nint = uint
  
  elif udiff == ldiff:
    nint = lint # x.5 will be round down to x.
    
  else:
    print(this_func, 'Error. Weird case.')
    quit()
  
  #print this_func + 'END'
  return nint


# -------------------------------------------------------------------------------
# DATA-STRUCTURES HANDLERS
# -------------------------------------------------------------------------------


def Array_Interleave(Z1, Z2, **kwargs):
  """ 
  Create an array composed of 
  interleaved arrays Z1 & Z2.
  
  Parameters
  ----------
  
  Z1, Z2 : 2D arrays
    Arrays to interleave.
  
  **chunk_size : int 
    No. of columns of 1 array
    before being proceeded 
    by 2nd array etc.
    
  Returns
  -------
  Z : 2D array
    Result array.
  
  Notes
  -----
  ASSUMPTIONS:
  - Z1, Z2 are already sliced => 2D
  - they also have the same shapes NOTE
  
  """
  this_func = this_lib + 'Array_Interleave: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START') 
    
  chunk_size = Kwarg('chunk_size', 10, kwargs)    

  Z = np.array(Z1)
  ncols = Z.shape[0]
  
  if ncols < 2 * chunk_size:
    eprint(this_func + 'No. of columns=' + str(ncols) + 
           ' < 2 * chunk_size! Outputting empty array\n')
    return []
  
  nchunks = ncols / chunk_size / 2
  
  for i, Zi in enumerate([Z1, Z2]):
    i_start = i * chunk_size
    for j in range(nchunks):
      i1 = i_start + j * 2 * chunk_size
      i2 = i_start + j * 2 * chunk_size + (chunk_size - 1) 
      Z[i1 : i2] = Zi[i1 : i2]

  Z = np.array(Z)

  #print this_func + 'END'
  return Z


# -------------------------------------------------------------------------------


def Array_Interleave_Split(Z1, Z2, **kwargs):
  """ 
  Similar to Array_Interleave but 
  it creates 2 separate arrays 
  completing each other. Each 
  has gaps (zero values) meant to
  be filled by the other array.
  This is useful if you want to 
  different plotting parameters 
  (e.g. colormaps) for each.
  
  Parameters
  ----------
  
  Z1, Z2 : 2D arrays
    Arrays to interleave.
  
  **chunk_size : int 
    No. of columns of 1 array
    before being proceeded 
    by 2nd array etc.
    
  Returns
  -------
  Z : 2D array
    Result array.
  
  Notes
  -----
  ASSUMPTIONS:
  - Z1, Z2 are already sliced => 2D
  - they also have the same shapes NOTE
  
  """
  this_func = this_lib + 'Array_Interleave_Split: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START') 
    
  chunk_size = Kwarg('chunk_size', 10, kwargs)    
  
  # SAME SHAPES ASSUMED
  #nZ1 = np.array(Z1)
  #nZ2 = np.array(Z1)
  ncols = Z1.shape[0]
  
  if ncols < 2 * chunk_size:
    raise ValueError('No. of columns must be >= 2 * chunk_size!')
  
  nchunks = ncols / chunk_size / 2
  
  # CREATE GAPS (ZEROES-PATCHES) IN EACH ARRAY
  for i, Zi in enumerate([Z2, Z1]): #NOTE: SWAPPED Z1, Z2 RELATIVE TO Array_Interleave 
    i_start = i * chunk_size
    for j in range(nchunks):
      i1 = i_start + j * 2 * chunk_size
      i2 = i_start + j * 2 * chunk_size + (chunk_size - 1) 
      Zi[i1 : i2] *= 0

  #print this_func + 'END'
  return Z1, Z2


# -------------------------------------------------------------------------------


def Array2D_Convert_To_3D(ZZ, **kwargs):
  """
  Convert 2D-surface-embedded-in-3D
  data (nx, ny) into 3D array (nx, ny, 1)
  
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
  
  this_func = this_lib + 'Array2D_Convert_To_3D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  nx, ny = ZZ.shape
  ZZZ = np.zeros([nx, ny, 1])
  
  ZZZ[:, :, 0] = ZZ
  
  #print ZZZ.shape

  #print this_func + 'END'
  return ZZZ


# -------------------------------------------------------------------------------


def Array1D_Convert_To_3D(Z, **kwargs):
  """
  Convert 1D array (e.g. time series)
  data (nx,) into 3D array (nx, 1, 1)
  
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
  
  this_func = this_lib + 'Array1D_Convert_To_3D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  nx = len(Z)
  ZZZ = np.zeros([1, 1, nx]) # DOUBLE-CHECKED - OK, COMPATIBLE!
  
  ZZZ[0, 0, :] = Z
  
  #print ZZZ.shape

  #print this_func + 'END'
  return ZZZ


# -------------------------------------------------------------------------------


def Equalize_Series_Lengths(series1, series2, mode, **kwargs):
  """
  
  Parameters
  ----------

  
  Returns
  -------

  
  Notes
  -----
  It assumes that first samples of both series
  correspond to the same x.
  
  """
  this_func = this_lib + 'Equalize_Series_Lengths: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  l1, l2 = len(series1), len(series2)
  
  if mode == 'trim':
    if l1 > l2:
      series1 = series1[ :l2]
    else:
      series2 = series2[ :l1]
  
  elif mode == 'pad':
    if l1 > l2:
      nseries2 = np.zeros(len(series1))
      nseries2[ :len(series2)] = series2
      series2 = nseries2
    else:
      nseries1 = np.zeros(len(series2))
      nseries1[ :len(series1)] = series1
      series1 = nseries1
  
  else:
    eprint(this_func, 'Mode ' + mode + 'not implemented yet.')
    quit()  
    
  #print this_func + 'END'
  return [series1, series2]


# -------------------------------------------------------------------------------


def Get_Chunks(l, n, **kwargs):
  # DIVIDE LIST l INTO n-ELEMENT PIECES RETURNE AS AN ITERABLE
  for i in range(0, len(l), n):
    yield l[i : i + n]    


# -------------------------------------------------------------------------------


def Flatten_List(list2D, **kwargs):
  # TURN A LIST OF SUBLISTS INTO A FLAT (NOT NESTED) LIST
  import itertools
  return list(itertools.chain.from_iterable(list2D))


# -------------------------------------------------------------------------------


def Delete_List_Midpoint(li, **kwargs):
  """
  Delete the middle point of the 
  list of odd no. of elements.
  
  Parameters
  ----------
  li : list
    List to treat. It can have 
    any length.
  
  Returns
  -------
  nli : list
    Processed list.
  
  Notes
  -----
  It deals with both possible list-lengths.
  
  """
  this_func = this_lib + 'Delete_List_Midpoint: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START') 
  
  lngth = len(li)
  
  # EVEN
  if lngth%2 == 0:
    eprint(this_func + 'Even no. of elements. No modif. applied.') 
    nli = li 
  
  # ODD  
  elif lngth%2 == 1:
    floor_half_len = int(lngth / 2)
    # 1ST HALF
    nli = list(li[ :floor_half_len])
    # 2ND HALF
    nli += list(li[floor_half_len + 1: ])
  
  #print this_func + 'END'
  return nli


# -------------------------------------------------------------------------------


def Split_Tuples(tuples, **kwargs):
  """
  Convert a list of 3D tuples into 3 lists
  of x, y, z coordinates.
  
  Parameters
  ----------
  tuples : list
    List of tuples of the form [x, y, z].
  
  Returns
  -------
  X, Y, Z : list 
    Lists of coordinates.
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Split_Tuples: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  X, Y, Z = [], [], []
  for tupl in tuples:
    X.append(tupl[0])
    Y.append(tupl[1])
    Z.append(tupl[2])
  
  #print this_func + 'END'
  return X, Y, Z


# -------------------------------------------------------------------------------


def Key_Remove(dictionary, key, **kwargs):
  """
  
  """
  this_func = this_lib + 'Key_Remove: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  r = dict(dictionary)
  
  if key in r:
    del r[key]
  
  #print this_func + 'END'
  return r


# -------------------------------------------------------------------------------
# FUNCTION EVALUATORS
# -------------------------------------------------------------------------------


def Vec(func, X, *args, **kwargs):
  """
  Make a scalar-function act on 
  an N-D array.
  
  Parameters
  ----------
  func : function
    It HAS TO to take arguments in
    the following order:
    func(x, *args, **kwargs) where 
    x is a scalar counterpart of X.
    
    NOTE: for func(x, y, z, **kwargs) 
    y and z will be read automatically
    from *args. NO NEED for Vec_2D etc.
    
  X : array
    Array counterpart of x.
  
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  X_out : array 
    The array A transformed element-wise 
    by the function func.
  
  Notes
  -----
  It is required e.g. when the function
  contains if statements.
  
  It should work in ND (Y, Z, ... are just 
  args).
  
  
  """
  
  this_func = this_lib + 'Vec: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  #print this_func, 'args', args
  #print this_func, 'kwargs', kwargs
  
  v_func = np.vectorize(func, excluded=list(kwargs.keys()))
  #print 'VECTORIZED'
  X_out = v_func(X, *args, **kwargs)
  
  #print this_func + 'END'
  return X_out


# -------------------------------------------------------------------------------





##########################################################################################









# FIXME: MAKE IT GENERIC
def Find_Max(depths, R):
  """
  Find max values of each row of the matrix.
  Now used only in resolution matrix.
  
  """
  x_max, y_max = [], []
  
  glob_max = 0
  
  for row in np.arange(len(R)):
    if row != 0:
      row_max = max(R[row])
      no_layer = list(R[row]).index(row_max)
      if row_max > glob_max:
        glob_max = row_max
      #print depths[row], depths[no_layer]
      x_max.append(depths[row])
      y_max.append(depths[no_layer])
  return glob_max, x_max, y_max

# FIXME
def convert_2D_to_3D_array(Z):
  print('lib_generic/convert_2D_to_3D_array: START')
  
  mode = 1 
  
  if mode == 1:
    nx, nz = len(Z), len(Z[0])
    print(nx, nz)
    #quit()
    Z = np.array(Z)
    model = []
    for x in np.arange(nx):
      trace = []
      for z in np.arange(nz):
        trace.append(Z[x, z])
      model.append([trace]) # BRACKETS ARE CRUCIAL HERE
  
  elif mode == 2:
    print('for fwi_xyz2model')
    #  model = [] 
    #  for x in range(int(nx)):
    #    trace_x = []
    #    for y in range(int(ny)):
    #      trace_x.append([Z[x, y]]) # BRACES ARE CRUCIAL TO BE GENERIC
    #    model.append(trace_x)
    #  
    #  model = np.array(model) 
    #  return model
  
  else:
    eprint('lib_generic/convert_2D_to_3D_array: ERROR! Unknown mode\n')
    quit()
    
  print('lib_generic/convert_2D_to_3D_array: END')
  return np.array(model)
  
# EITHER X OR Y (NOT OBLIQUE SO FAR)
def prepare_2D_projection(data3D, x_or_y, coord):
  print('lib_generic/prepare_2D_projection: START')
  
  print('Assuming: data(x, y, z)')
  nx1, nx2, nx3 = len(data3D), len(data3D[0]), len(data3D[0][0])
  print(nx1, nx2, nx3)
  data3D = np.array(data3D)
  
  if x_or_y == 'x':
    data = data3D[coord]
    Y = list(range(nx2))
    Z = list(range(nx3))
    YY, ZZ = np.meshgrid(Y, Z)
    data = data.transpose()
    XX = YY * 0 + coord # SAME SHAPE AS YY, BUT ONLY ONE COORDINATE PICKED

  elif x_or_y == 'y':
    data = data3D[:, coord] # ONLY FOR ARRAYS (NOT LIST OF LISTS)
    X = list(range(nx1))
    Z = list(range(nx3))
    XX, ZZ = np.meshgrid(X, Z)
    #print len(XX), len(ZZ), len(data)
    #print 'len(XX), len(XX[0])', len(XX), len(XX[0])
    #print 'len(data), len(data[0])', len(data), len(data[0])
    data = data.transpose()
    #print 'len(data), len(data[0])', len(data), len(data[0])
    #quit()
    #data = np.reshape(data, ZZ.shape)
    YY = XX * 0 + coord # SAME SHAPE AS XX, BUT ONLY ONE COORDINATE PICKED
  
  else:
    eprint('lib_generic/prepare_2D_projection: ERROR! WRONG COORDINATE\n')
    quit()
  
  print('lib_generic/prepare_2D_projection: END')
  return XX, YY, ZZ, data


# FROM NUMBER OF MULTIDIM. LISTS PICK OUT ELEMENTS WITH x AS A FIRST INDEX
def cherrypick_x(x, *lists):
  nlists = []
  for l in lists:
    xy = l[x]
    nlists.append(xy)
  return nlists

# FROM NUMBER OF MULTIDIM. LISTS PICK OUT ELEMENTS WITH (x, y) AS 2 FIRST INDICES 
def cherrypick_xy(x, y, *lists):
  nlists = []
  for l in lists:
    xy = l[x][y]
    nlists.append(xy)
  return nlists



# ASSUMES FLOATS ONLY
def read_bin(file_name):
  from scipy.io import FortranFile
  f = FortranFile(file_name, 'r')
  print('File ', file_name, 'is about to be read.')
  
  content = []
  i_start = 0
  i = i_start
  
  while True:
    record = [] # EMPTY JUST IN CASE
    try:
      record = f.read_reals(dtype = np.float32)
    except TypeError:
      eprint('EOF reached.\n')
      break
    if i > i_start: # SKIP FIRST LINE WHICH IS ONLY 4-FLOATS LONG
      content.append(record)
    i += 1
  f.close()
  return content



#def read_wavefield_vtr(file_name):
  #from scipy.io import FortranFile
  #f = FortranFile(file_name, 'r')
  #print 'File ', file_name, 'is about to be read.'
  





# CHECK IF THE STRING STANDS FOR A FLOAT NUMBER AND RETURN TRUE/FALSE
def is_number(s):
  try:
    float(s)
    return True
  except ValueError:
    return False

# CHECK IF THE LIST CONSISTS OF INTs ONLY FIXME DOUBLE-CHECK IT
def str_is_int_list(l):
  l = np.ravel(l)
  for ll in l:
    if not is_int(ll):
      eprint('ERROR! ELEMENT OF THIS LIST IS NOT AN INT\n')
      quit()

# CHECK IF THE STRING STANDS FOR A INT NUMBER AND RETURN TRUE/FALSE
def str_is_int(s):
  try: 
    int(s)
    return True
  except ValueError:
    return False





# SEARCH THE LIST OF "NUMBERS" TO REPLACE NAN WITH ZEROS AND RETURN CORRECTED LIST
def NaN2zero(content):
  new_content = []
  for c in content:
    if isnan(c):
      c = 0.0
    new_content.append(c)
  return new_content

# SEARCH THE LIST OF "NUMBERS" TO REPLACE NAN WITH ZEROS AND RETURN CORRECTED LIST
def NaN2number(content, number):
  new_content = []
  for c in content:
    if isnan(c):
      c = number
    new_content.append(c)
  return new_content



# SEARCH THE LIST OF NUMBERS FROM THE END AND RETURN THE LENGTH OF ZERO-TAIL
def empty_ending(li):
  count = 0
  for i in reversed(li):
    if i == 0.:
      count += 1
    else:
      break
  return count

def Read_File_1st_char(file_name):
  f = open(file_name, 'r')    
  for line in f:
    line = line.split(None)
    first_char = line[0]
    break
  f.close()
  return first_char

def read_n_columns(file_name, n):
  f = open(file_name, 'r')
  content = []
  for line in f:
    content.append(line.split(None))
  f.close()  
  columns = []
  for i in np.arange(n):
    col = [line[i] for line in content]
    columns.append(col)
  return columns






def lengthen_matrix(M, xmin0, xmax0, xmin, xmax):
  print(xmin0, xmax0) # INITIAL LIMITS OF EACH MATRIX ROW
  print(xmin, xmax) # DESIRED LIMITS OF EACH MATRIX ROW (OBJECTIVE OF THIS FUNCTION)
  #print len(M[0]), len(M[-1])
  step = (xmax0 - xmin0) / len(M[0]) # CONSTANT DIFFERENCE BETWEEN SUBSEQUENT ELEMENTS OF THE ROW
  n_add = int(ceil((xmax - xmax0) / step)) # NUMBER OF ELEMENTS TO ADD AT THE END OF EACH ROW
  print('n_add: ', n_add)
  exact_xmax = xmax0 + n_add * step 
  #n_cut = floor((xmin - xmin0) / step) # NUMBER OF ELEMENTS TO SUBTRACT AT THE BEGINNING OF EACH ROW
  end = np.zeros(n_add) # NEW END OF EACH ROW
  new_M = []
  for row in M:
    new_M.append(np.concatenate((row, end))) # GLUE 2 NP.ARRAYS
  new_M = np.array(new_M) # MAKE NP.ARRAY OUT OF LIST
  
  return new_M, exact_xmax

# iTH COLUMN OF THE LIST
def ith_col(i, li):
  return [l[i] for l in li]



# PRINT UNKNOWN NUMBER OF ARGUMENTS  
def all_args(*arg):
  print("I was called with", len(arg), "arguments:", arg)


def replace_pattern(line, pattern_dict):
  print('Line before: ', line)
  for key in pattern_dict:
    print('Looking for a keyword: ', key)
    chunks = [_f for _f in line.split(key) if _f]
    if len(chunks) == 1:
      eprint('Line doesnt contain keyword ' + key + '\n')
    elif len(chunks) == 2:
      prefix, suffix = chunks
      line = prefix + pattern_dict[key] + suffix
    else:
      eprint('ERROR. LINE CONTAINED MORE THAN ONE OCCURENCE OF THE KEYWORD: ' + key + '\n')
      quit()
  print('Line after: ', line)
  return line



#
def my_formatter(x, pos):
  if x.is_integer():
    return str(int(x))
  else:
    return str(x)

#
def cm2inch(value):
  print('')
  #def cm2inch(*tupl):
  #inch = 2.54
  #if isinstance(tupl[0], tuple):
  #  return tuple(i/inch for i in tupl[0])
  #else:
  #  return tuple(i/inch for i in tupl)
  return value / 2.54
