"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permision writing to k.chrapkiewicz17@imperial.ac.uk.

This library provides framework procedures for...
At the moment fully implemented is:

1. ...

NOTE: Python devs usually don't write nested classes
and neither do I.

"""

## MODULES
import numpy as np
import matplotlib.pyplot as plt
import cmocean

## MY MODULES
from lib_generic import *
from lib_io_fullwave import Read_vtr, Save_vtr
from lib_generic_PLOTT import Plot, Animate
##

## CONSTS
this_lib = 'lib_fwi_project_CLASS.py/'

## CLASSES


# -------------------------------------------------------------------------------


class New_Class(object):
  """
  
  Parameters
  ----------
  
  Returns
  -------
  
  Notes
  -----
  
  """  
  
  # -----------------------------------------------------------------------------
  
  def __init__(self, **kwargs):
    """
    
    Parameters
    ----------
    
    Returns
    -------
    
    Notes
    -----
    
    """          
    raise NotImplementedError('This method needs to be overwitten in a subclass')
  
  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------
# 3 MAJOR PROJECT-TYPES: SYNTHETIC RUN, INVERSION & BOTH
# -------------------------------------------------------------------------------


class Proj(object):
  """
  Project: FWI or related.
  
  """
  
  this_class = 'Proj.'
  
  # -----------------------------------------------------------------------------
  
  def __init__(self, name, path, **kwargs):
    """
    Initialize a Fullwave project:
    - define the project
    - sets all variables defining i/o
    
    Parameters
    ----------
    
    Returns
    -------
    
    Notes
    -----
    It is Project_Suffices that picks the right 
    subset of Fullwave files.
    
    # NOTE: qp, qs ARE NOT RUNFILE PARAMS => THEY ARE PROBABLY 
    # SWITCHED ON IF TrueQp AND TrueQs FILES ARE PRESENT    
    """
    this_func = self.this_class + '__init__: '
    verbos = Kwarg('verbos', 1, kwargs)
    
    if path[-1] != '/': # ADD A TRAILING SLASH IF FORGOTTEN
      path += '/'
      
    self.path = path
    self.name = name
    self.proj_type = 'generic'
    
    proj_def = Proj_Def(self, **kwargs)
    paths = Paths(self, **kwargs)
    
    # MAKE DIRECTORIES
    if not Exists(path, **kwargs):
      o, e = Bash2('mkdir ' + path)
    if not Exists(path + '/inp/', **kwargs):
      o, e = Bash2('mkdir ' + path + '/inp/')
    if not Exists(path + '/out/', **kwargs):
      o, e = Bash2('mkdir ' + path + '/out/')
    
    try:
      # ENVIRONMENTAL VARIABLES
      self.env = Env_Vars(**kwargs) 
      self.geom = Geometry(self, **kwargs)
    except IOError:
      eprint(this_func + 'Make sure you have run inp.Prepare()\n')
      
      # THROUGHPUT    
    try:
      self.inp = Proj_Input(self, **kwargs)
      self.inp.Init(**kwargs)
    except IOError:
      eprint(this_func + 'Make sure you have run inp.Prepare()\n')
    
    try:
      self.out = Proj_Output(self, **kwargs)
      self.out.Init(**kwargs)
    except IOError:
      eprint(this_func + 'Make sure you have run inp.Prepare()\n')
      
    # LIST OF EXISTING FILES
    self.Fnames(**kwargs)  
    
  # -----------------------------------------------------------------------------

  def Duplicate(self, old_name, old_path, **kwargs):
    """
    Make the input of this project 
    the same as in old_proj.
    
    Parameters
    ----------
    
    Returns
    -------
    
    Notes
    -----
    # NOTE io OF OLD AN NEW ARE ASSUMED 
    TO BE THE SAME OR AT LEAST OVERLAP
    
    """
    #from lib_fwi_project import Project_Duplicate
    
    ##self.files2dupl = old_proj.inp.fnames_actual #NOTE
    #files2dupl = old_proj.inp.fnames_actual
    #Project_Duplicate(self, old_proj, files2dupl, **kwargs) 
    
    ## TRANSCRIBE MOST IMPORTANT ATTRIBUTES
    #self.dx = old_proj.dx
    #self.dims = old_proj.dims
    #self.nx = old_proj.nx
    #self.z_sea = old_proj.z_sea
    
    fnames = Get_Files(old_path, old_name + '*', **kwargs)
    for fname in fnames:
      print(fname)
    

  # -----------------------------------------------------------------------------
  
  def Fnames(self, **kwargs):
    """
    Get current list of files.
    
    Parameters
    ----------
    
    Returns
    -------
    
    Notes
    -----
    
    """
    fnames_inp = self.inp.Fnames(**kwargs)
    fnames_out = self.out.Fnames(**kwargs)
    
    #fnames_inp = Get_Files(self.path + '/inp/', self.name + '-*', **kwargs) # NOTE -, OTHERWISE p1==p11
    #fnames_out = Get_Files(self.path + '/out/', self.name + '-*', **kwargs)
    fnames = fnames_inp + fnames_out
    #fnames = [Path_Leave(fname) for fname in fnames]
    
    self.fnames_actual = fnames
    
    return fnames
    
  # -----------------------------------------------------------------------------
  
  def Ls(self, **kwargs):
    self.inp.Ls(**kwargs)
    self.out.Ls(**kwargs)

  # -----------------------------------------------------------------------------    
    
  def Tidy(self, **kwargs):
    self.inp.Tidy(**kwargs)
    self.out.Tidy(**kwargs)

  # -----------------------------------------------------------------------------  
  
  def Rm(self, **kwargs):
    raise NotImplementedError('Disabled (safety). Call inp.Rm & out.Rm instead.') 
  
  # -----------------------------------------------------------------------------
  # PRIVATE METHODS
  # -----------------------------------------------------------------------------  
  
  def _Init_Input(self, **kwargs):
    """
    Input files common for all projects.
    
    Notes
    -----
    Leading underscore - convention 
    for 'private' methods.
    
    """
    from lib_fwi_project import Project_Filenames
    
    if 'io' in kwargs:
      del kwargs['io']
    
    self.inp.files = Project_Filenames(self.name, self.io, 'synth',
                                       in_out='in', **kwargs)
    
    self.inp.fs = FS_File('FS', self, self.inp.path, **kwargs)
    self.inp.sr = SR_Files('S', 'R', self, self.inp.path, **kwargs)
    
    # SAMPLED IN TIME    
    self.inp.wavelet = Wavelet_File('RawSign', self, self.inp.path, **kwargs)
    
    #FIXME
    self.inp.rawdata = Raw_Data_File('RawSeis.sgy', self, self.inp.path, **kwargs)
    
    # META DATA
    self.inp.rawseis = RawSeis_File('RawSeis.txt', self, self.inp.path, **kwargs)
    self.inp.picks = Pick_Files(self, self.inp.path, **kwargs)
    
    # CONTROLS
    self.inp.sp = Segyprep_File('SegyPrep', self, self.inp.path, **kwargs)
    self.inp.runfile = Runfile('runfile', self, self.inp.path, **kwargs)
    self.inp.pbs = PBS_File('pbs', self, self.inp.path, **kwargs)

  # -----------------------------------------------------------------------------  

  def _Init_Output(self, **kwargs):
    pass

  # ----------------------------------------------------------------------------- 
  
  def _Prepare_Input(self, **kwargs):
    """
    
    """
    this_func = self.name + '.Prepare_Input: '
    
    verbos = Kwarg('verbos', 1, kwargs)
    
    #mods = Kwarg('mods', True, kwargs)
    wavelet = Kwarg('wavelet', True, kwargs)
    rawseis = Kwarg('rawseis', True, kwargs)
    sp = Kwarg('sp', True, kwargs)
    runfile = Kwarg('runfile', True, kwargs)
    
    rm = Kwarg('rm', True, kwargs)
    run = Kwarg('run', True, kwargs)
    #plot = Kwarg('plot', True, kwargs)
    cat = Kwarg('cat', True, kwargs)
    #anim = Kwarg('anim', False, kwargs) #NOTE
    check = Kwarg('check', True, kwargs)
    
    if rm:
      self.inp.Rm(**kwargs)
      self.out.Rm(**kwargs)
    
    # -----------------------------------------------------------------------------
    # SOURCE WAVELET
    # -----------------------------------------------------------------------------
    if wavelet:
      #print('\n\n' + 'Creating wavelet...' + '\n\n')      
      #self.inp.wavelet.Create(**kwargs)
      print(('\n\n' + 'Duplicating wavelet...' + '\n\n')) #FIXME: SPLIT INTO CASES   
      fname = '/home/kmc3817/heavy_PhD/source_wavelets/jm_source_Xd_2-6Hz_17_01_01.sgy'
      self.inp.wavelet.Duplicate(fname)
      #if plot > 0:
        #self.inp.wavelet.Plot(**kwargs)
   
   
    # -----------------------------------------------------------------------------
    # RAWSEIS TEXT FILE
    # -----------------------------------------------------------------------------    
    if rawseis:
      print(('\n\n' + 'Preparing RawSeis.txt...' + '\n\n'))      
      self.inp.rawseis.Create(**kwargs)
      if cat > 0:
        self.inp.rawseis.Cat(**kwargs)      

    # -----------------------------------------------------------------------------
    # SEGYPREP RUNFILE & RUN
    # -----------------------------------------------------------------------------    
    if sp:
      print('\n\n' + 'Preparing SegyPrep.key... (recipr=1)' + '\n\n')
      self.inp.sp.Create(reciprocity=1, **kwargs)
      if cat > 0:
        self.inp.sp.Cat(**kwargs)
      if run > 0:
        self.inp.sp.Run(**kwargs)

    
    # -----------------------------------------------------------------------------
    # SOURCES & RECEIVERS QC
    # -----------------------------------------------------------------------------          
    if cat > 0:
      self.inp.sr.Cat(**kwargs)
    #if plot > 0:
      #kwargs['dx'] = self.dx
      #self.inp.sr.Plot_All(**kwargs)
      #self.inp.sr.Show_Offsets(**kwargs)
    
    
    if runfile:
      print(('\n\n' + 'Preparing the runfile...' + '\n\n'))
      self.inp.runfile.Create(**kwargs)
      if cat > 0:
        self.inp.runfile.Cat(**kwargs)
     
    if verbos > 0:
      print('Note, PBS must be created manually')
    
    if check > 0:
      print('Checking all the input with Fullwave3D...')
      self.inp.Check()

  # -----------------------------------------------------------------------------
  
  def _Prepare_Output(self, **kwargs):
    pass

  # -----------------------------------------------------------------------------

  def _Check_Input(self, **kwargs):
    path = self.path + '/inp/'
    o, e = Bash2(self.paths['fullwave'] + ' -checkinput ' + self.name, path=path)
    print(o, e)

  # -----------------------------------------------------------------------------
  
  def _Check_Output(self, **kwargs):
    pass

  # -----------------------------------------------------------------------------
  
  def _Plot_Input(self, **kwargs):
    """
    
    """
    plot_wavelet = Kwarg('plot_wavelet', True, kwargs)
    plot_sr = Kwarg('plot_sr', True, kwargs)
    
    if plot_wavelet:
      self.inp.wavelet.Plot(**kwargs)
      
    #eprint('\n\nplot_sr disabled for now\n\n')
    if plot_sr:
      self.inp.sr.Plot(**kwargs)
      #self.inp.sr.Show_Offsets(dx=self.dx, **kwargs)

  # -----------------------------------------------------------------------------

  def _Plot_Output(self, **kwargs):
    pass

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


class Proj_Syn(Proj):
  """
  Generation of synthetics.
  
  
  
  """
  this_class = 'Proj_Syn.'
  
  # -----------------------------------------------------------------------------  
  
  def __init__(self, name, path, **kwargs):
    """
    
    """
    this_func = self.this_class + '__init__: '
    
    super(Proj_Syn, self).__init__(name, path, **kwargs)
    
    self.problem = 'synthetic' # OVERWRITE THE DEFAULT
    self.proj_type = 'synth'
    #print 'synth'
    
  # -----------------------------------------------------------------------------
  # PRIVATE METHODS
  # -----------------------------------------------------------------------------  
  
  def _Init_Input(self, **kwargs):
    """
    
    Notes
    -----
    Leading underscore - convention 
    for 'private' methods.
    
    self.inp.vp.true = self.vp.true #NOTE: ADD THESE OR NOT?    
    
    """
    super(Proj_Syn, self)._Init_Input(**kwargs)
    
    # MODELS
    
    # THIS SERVES AS STH EASY TO ITERATE OVER (FOR PLOTS ETC.)
    self.mods_true = [] # NEEDS TO BE DIFFERENT THAN mods_it BECAUSE DOESN'T HAVE .it
    
    self.vp = Model_Files(self, 'Vp', **kwargs)
    self.vp.true = Model_File('TrueVp', self, self.inp.path, **kwargs)
    self.mods_true.append(self.vp.true)
    
    if self.anisotropy != 'none':
      self.delta = Model_Files(self, 'Delta', **kwargs)
      self.delta.true = Model_File('TrueDelta', self, self.inp.path, **kwargs)
      self.mods_true.append(self.delta.true)
      
      self.epsil = Model_Files(self, 'Epsilon', **kwargs)
      self.epsil.true = Model_File('TrueEpsil', self, self.inp.path, **kwargs)
      self.mods_true.append(self.epsil.true)
      
    if self.qp:
      self.qp = Model_Files(self, 'Qp', **kwargs)
      self.qp.true = Model_File('TrueQp', self, self.inp.path, **kwargs)
      self.mods_true.append(self.qp.true)
      
    if self.qs:
      self.qs = Model_Files(self, 'Qs', **kwargs)
      self.qs.true = Model_File('TrueQs', self, self.inp.path, **kwargs)
      self.mods_true.append(self.qs.true)
      
    if self.equation == 'elastic':
      self.vs = Model_Files(self, 'Vs', **kwargs)
      self.vs.true = Model_File('TrueVs', self, self.inp.path, **kwargs)
      self.mods_true.append(self.vs.true)
    
    for mod in self.mods_true:  
      mod.Convert(**kwargs)
    
  # -----------------------------------------------------------------------------  

  def _Init_Output(self, **kwargs):
    """
    This should be called by
    proj.out.Plot(...)

    Notes
    -----
    Leading underscore - convention 
    for 'private' methods.
    
    """
    super(Proj_Syn, self)._Init_Output(**kwargs)
    
    self.out.synout = Log_File('SynOut', self, self.out.path, **kwargs)
    self.out.synerr = Log_File('SynErr', self, self.out.path, **kwargs) 
 
    # OUTPUT SPECIFIC FOR THIS PROJECT TYPE
    self.out.synth = Data_File('Synth', self, self.out.path, **kwargs)
    self.out.synth_idx = Data_File('SynthIdx', self, self.out.path, **kwargs)
    #self.out.synth.Files(**kwargs) # IN Prepare_Output, LIKE fw
    
    
    snap_step = -500 #FIXME: READ IT FROM SOMEWHERE
    self.out.fw = Wavefield_Files('fw', snap_step, self, self.out.path, **kwargs)

  # ----------------------------------------------------------------------------- 
  
  def _Prepare_Input(self, **kwargs):
    """
    
    """
    prep_mods = Kwarg('prep_mods', True, kwargs)
    rm = Kwarg('rm', True, kwargs)
    
    if rm: # FIXME: ALREADY DONE SOMEWHERE?
      self.inp.Rm(**kwargs)
      self.out.Rm(**kwargs)
      
    # -----------------------------------------------------------------------------
    # TRUE MODELS
    # -----------------------------------------------------------------------------
    
    if prep_mods:
      for mod in self.mods_true:
        print(('\n\n' + 'Preparing ' + mod.name + '...' + '\n\n')) #FIXME: VP SPECIFIC
        file2dupl = '/home/kmc3817/heavy_PhD/start_mods/Ben_whole_model_24-04-18_x-30000_30000_y-5000_15000.sgy'
        #NOTE: IT HAS z_origin AND dims AS BELOW
        
        mod.Duplicate(file2dupl, dupl_cmd='cp', 
                               z_origin=-1500.0, dims=[1201, 401, 131], 
                               **kwargs)
        
        mod.Resize(**kwargs)
        mod.Convert(**kwargs)
      
        #if plot > 0:
          #self.vp.true.Plot(**kwargs)
        #if anim > 0:
          #self.vp.true.Animate(**kwargs)    

    # -----------------------------------------------------------------------------
    # REMAINDER
    # -----------------------------------------------------------------------------    
    
    # NOTE: AFTER MODELS BECAUSE SP NEEDS THEM
    kwargs['rm'] = 0 #NOTE: IMPORTANT
    super(Proj_Syn, self)._Prepare_Input(**kwargs)
   
  # -----------------------------------------------------------------------------
  
  def _Prepare_Output(self, **kwargs):
    super(Proj_Syn, self)._Prepare_Output(**kwargs)
    
    kwargs['convert'] = True
    self.out.synth.Split(**kwargs)
    
    self.out.fw.Files(**kwargs)

  # -----------------------------------------------------------------------------

  def _Check_Input(self, **kwargs):
    super(Proj_Syn, self)._Check_Input(**kwargs)    

  # -----------------------------------------------------------------------------
  
  def _Check_Output(self, **kwargs):
    super(Proj_Syn, self)._Check_Output(**kwargs)

  # -----------------------------------------------------------------------------
  
  def _Plot_Input(self, **kwargs):
    """
    This should be called by
    proj.out.Plot(...)

    Notes
    -----
    Leading underscore - convention 
    for 'private' methods.
    
    """
    this_func = self.this_class + '_Plot_Input: '    
    verbos = Kwarg('verbos', 1, kwargs)
    if verbos == verbos_func:
      print(this_func + 'START')
    
    super(Proj_Syn, self)._Plot_Input(**kwargs)
    
    plot_mods = Kwarg('plot_mods', True, kwargs)
    # COMMON ERROR (OVERWRITTING OF A kwargs)
    # NOTE: MUST BE AFTER super(
    if 'plot_sr' in kwargs:
      del kwargs['plot_sr']
    
    if plot_mods:
      for mod in self.mods_true:
        mod.Plot_Various_Slices(**kwargs)
    
    if verbos == verbos_func:
      print(this_func + 'END')

  # -----------------------------------------------------------------------------

  def _Plot_Output(self, **kwargs):
    """
    This should be called by
    proj.out.Plot(...)

    Notes
    -----
    Leading underscore - convention 
    for 'private' methods.
    
    """
    """
    This should be called by
    proj.out.Plot(...)
    
    Noteskwargs['
    -----
    Leading underscore - convention 
    for 'private' methods.
    
    """
    super(Proj_Syn, self)._Plot_Output(**kwargs)
    
    plot_fw = Kwarg('plot_fw', True, kwargs)
    plot_synth = Kwarg('plot_synth', True, kwargs)    
    
    kwargs['save'] = Kwarg('save', 1, kwargs)
    
    if plot_fw:
      kwargs['animate'] = Kwarg('animate', 1, kwargs)
      
      kwargs['plot_mod'] = 1
      self.out.fw.Plot(**kwargs)
      kwargs['plot_mod'] = 0
      self.out.fw.Plot(**kwargs)

    if plot_synth:
      kwargs['juxtapose'] = 0
      self.out.synth.Plot(**kwargs)

  # -----------------------------------------------------------------------------
    

# -------------------------------------------------------------------------------


class Proj_Inv(Proj):
  """
  Inversion.
  
  """

  # -----------------------------------------------------------------------------  
  
  def __init__(self, name, path, **kwargs):
    """
    NOTE: CHECKING IF OBJECTS EXIST FOR Proj_Both
    (IT MIGHT HAVE BEEN CREATED BY Proj_Syn)
  
    """
    super(Proj_Inv, self).__init__(name, path, **kwargs)
    self.problem = 'tomography' # OVERWRITE THE DEFAULT
    self.proj_type = 'inv'
    ####print 'inv'
    
  # -----------------------------------------------------------------------------
  # PRIVATE METHODS
  # -----------------------------------------------------------------------------  
  
  def _Init_Input(self, **kwargs):
    """
    
    Notes
    -----
    Leading underscore - convention 
    for 'private' methods.
    
    """
    super(Proj_Inv, self)._Init_Input(**kwargs)
    
    if 'dump_fw' in self.env.var:
      del self.env.var['dump_fw'] # FIXME
    
    # PREPARE A LIST USEFUL FOR Set_Iterations (append IN IFs)
    self.mods_start = []
    self.mods_it = []
    
    #if not self.vp: 
    self.vp = Model_Files(self, 'Vp', **kwargs)
    self.vp.start = Model_File('StartVp', self, self.inp.path, **kwargs)
    self.mods_it.append(self.vp)
    self.mods_start.append(self.vp.start)
    
    
    if self.anisotropy != 'none':
      #if not self.delta: 
      self.delta = Model_Files(self, 'Delta', **kwargs)      
      self.delta.start = Model_File('StartDelta', self, self.inp.path, **kwargs)
      self.mods_it.append(self.delta)
      self.mods_start.append(self.delta.start)

      #if not self.epsil: 
      self.epsil = Model_Files(self, 'Epsilon', **kwargs)      
      self.epsil.start = Model_File('StartEpsil', self, self.inp.path, **kwargs)
      self.mods_it.append(self.epsil)
      self.mods_start.append(self.epsil.start)
      
    if self.qp: # FIXME????!!!
      #if not self.qp: # FIXME????!!!
      self.qp = Model_Files(self, 'Qp', **kwargs)      
      self.qp.start = Model_File('StartQp', self, self.inp.path, **kwargs)
      self.mods_it.append(self.qp)
      self.mods_start.append(self.qp.start)
      
    if self.qs:
      #if not self.qs: 
      self.qs = Model_Files(self, 'Qs', **kwargs)      
      self.qs.start = Model_File('StartQs', self, self.inp.path, **kwargs)
      self.mods_it.append(self.qs)
      self.mods_start.append(self.qs.start)
      
    if self.equation == 'elastic':
      #if not self.vs: 
      self.vs = Model_Files(self, 'Vs', **kwargs)      
      self.vs.start = Model_File('StartVs', self, self.inp.path, **kwargs)
      self.mods_it.append(self.vs)   
      self.mods_start.append(self.vs.start)
    
    # CONVERT ALL
    for mod in self.mods_start: 
      mod.Convert(**kwargs)     
    
    
    # INIT DERIVATIVES 
    if self.env.var['dump_grad']:
      self.grad = Model_Files(self, 'Grad', **kwargs)
      self.mods_it.append(self.grad) 

    if self.env.var['dump_rawgrad']:
      self.rawgrad = Model_Files(self, 'RawGrad', **kwargs)
      self.mods_it.append(self.rawgrad)

    if self.env.var['dump_prec']:
      self.prec = Model_Files(self, 'Prec', **kwargs)
      self.mods_it.append(self.prec)     

    if self.env.var['dump_rawprec']:
      self.rawprec = Model_Files(self, 'RawPrec', **kwargs)
      self.mods_it.append(self.rawprec)
    
    
    # INIT OBSERVED DATA ?
    self.inp.obser = Data_File('Obser', self, self.inp.path, **kwargs)
    self.inp.obser_idx = Data_File('ObserIdx', self, self.inp.path, **kwargs)
    
    
  # -----------------------------------------------------------------------------  

  def _Init_Output(self, **kwargs): 
    """
    This should be called by
    proj.out.Plot(...)

    Notes
    -----
    Leading underscore - convention 
    for 'private' methods.
    
    FORMER SET_ITERATIONS!!!!
    
    """
    verbos = Kwarg('verbos', 1, kwargs)
    
    super(Proj_Inv, self)._Init_Output(**kwargs)
    
    self.out.invout = Log_File('InvOut', self, self.out.path, **kwargs)
    self.out.inverr = Log_File('InvErr', self, self.out.path, **kwargs)    
    
    blocks = Kwarg('blocks', [], kwargs) #FOR RUNFILE
    self.blocks = blocks
    
    niters = 0
    for block in blocks:
      niters += block.niters
    
    self.niters = niters # REDUNDANCY
    if verbos > 0:
      print(('Proj_Inv._Init_Output: Total no. of inversion ' + 
            'iterations: ' + str(niters)))
    
    self.out.functional = Functional(self, **kwargs)
    
    for model in self.mods_it: # model = vp ETC.
      
      ID = model.ID
      suffix = ID
      
      if ('Grad' in ID) or ('Prec' in ID): # ANY MORE?
        ext = 'vtr' # THAT'S WHAT FULLWAVE'S DUMP FORMAT IS
        model.it = list(np.zeros(niters))
        i_min = 0
        i_add = 1        
        
      else:
        ext = Ext(model.start.name)
        model.it = list(np.zeros(niters + 1))
        model.it[0] = model.start
        i_min = 1
        i_add = 0
        
      for i in range(i_min, len(model.it)):
        
        prefix = self.name + '-CP' + str(i+i_add).rjust(5, '0')
        fname = prefix + '-' + suffix + '.' + ext
        
        if 'Grad' in ID:
          model.it[i] = Grad_File(fname, self, self.out.path, **kwargs) 
        elif 'Prec' in ID:
          model.it[i] = Prec_File(fname, self, self.out.path, **kwargs) 
        else:
          model.it[i] = Model_File_CP(fname, self, self.out.path, **kwargs) 
  
  # ----------------------------------------------------------------------------- 
  
  def _Prepare_Input(self, **kwargs):
    """
    Very similar to Proj_Syn (true -> start)
    
    """
    rm = Kwarg('rm', True, kwargs)
    
    if rm:
      self.inp.Rm(**kwargs)
      self.out.Rm(**kwargs)    
    
    prep_mods = Kwarg('prep_mods', True, kwargs)
    
    # -----------------------------------------------------------------------------
    # START MODELS
    # -----------------------------------------------------------------------------
    
    if prep_mods:
      for mod in self.mods_start:
        print(('\n\n' + 'Preparing ' + mod.name + '...' + '\n\n')) #FIXME: VP SPECIFIC
        file2dupl = '/home/kmc3817/heavy_PhD/start_mods/Ben_whole_model_24-04-18_x-30000_30000_y-5000_15000.sgy'
        #NOTE: IT HAS z_origin AND dims AS BELOW
        
        mod.Duplicate(file2dupl, dupl_cmd='cp', 
                               z_origin=-1500.0, dims=[1201, 401, 131], 
                               **kwargs)
        mod.Convert(**kwargs)
        mod.Resize(**kwargs)
      
        #if plot > 0:
          #self.vp.true.Plot(**kwargs)
        #if anim > 0:
          #self.vp.true.Animate(**kwargs)    


    # -----------------------------------------------------------------------------
    # SPLIT OBSERVED DATA
    # -----------------------------------------------------------------------------
    self.inp.obser.Split(**kwargs) #FIXME: HERE?
     
     
    # -----------------------------------------------------------------------------
    # REMAINDER
    # -----------------------------------------------------------------------------    
    
    # NOTE: AFTER MODELS BECAUSE SP NEEDS THEM
    kwargs['rm'] = 0 #NOTE: IMPORTANT
    super(Proj_Inv, self)._Prepare_Input(**kwargs)    

  # -----------------------------------------------------------------------------
  
  def _Prepare_Output(self, **kwargs):
    """
    
    """
    super(Proj_Inv, self)._Prepare_Output(**kwargs)  
    
    prep_functional = Kwarg('prep_functional', True, kwargs)
    prep_mods = Kwarg('prep_mods', True, kwargs)
    
    if prep_functional:
      eprint('Note: is timestamp correct?\n')
      self.out.functional.Read(**kwargs)
    
    if prep_mods:
      for mod_it in self.mods_it:
        mod_it.Convert(**kwargs)
  
  # -----------------------------------------------------------------------------

  def _Check_Input(self, **kwargs):
    super(Proj_Inv, self)._Check_Input(**kwargs)    

  # -----------------------------------------------------------------------------
  
  def _Check_Output(self, **kwargs):
    super(Proj_Inv, self)._Check_Output(**kwargs)

  # -----------------------------------------------------------------------------
  
  def _Plot_Input(self, **kwargs):
    """

    Notes
    -----
    Leading underscore - convention 
    for 'private' methods.
    
    """
    plot_mods = Kwarg('plot_mods', True, kwargs)
    
    super(Proj_Inv, self)._Plot_Input(**kwargs)
    
    # COMMON ERROR (OVERWRITTING OF A kwargs)
    # NOTE: MUST BE AFTER super(
    if 'plot_sr' in kwargs:
      del kwargs['plot_sr']
    
    if plot_mods:
      for mod in self.mods_start:
        mod.Plot_Various_Slices(**kwargs)
        
    self.inp.obser.Plot(**kwargs)

  # -----------------------------------------------------------------------------

  def _Plot_Output(self, **kwargs):
    """
    This should be called by
    proj.out.Plot(...)

    Notes
    -----
    Leading underscore - convention 
    for 'private' methods.
    
    FIXME: FOR MULTI-PARAMETER FWI - ARE THERE 
    MORE GRADIENT FILES?
    
    """
    super(Proj_Inv, self)._Plot_Output(**kwargs)
    
    plot_functional = Kwarg('plot_functional', True, kwargs)
    plot_mod_it = Kwarg('plot_mod_it', True, kwargs)
    recover = Kwarg('recover', True, kwargs)
    recover_only = Kwarg('recover_only', False, kwargs)
    
    if plot_functional:
      self.out.functional.Plot(**kwargs)  
    
    if plot_mod_it:
      for mod in self.mods_it: # THIS INCLUDES GRAD, PREC, ETC.
        ID = mod.ID
        
        for model in mod.it:
          # GUARD AGAINST MISSING FILES
          if not Exists(model.fname, **kwargs):
            eprint('File ' + model.fname + 
                   ' not found. Skipping this figure.\n')          
            continue
          
          # PLOT THE MODEL
          if not recover_only:
            model.Plot_Various_Slices(**kwargs)
          
          # PLOT THE MODEL MINUS THE START MODEL 
          if recover and not (('Grad' in ID) or ('Prec' in ID)):
            model.Plot_Recovered_Anom(**kwargs)
          
          #break
        #break
      #eprint('brrreaaak\n')
      
  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


class Proj_Both(Proj_Inv, Proj_Syn): 
  """
  NOTE: SINCE IT IS NOW JUST MADE 
  TO BE A CONTAINER OF 2 SUBPROJECTS
  (syn & inv), ALL OF THE MULTI-INHERITANCE
  DOESN'T MATTER AND IT ACTUALLY SHOULDN'T 
  BE CALLED. BUT I AM NOT DELETING IT
  BECAUSE I MIGHT REVERT TO IT IN THE FUTURE.
  
  Inversion of synthetics.
  
  Notes
  -----
  Multiple inheritance - for a child class:
    
    class C(B, A):
       ...
  
  Python will start by looking at A, and, if A doesn't have
  the attribute, then it will look at B.
  
  NOTE: Be mindful of 'diamond inheritance'  problem!
  
  """
  
  # -----------------------------------------------------------------------------
  
  def __init__(self, name, path='.', **kwargs):
    """
    """
    #super(Proj_Both, self).__init__(name, path, **kwargs) # NEEDS TO STAY
    self.problem = 'both' # OVERWRITE THE DEFAULT
    self.proj_type = 'both'

  # -----------------------------------------------------------------------------
  
  def Ls(self, **kwargs):
    self.syn.Ls(**kwargs)
    self.inv.Ls(**kwargs)

  # -----------------------------------------------------------------------------

  def Tidy(self, **kwargs):
    self.syn.Tidy(**kwargs)
    self.inv.Tidy(**kwargs)
    
  # -----------------------------------------------------------------------------
  # PRIVATE METHODS
  # -----------------------------------------------------------------------------  
  
  def _Init_Input(self, **kwargs):
    pass
    #super(Proj_Both, self)._Init_Input(**kwargs)
    
  # -----------------------------------------------------------------------------  

  def _Init_Output(self, **kwargs):
    pass
    #super(Proj_Both, self)._Init_Output(**kwargs)

  # ----------------------------------------------------------------------------- 
  
  def _Prepare_Input(self, **kwargs):
    #super(Proj_Both, self)._Prepare_Input(**kwargs)
    pass

    
  # -----------------------------------------------------------------------------
  
  def _Prepare_Output(self, **kwargs):
    pass
    #super(Proj_Both, self)._Prepare_Output(**kwargs)

  # -----------------------------------------------------------------------------

  def _Check_Input(self, **kwargs):
    super(Proj_Both, self)._Check_Input(**kwargs)

  # -----------------------------------------------------------------------------
  
  def _Check_Output(self, **kwargs):
    pass
    #super(Proj_Both, self)._Check_Output(**kwargs)

  # -----------------------------------------------------------------------------
  
  def _Plot_Input(self, **kwargs):
    pass
    #super(Proj_Both, self)._Plot_Input(**kwargs)

  # -----------------------------------------------------------------------------

  def _Plot_Output(self, **kwargs):
    pass
    #super(Proj_Both, self)._Plot_Output(**kwargs)

  # -----------------------------------------------------------------------------

          
# -------------------------------------------------------------------------------
# DIFFERENT SUBTYPES OF ABOVE PROJECT-TYPES
# -------------------------------------------------------------------------------


class Proj_Syn_vs_Obs(Proj_Syn):
  """
  Generation of synthetics and 
  juxtaposition against observed data
  to track cycle-skipping and general
  QC.
  
  """
  this_class = 'Proj_Syn_vs_Obs.'
  
  # -----------------------------------------------------------------------------   
  
  def __init__(self, name, path, **kwargs):
    super(Proj_Syn_vs_Obs, self).__init__(name, path, **kwargs)
    self.proj_type = 'syn_vs_obs'
  
  # -----------------------------------------------------------------------------
  # PRIVATE METHODS
  # -----------------------------------------------------------------------------  
  
  def _Init_Input(self, **kwargs):
    super(Proj_Syn_vs_Obs, self)._Init_Input(**kwargs)
    self.inp.outseis = Data_File('OutSeis', self, self.inp.path, 
                                 against='synth', **kwargs)
    
  # -----------------------------------------------------------------------------  

  def _Init_Output(self, **kwargs):
    super(Proj_Syn_vs_Obs, self)._Init_Output(**kwargs)

  # ----------------------------------------------------------------------------- 
  
  def _Prepare_Input(self, **kwargs):
    """
    
    """
    super(Proj_Syn_vs_Obs, self)._Prepare_Input(**kwargs)

    split = Kwarg('split', True, kwargs)
    check = Kwarg('check', True, kwargs) 
    
    kwargs['convert'] = Kwarg('convert', True, kwargs)
    
    if split > 0:
      self.inp.outseis.Split(**kwargs) 

  # -----------------------------------------------------------------------------
  
  def _Prepare_Output(self, **kwargs):
    super(Proj_Syn_vs_Obs, self)._Prepare_Output(**kwargs)

  # -----------------------------------------------------------------------------

  def _Check_Input(self, **kwargs):
    super(Proj_Syn_vs_Obs, self)._Check_Input(**kwargs)

  # -----------------------------------------------------------------------------
  
  def _Check_Output(self, **kwargs):
    super(Proj_Syn_vs_Obs, self)._Check_Output(**kwargs)

  # -----------------------------------------------------------------------------
  
  def _Plot_Input(self, **kwargs):
    super(Proj_Syn_vs_Obs, self)._Plot_Input(**kwargs)
    
    kwargs['juxtapose'] = Kwarg('juxtapose', 0, kwargs)
    kwargs['plot_wiggles'] = Kwarg('plot_wiggles', 0, kwargs) #FIXME
    kwargs['save'] = Kwarg('save', 1, kwargs)
    
    self.inp.outseis.Plot(**kwargs)
  
  # -----------------------------------------------------------------------------

  def _Plot_Output(self, **kwargs):
    kwargs['plot_synth'] = 0 # IT WILL BE PLOTTED BY self.out.synth.Plot(**kwargs)
    super(Proj_Syn_vs_Obs, self)._Plot_Output( **kwargs)

    kwargs['juxtapose'] = Kwarg('juxtapose', 1, kwargs)
    kwargs['plot_wiggles'] = Kwarg('plot_wiggles', 0, kwargs) #FIXME
    kwargs['save'] = Kwarg('save', 1, kwargs)
    kwargs['plot_synth'] = 1

    self.out.synth.Plot(**kwargs)

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


class Proj_Inv_Synth(Proj_Both): 
  """
  Inversion of synthetic data, e.g. 
  checkerboard study.
  
  # CAUTION
  SEE PROJ_BOTH NOTES
  
  """
  
  this_class = 'Proj_Inv_Synth.'
  
  # -----------------------------------------------------------------------------
  
  def __init__(self, name, path, **kwargs):
    """
    
    """
    this_func = self.this_class + '__init__: '
    
    #super(Proj_Inv_Synth, self).__init__(name, path, **kwargs)
    self.proj_type = 'inv_synth'
    
    # MUST CREATE BEFORE INIT. Syn AND Inv
    if not Exists(path, **kwargs):
      o, e = Bash2('mkdir ' + path)
    if not Exists(path + '/syn/', **kwargs):
      o, e = Bash2('mkdir ' + path + '/syn/')
    if not Exists(path + '/inv/', **kwargs):
      o, e = Bash2('mkdir ' + path + '/inv/')
    
    self.syn = Proj_Syn_vs_Obs(name + 's', path=path + '/syn/', **kwargs)
    kwargs['plot'] = 0 # PLOT ONLY ONCE (IF ANY)
    self.inv = Proj_Inv(name + 'i', path=path + '/inv/', **kwargs)
      
    # THIS ALLOWS SOME SUPERVISION OVER SUBPROJECTS
    #try:
    self._Init_Input(**kwargs)
    self._Init_Output(**kwargs)
    #except AttributeError:
      #eprint(this_func + 'Make sure you have run inp.Prepare()\n')

  # -----------------------------------------------------------------------------
  # PRIVATE METHODS
  # -----------------------------------------------------------------------------  
  
  def _Init_Input(self, **kwargs):
    """
    
    """
    # WE CAN'T CALL THIS BECAUSEE OF THE SPECIAL CONSTRUCTOR
    #super(Proj_Inv_Synth, self)._Init_Input(**kwargs)
    
    self.syn._Init_Input(**kwargs)
    self.inv._Init_Input(**kwargs)
    
    
    # SYNERGISTIC BIT
    for mod in self.syn.mods_true:
      fname_anom = Strip(mod.name) + '_Anom.vtr' #NOTE: ONLY vtr
      mod.anom = Anomaly_File(fname_anom, self.syn, 
                              self.syn.inp.path, **kwargs)
    
  # -----------------------------------------------------------------------------  

  def _Init_Output(self, **kwargs):
    # WE CAN'T CALL THIS BECAUSEE OF THE SPECIAL CONSTRUCTOR    
    #super(Proj_Inv_Synth, self)._Init_Output(**kwargs)
    self.syn.out.Init(**kwargs)
    self.inv.out.Init(**kwargs)

  # ----------------------------------------------------------------------------- 
  
  def Prepare_Input(self, **kwargs): # NOTE NOT PRIVATE
    """
    
    """
    verbos = Kwarg('verbos', 1, kwargs)
    
    prep_syn = Kwarg('prep_syn', True, kwargs)
    prep_inv = Kwarg('prep_inv', True, kwargs)
    prep_anom = Kwarg('prep_anom', True, kwargs)

    if prep_syn:
      self.syn.inp.Prepare(**kwargs)

    if prep_inv:
      self.inv.inp.Prepare(**kwargs)
    
    # SYNERGISTIC BIT
    if prep_anom:
      #FIXME THIS ASSUMES THE SAME ORDER OF APPENDING 
      for i, m in enumerate(self.syn.mods_true):
        print(m.name)
        m.Create_Anomaly(
          shape='rect', 
          sizes=np.array([10, 10, 10]),
          ampl=0.05,
          sparsity='dense',
          pads=[[5,5], [5,5], [20,5]], 
          **kwargs)            
        m.Add_Anomaly(**kwargs)
        
        # MOVING TO INV
        #path = self.syn.inp.path
        #source = self.syn.name + '-' + self.inv.mods_start[i].suffix
        #destination = '../../inv/inp/' + self.inv.mods_start[i].name
        #if verbos > 0:
        #  print('Moving ' + source + ' from path ' + path + ' to ' +
        #        destination)
        #
        #o, e = Bash2('mv ' + source + ' ' + destination, path=path,
        #             **kwargs)
        #
        ## AND CONVERTING AT DESTINATION #FIXME: WILL NOT OVERWRITE
        #self.inv.mods_start[i].Convert()
    
    #  if plot > 0:
    #    xyz = [25,25,20]
    #    self.vp.start.Plot(xyz=xyz, **kwargs)
    #    self.vp.true.anom.Plot(xyz=xyz, **kwargs)
    #    self.vp.true.Plot(xyz=xyz, **kwargs)
        
      
    # self.inv.inp.Prepare(**kwargs)
      
      
    #NOTE: READ ttr INSTEAD!!! OTHERWISE TOO EXPENSIVE
    #if subproj:
    #  # QC OF THE STARTMOD (CYCLE-SKIPPING) - SYNTHETICS
    #  it = 0
    #  self.out.Init_Subproject(it, **kwargs)
    #       
    #  if check > 0:
    #    self.out.synCP.inp.Check(**kwargs)
    #  
    #  
    #  if run > 0:
    #    self.out.Create_Subproject(it, **kwargs) 
    #    print 'Running synthetics from the start model...'
    #    
    #    o, e = Bash2('fwi_run_locally.sh synth 690 ' +
    #                 self.out.synCP.name  + ' 8', **kwargs)
    #  if split > 0:
    #    self.out.synCP.out.synth.Split(**kwargs)
    #    #self.out.synCP.out.synth.Files()
    #  
    #  if plot > 0:
    #    pass
    #    #self.out.synCP.out.synth.split['4136']['27'].Plot()
    #  
    #  
    #  # 'OBSERVED' DATA (SYNTHETICS FROM THE TRUE MODEL)
    #  if run > 0:
    #    print 'Running synthetics from the true model...'
    #    o, e = Bash2('fwi_run_locally.sh synth 690 ' + self.name  + ' 8', **kwargs)
    #  
    #  if split > 0:
    #    self.out.synth.Split(**kwargs)
    #    # self.out.synth.Files()
    #  if plot > 0:
    #    pass
    #    #obs = '4136'
    #    #line = '8'
    #    #proj7.out.synth.split[obs][line].Plot()
    #    #proj7.out.synth.split[obs][line].Juxtapose('CP', juxt='diff', yflip=1)  

  # -----------------------------------------------------------------------------
  
  def Prepare_Output(self, **kwargs):
    """
    
    """
    prep_syn = Kwarg('prep_syn', True, kwargs)
    prep_inv = Kwarg('prep_inv', True, kwargs)
    
    synth = self.syn.out.synth.fname
    synth_idx = self.syn.out.synth_idx.fname
    obser = self.inv.inp.obser.fname
    obser_idx = self.inv.inp.obser_idx.fname
    
    if prep_syn:
      self.syn.out.Prepare(**kwargs)
  
      # SYNERGISTIC BIT
      if not (Exists(synth) and Exists(synth_idx)):
        raise IOError('Either ' + synth + ' or ' +
                      synth_idx + ' does not exist.')
      
      sources = [synth, synth_idx]
      destins = [obser, obser_idx]
      
      for source, destin in zip(sources, destins):
        o, e = Bash2('cp ' + source + ' ' + destin) # IT WILL OVERWRITE

    if 0: #prep_inv: #FIXME
      self.inv.out.Prepare(**kwargs)
      # SYNERGISTIC BIT
      self.vp.Recover_Anomaly()

  # -----------------------------------------------------------------------------

  def Check_Input(self, **kwargs):
    self.syn.inp.Check(**kwargs)
    self.inv.inp.Check(**kwargs)

  # -----------------------------------------------------------------------------
  
  def _Check_Output(self, **kwargs):
    super(Proj_Inv_Synth, self)._Check_Output(**kwargs)

  # -----------------------------------------------------------------------------
  
  def Plot_Input(self, **kwargs):
    verbos = Kwarg('verbos', 1, kwargs)
    
    plot_syn = Kwarg('plot_syn', True, kwargs)
    plot_inv = Kwarg('plot_inv', True, kwargs)
    plot_anom = Kwarg('plot_anom', True, kwargs)

    if plot_syn:
      self.syn.inp.Plot(**kwargs)

    if plot_inv:
      self.inv.inp.Plot(**kwargs)
      
    if plot_anom:
      for mod in self.syn.mods_true:
        print('Plotting ' + mod.anom.fname) 
        kwargs['save'] = Kwarg('save', 1, kwargs)
        mod.anom.Plot_Various_Slices(**kwargs)

  # -----------------------------------------------------------------------------

  def Plot_Output(self, **kwargs):
    plot_syn = Kwarg('plot_syn', True, kwargs)
    plot_inv = Kwarg('plot_inv', True, kwargs)
    #plot_anom = Kwarg('plot_anom', True, kwargs)    
    
    if plot_syn:
      self.syn.out.Plot(**kwargs)    
    
    if plot_inv:
      self.inv.out.Plot(**kwargs)    
    

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# PROJECT INPUT/OUTPUT
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


class Proj_Throughput(object):
  """
  Parent class of input and output.
  
  Notes
  -----
  Not much implemented - it's left 
  for subclasses to define their own
  versions of methods.
  
  """  
  
  # -----------------------------------------------------------------------------
  
  def __init__(self, proj, **kwargs):
    """
    Inherit most frequently used 
    project-attributes.
    
    Parameters
    ----------
    
    Returns
    -------
    
    Notes
    -----
    
    """          
    self.proj = proj # ALL INFO!
    #NOTE: WHAT'S BELOW IS REDUNDANT
    self.name = proj.name
    self.pname = proj.name #NOTE: REDUNDANCY
    self.io = proj.io 
    
  # -----------------------------------------------------------------------------
  
  def Fnames(self, **kwargs):
    """
    List existing files with proj_name prefix. 
    
    # NOTE -, OTHERWISE p1==p11
    """
    fnames = Get_Files(self.path, self.proj.name + '-*', **kwargs) 
    self.fnames_actual = fnames
    return fnames
  
  # -----------------------------------------------------------------------------
  
  def Ls(self, **kwargs):
    o, e = Bash2('pwd', path=self.path)
    print('Content of ', o)

    o, e = Bash2('ls -lth ' + self.path)
    print(o)
    
  # -----------------------------------------------------------------------------

  def Tidy(self, **kwargs):
    o, e = Bash2('pwd', path=self.path)
    print('Tidying up (deleting *_*_*)', o)    
    
    o, e = Bash2('rm ' + self.path + '/*_*_*')
  
  # -----------------------------------------------------------------------------  
  
  def Rm(self, **kwargs):
    o, e = Bash2('pwd', path=self.path, **kwargs)
    print('Removing content of ', o)  
    
    o, e = Bash2('rm ' + self.path + '/*')
  
  # -----------------------------------------------------------------------------
  
  
# -------------------------------------------------------------------------------


class Proj_Input(Proj_Throughput):
  """
  
  """  
  
  # -----------------------------------------------------------------------------
  
  def __init__(self, proj, **kwargs):
    """

    """  
    super(Proj_Input, self).__init__(proj, **kwargs)
    self.path = proj.path + '/inp/'
  
  # ----------------------------------------------------------------------------- 
  
  def Init(self, **kwargs):
    """
    This is a necessary work-around 
    because proj._Init_Input refer to
    proj.inp so they must be called 
    AFTER the constructor, not from it.
    
    """
    self.proj._Init_Input(**kwargs)

  # ----------------------------------------------------------------------------- 
  
  def Prepare(self, **kwargs):
    """

    """
    self.proj._Prepare_Input(**kwargs)
    
  # ----------------------------------------------------------------------------- 

  def Check(self, **kwargs):
    """
    Checkrun is done by running Fullwave with 
    '-check' flag.
    
    """
    self.proj._Check_Input(**kwargs)

  # -----------------------------------------------------------------------------

  def Plot(self, **kwargs):
    """
    
    """  
    kwargs['save'] = Kwarg('save', True, kwargs)
    self.proj._Plot_Input(**kwargs)
    
  # ----------------------------------------------------------------------------- 


# -------------------------------------------------------------------------------


class Proj_Output(Proj_Throughput):
  """

  """  
  
  # -----------------------------------------------------------------------------
  
  def __init__(self, proj, **kwargs):
    """
    
    Parameters
    ----------
    
    Returns
    -------
    
    Notes
    -----
    NOTE: CP-Runfiles GIVE step-length VALUES
    
    """   
    super(Proj_Output, self).__init__(proj, **kwargs) 
    self.path = proj.path + '/out/'
    
  # ----------------------------------------------------------------------------- 

  def Init(self, **kwargs):
    """
    This is a necessary work-around 
    because proj._Init_Input refer to
    proj.inp so they must be called 
    AFTER the constructor, not from it.
    """    
    self.proj._Init_Output(**kwargs)

  # -----------------------------------------------------------------------------     

  def Prepare(self, **kwargs):
    """

    """
    self.proj._Prepare_Output(**kwargs)
    
  # ----------------------------------------------------------------------------- 

  def Check(self, **kwargs):
    """
    Checkrun is done by running Fullwave with '-check' flag.
    
    """
    self.proj._Check_Output(**kwargs)

  # -----------------------------------------------------------------------------

  def Plot(self, **kwargs):
    """
    
    """  
    kwargs['save'] = Kwarg('save', True, kwargs)
    self.proj._Plot_Output(**kwargs)
    
  # -----------------------------------------------------------------------------
  # GENERAL
  # -----------------------------------------------------------------------------
  
  def Files(self, **kwargs): # NOTE ?
    """
    """
    from lib_fwi_project import Project_Filenames
    proj = self.proj
    fnames_dict = Project_Filenames(proj.name, proj.io, proj.problem, in_out='out', **kwargs)
    fnames = []
    for key in fnames_dict:
      fnames.append(fnames_dict[key])
    
    self.files = fnames
    return fnames

  # -----------------------------------------------------------------------------       
  
  def Remove(self, **kwargs):
    """
    
    """
    fnames = self.Files()
    for fname in fnames:
      core = Strip(fname)
      print('Removing files: ')
      o, e = Bash2('ls ' + core + '*')
      print(o)
      o, e = Bash2('rm ' + core + '*')

  # ----------------------------------------------------------------------------- 
  # -----------------------------------------------------------------------------
  
  def Check_Synthetics(self, **kwargs): # FIXME: DEL/MOVE
    pattern = self.proj.name + '-Synthetic*vtr'
    fnames = Get_Files(self.proj.path, pattern, **kwargs)
    for i, fname in enumerate(fnames):
      n1, n2, n3, d = Read_vtr(fname)
      maxx = np.max(d)
      tiny = 1e-3
      if maxx < tiny:
        eprint('Max amplitude in file ' + fname + ' is only: ' + str(maxx) + '\n')

  # -----------------------------------------------------------------------------       
  # SUBPROJECT
  # -----------------------------------------------------------------------------   
  
  def Init_Subproject(self, it, **kwargs): 
    """
    Useful when you reload cells
    and don't want to create it again.
    
    """
    
    proj_name = self.proj.name + '_synCP' + str(it) 
    print(('Initializing a new project at iter. ' +
          str(it) + ': ' + proj_name + ' - run ' + 
          'self.Create_Subproject(it) to endow it with files'))
    print('(project "hook": self.out.synCP)')
    
    # CREATE A NEW PROJECT
    dir_name = 'synCP' + str(it).rjust(5, '0') 
    self.synCP = Proj_Syn(proj_name, self.path + dir_name, **kwargs)
    print(self.synCP.name)
    old_proj = self.proj  
    
  # -----------------------------------------------------------------------------       
  
  def Create_Subproject(self, it, **kwargs): 
    """
    Creates a subproject for synthetic
    calculation using one of the output 
    models (it=2 => CP00002 etc.).
    
    """
    #NOTE model_fname IS ASSUMED TO CONTAIN PROJECT NAME
    #core = model_fname[len(self.proj.name)+1 :-len(Ext(model_fname))-1]
    
    proj_name = self.proj.name + '_synCP' + str(it) 
    print('Creating a new project', proj_name)
    print('(for iteraton ' + str(it) + ').')
    
    # CREATE A NEW PROJECT
    self.synCP = Proj_Syn(proj_name)
    print(self.synCP.name)
    old_proj = self.proj
    #print old_proj.name
    
    # DUPLICATE ALL THE INPUT
    self.synCP.Duplicate(old_proj, dupl_cmd='cp', **kwargs)
    
    # OVERWRITE ALL THE MODELS WITH CP ONES
    file2dupl = self.proj.vp.it[it].name
    print('file2dupl', file2dupl) 
    self.synCP.vp.true.Duplicate(file2dupl, dupl_cmd='cp', **kwargs)
    #NOTE: ADD OTHER MODEL:H DELTA, EPSILON ETC.
    
  # ----------------------------------------------------------------------------- 


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# PROJECT FILES
# -------------------------------------------------------------------------------  
# -------------------------------------------------------------------------------


class Proj_File(object):
  """
  Project file.
  
  """
  
  # -----------------------------------------------------------------------------
  
  def __init__(self, file_id, proj, path, **kwargs):
    """
    
    Notes
    -----
    proj_name AND io WILL BE INHERITED FROM Proj
    
    """    
    from lib_fwi_project import Project_Suffices
    
    verbos = Kwarg('verbos', 1, kwargs)
    
    self.ID = file_id # USED BY Prepare(...)
    self.pname = proj.name
    self.path = path
    
    self.io = proj.io
    self.proj = proj
    
    # FILE NAME
    suffices = Project_Suffices(self.io)
    suffix = suffices[self.ID]
    self.ext = Ext(suffix)
    self.core = Strip(suffix)[1: ]
    self.prefix = self.pname + '-' + self.core
    self.suffix = self.core + '.' + self.ext
    self.name = self.pname + suffix 
    self.fname = self.path + self.name # FULL NAME
    
    if verbos > 5:
      print('Initialized file: ' + self.name)
      print('Prefix: ' + self.prefix)
      print('Core: ' + self.core)      
      print('Suffix: ' + self.suffix)
      print('Ext: ' + self.ext)
  
  # -----------------------------------------------------------------------------

  def Create(self, **kwargs): #EMPTY
    raise NotImplementedError('This method needs to be overwitten in a subclass')
  
  # -----------------------------------------------------------------------------

  def Duplicate(self, file2dupl, **kwargs):
    """
    
    Parameters
    ----------
    
    Returns
    -------
    
    
    Notes
    -----
    
    """  
    Duplicate(file2dupl, self.path+self.name, **kwargs)
     
  # -----------------------------------------------------------------------------
  
  def Set_Geometry(self, x_origin, y_origin, z_origin, **kwargs):
    # REDUNDANT?
    self.x_origin = x_origin
    self.y_origin = y_origin
    self.z_origin = z_origin
  
  # -----------------------------------------------------------------------------  
  
  def Read(self, **kwargs):
    """
    
    Parameters
    ----------
    
    Returns
    -------
    
    
    Notes
    -----
    
    """  
    fname = self.fname
    
    if not Exists(fname):
      raise IOError(fname + ' does not exist.')
    
    ext = Ext(fname)
    if  ext == 'vtr':
      from lib_fwi_generic import Read_vtr
      n1, n2, n3, Z = Read_vtr(fname, **kwargs)
    
    else:
      raise NotImplementedError('Extension: ' + ext)
    
    return Z  

  # -----------------------------------------------------------------------------  
  
  def Check(self, **kwargs): #EMPTY
    raise NotImplementedError('This method needs to be overwitten in a subclass')     

  # -----------------------------------------------------------------------------  

  def Modify(self, **kwargs): #EMPTY 
    raise NotImplementedError('This method needs to be overwitten in a subclass')  

  # -----------------------------------------------------------------------------

  def Plot(self, **kwargs):
    """
    We assume the file is plottable,
    if not Plot will raise an error.
    
    """
    from lib_generic_PLOTT import Plot, Save
    
    save = Kwarg('save', False, kwargs)
    
    Plot(self.fname, **kwargs)   
    
    if save:
      print('save', self.fname)
      fname = Save(self.fname, **kwargs)
    
    return 0

  # -----------------------------------------------------------------------------  
  
  def Cat(self, **kwargs):
    """
    NOTE: Should be overwritten by NotIMplementedError 
    in binaries etc.
    
    """
    fname = self.fname
    o, e = Bash2('cat ' + fname)
    print('Content of ', fname, ': ')
    print(o, e)    
  
  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------  
# MULTI-FILE CONTAINERS
# -------------------------------------------------------------------------------


class Model_Files(object):
  """
  Created only to handle it[1] etc.
  
  """

  # -----------------------------------------------------------------------------  
  
  def __init__(self, proj, ID, **kwargs):
    self.proj = proj
    self.ID = ID # IMPORTANT NEW ATTRIBUTE - ALLOWS TO DISTINGUISH GRAD ETC.

  # -----------------------------------------------------------------------------
  
  def Convert(self, **kwargs):
    for m in self.it:
      m.Convert(**kwargs)


# -------------------------------------------------------------------------------


class Wavefield_Files(Proj_File):
  """
  
  """

  # -----------------------------------------------------------------------------    
  
  def __init__(self, interfix, snap_step, proj, path, **kwargs):
    """
    
    
    Parameters
    ----------
    interfix : str 
      fw - for forward waviefield
      bw - for backpropagated wavefield
    
    Returns
    -------
    0
    
    Notes
    -----
    proj_name AND io WILL BE INHERITED FROM Proj
    
    """    
    self.pname = proj.name
    self.path = path
    self.io = proj.io
    self.proj = proj
    
    self.interfix = interfix
    self.snap_step = snap_step
    
  # -----------------------------------------------------------------------------  

  def Files(self, **kwargs):
    """
    FIXME; HANDLE DIFFERENT ITERATIONS AND taskIDs
    
    """
    snap_step = self.snap_step # Kwarg('snap_step', 0, kwargs)
    nsteps = self.proj.ns
    interfix = self.interfix
    
    if snap_step < 0:
      snap_step = -snap_step
      snaps = []
      for i in range(nsteps):
        snap = (i + 1) * snap_step
        if snap <= nsteps:
          snaps.append(snap)
        else:
          break
    
    elif snap_step == 0:
      eprint('snap_step set to: ' + str(snap_step) + ' - not preparing.\n')
      snaps = []
      
    else:
      raise NotImplementedError('Positive snap_step ' + str(snap_step))
      
    self.proj.inp.sr.Read(dx=self.proj.dx)
    sources = self.proj.inp.sr.s
    
    if interfix == 'fw':
      self.snap = {} 
      for key in sources:
        self.snap[key] = []
        for i_snap, snap in enumerate(snaps):
          #FIXME; HANDLE DIFFERENT ITERATIONS AND task IDs
          fname = self.pname + '-' + interfix + '-' + str(snap).rjust(6, '0') + '-csref' + str(int(key)).rjust(5, '0') + '-iter*-taskid*.vtr'
          self.snap[key].append(Snapshot_File(fname, key, snap, i_snap, 
                                              self.proj, self.proj.out.path, **kwargs))
    
    elif interfix == 'bw':
      pass
    
    else:
      raise NotImplementedError('Wavefield of type: ' + interfix + ' not supported.')
        
  # -----------------------------------------------------------------------------  

  #def Plot(self, **kwargs):
    #"""
    #NOTE: IN FUTURE REPLACE Plot WITH Plot_Various_Slices 
    #(NEEDS SOME COSMETICS WITH fw)
    
    #"""
    #try:
      #self.snap
    #except AttributeError:
      #self.Files(**kwargs)
    
    #for sid in self.snap:
      #for snapshot in self.snap[sid]:
        #snapshot.Plot(**kwargs) 
  
  # -----------------------------------------------------------------------------         

  def Plot(self, **kwargs):
    """
    NOTE: IN FUTURE REPLACE Plot WITH Plot_Various_Slices 
    (NEEDS SOME COSMETICS WITH fw)
    
    """
    plot_mod = Kwarg('plot_mod', False, kwargs) # HERE FOR SUFFIX
    save = Kwarg('save', False, kwargs)
    animate = Kwarg('animate', True, kwargs)
    try:
      self.snap
    except AttributeError:
      self.Files(**kwargs)
    
    for sid in self.snap:
      fnames_png = ''
      for snapshot in self.snap[sid]:
        
        # PREPARE A FNAME TO SAVE
        fname_png = Strip(snapshot.fname) 
        if plot_mod:
          fname_png += '_m' # IF WE WANT TO PLOT WITH/OUT MODELS
        #if plot_...:
          #fname_png += ...          
        fname_png += '.png ' #NOTE SPACE
        
        
        snapshot.Plot(save_as=fname_png, **kwargs)
        
        print(fname_png)
        fnames_png += fname_png
        
        #break
      
      if animate and save:
        delay = 1000. / len(self.snap[sid]) # FIXME: WORKS ONLY FOR SMALL len
        anim_name = Strip(fname_png)
        
        print('fnames', fname_png)
        anim_name += '.gif'
        o, e = Bash2('convert -delay ' + str(delay) + ' -loop 0 ' + 
                 fnames_png + ' ' + anim_name)
      #break
    #print 'brrrreakkk'

  # ----------------------------------------------------------------------------- 


# -------------------------------------------------------------------------------


class Data_File(Proj_File):
  """
  Files such as:
  - MGL1521_W188_4.sgy
  - Observed.sgy
  - RawSeis.sgy
  - OutSeis.sgy
  - Synthetic.sgy
  and their split, filtered,
  or altered in another way versions.
  
  Notes
  -----
  Formats other than .sgy not supported
  at the moment.
  
  Data-processing methods can be applied 
  to the whole 'lump' (file before splitting)
  or to gather files which are instances
  of a child class.
  
  """
  
  
  def __init__(self, file_id, proj, path, **kwargs):
    super(Data_File, self).__init__(file_id, proj, path, **kwargs)
    
    self.against = Kwarg('against', 'outseis', kwargs)
  
  # -----------------------------------------------------------------------------   
  
  def Files(self, **kwargs): #NOTE: PSEUDO-INIT. (SEE Notes)
    """
    
    Notes
    -----
    NOTE: it can't be in __init__
    because it's initalized in 
    Proj_Input.__init__ => can't 
    refer to self.proj.inp... which
    is not yet created completely.
    """
    kwargs['convert'] = Kwarg('convert', False, kwargs)
    
    # GET PROJECT'S SHOTS
    self.proj.inp.sr.Read(dx=1) #NOTE: dx WRONG <= DUMMY HERE
    sources = self.proj.inp.sr.s
    
    # DICT OF DICTS (SEE INNER LOOP)
    self.split = {}
    for sid in sources:
      fnames = Get_Files(self.path, self.prefix + '_' + sid + '_*.sgy', **kwargs)
      self.split[sid] = {}
      for fname in fnames:
        # EXTRACT THE LINE NUMBER NOTE: MAKE IT A FUNCTION?
        fname = Path_Leave(fname, **kwargs)
        core = Strip(fname, **kwargs)
        line_no = Split(core, '_')[-1]
        #NOTE: Gather_File
        self.split[sid][line_no] = Gather_File(fname, self.proj, self.path, **kwargs)  
  
  # -----------------------------------------------------------------------------       
  
  def Duplicate(self, **kwargs):
    """
    Let SegyPrep do it for you, 
      just create the RawSeis.txt
    
    It uses the files from RawSeis.txt
    => data is selected once.
    
    """
    
    #print self.proj.inp.rawseis.name
    #c = Read_File(self.proj.inp.rawseis.name, **kwargs)
    #for line in c:
      #print line
    raise NotImplementedError('Let SegyPrep do it for you, just create the RawSeis.txt!')

  # -----------------------------------------------------------------------------   

  def Split(self, **kwargs):
    """
    Into lines.
    
    Notes
    -----
    It collects file names
    using the list of project 
    shots and Get_Files().
    
    # NOTE: ONLY sgy FORMAT
    
    """
    o, e = Bash2('su_supergather_split_stations_n_lines.sh ' + 
                 self.proj.SEGY['station'] + ' ' + self.proj.SEGY['line'] + 
                 ' ' + self.name, path=self.path)
    
    self.Files(**kwargs)
    
  # -----------------------------------------------------------------------------

  def Plot(self, **kwargs):
    """
    Plot all of the split files.
    
    """
    try:
      self.split
    except AttributeError:
      self.Files(**kwargs)
    
    kwargs['plot_wiggles'] = Kwarg('plot_wiggles', True, kwargs)
    juxtapose = Kwarg('juxtapose', True, kwargs)
    


    #NOTE: Init_Input IS A MORE NATURAL PLACE
    # BUT WE DON'T WANT TO CALL IT EVERY TIME WE RELOAD LIBS AND proj
    try:
      self.split
    except AttributeError:
      self.Files(**kwargs)
    
    for sid in self.split:
      for lid in self.split[sid]:
        gather = self.split[sid][lid]
        
        # PLOT
        plt.figure()
        gather.Plot(**kwargs)  
        
        # JUXTAPOSE
        if juxtapose:
          against = Kwarg('against', self.against, kwargs)
          kwargs['against'] = against
          
          # FIXME: D-C IT WORKS IN Proj_Both
          if against == 'synth':
            try:
              self.proj.out.synth.split
            except AttributeError:
              self.proj.out.synth.Files(**kwargs)          
          elif against == 'outseis':
            try:
              self.proj.inp.outseis.split
            except AttributeError:
              self.proj.inp.outseis.Files(**kwargs)
          else:
            eprint('Cant split here for againtst: ' + 
                   against + '\n')
          
          for juxt in ['alternate_cmaps', 'wigg_on_cmap']: #, 'wigg_on_cmap', 'diff']: FIXME
            juxt_kwargs = dict(kwargs)
            if juxt == 'alternate_cmaps':
              juxt_kwargs['normalize'] = 'rms'
              juxt_kwargs['interp'] = None
              juxt_kwargs['suffix'] = '_AC'
              
            elif juxt == 'wigg_on_cmap':
              juxt_kwargs['normalize'] = None
              juxt_kwargs['interp'] = 'bilinear'
              juxt_kwargs['suffix'] = '_WOC'
          
            elif juxt == 'diff':
              juxt_kwargs['interp'] = 'bilinear'
              juxt_kwargs['suffix'] = '_DIFF'
            else:
              raise ValueError('juxt: '+ juxt)

            plt.figure()
            gather.Juxtapose(juxt=juxt, **juxt_kwargs)
        
        #break
      #break
    #print 'bbbbrk'
  
  # -----------------------------------------------------------------------------

  def Plot_2D_or_1D(self, **kwargs):
    """
    Plot data if it doesn't require
    splitting, i.e. it is a single 
    gather or even a trace.
    
    Otherwise, use self.Plot instead.
    
    """
    cmap = Kwarg('cmap', 'seismic', kwargs)
    yflip = Kwarg('yflip', 1, kwargs)
    
    kwargs['cmap'] = cmap
    kwargs['yflip'] = yflip
    kwargs['plot_type'] = 'slice'
    kwargs['slice_coord'] = 'y'
    kwargs['coord_value'] = 0
    super(Data_File, self).Plot(**kwargs)  

  # -----------------------------------------------------------------------------
  
  def Animate(self, **kwargs):
    pass
  
  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


class Raw_Data_File(Data_File):
  """
  Before processing.
  
  Notes
  -----
  This feeds on OutSeis, i.e. 
  it lets SegyPrep prepare the 
  relevant subset of data.
  
  Once processing is done it serves 
  it needs to be fed into 
  SegyPrep once again 
  (new RawSeis).
  
  SegyPrep will kill traces etc.
  
  
  """

  # -----------------------------------------------------------------------------
  
  def Create(self, **kwargs):
    """
    
    
    Notes
    -----
    Moving to RawSeis will prevent 
    new SegyPrep run (which DOES modify
    the data, e.g. kill the traces) 
    from overwriting.
    
    """
    verbos = Kwarg('verbos', 1, kwargs)
    
    if verbos > 0:
      print('Moving OutSeis.sgy to RawSeis.sgy...')
    
    o, e = Bash2('mv ' + self.proj.inp.outseis.name + ' ' + self.name)
  

# -------------------------------------------------------------------------------


class Pick_Files(Proj_File):
  """
  
  Notes
  -----
  
  
  """
  
  # -----------------------------------------------------------------------------    
  
  def __init__(self, proj, path, **kwargs):
    """
    This is yet another type of multi-file 
    containers - it contains multiple files 
    each of which (!) will be split.
    
    Parameters
    ----------
    
    Returns
    -------
    0
    
    Notes
    -----
    proj_name AND io WILL BE INHERITED FROM Proj
    
    """   
    verbos = Kwarg('verbos', 1, kwargs)
    
    self.proj = proj
    self.pname = proj.name
    self.path = path
    self.fmt = Kwarg('fmt', 'tl', kwargs)
   
   
    self.shot_file = Kwarg('shot_file', 
                           "/home/kmc3817/heavy_PhD/meta_data/shot_lines_ranges.txt", 
                           kwargs)   
    pick_path = '/home/kmc3817/heavy_PhD/picks/first_arriv/UO/'
    pick_type = Kwarg('pick_type', 'P', kwargs)
    self.pick_type = pick_type
    
    # ----------------------------------------------------------------------------- 
    # CHOOSE RELEVANT REPOSITORIES (SYNC WITH BITBUCKET kajetan_ch)
    # ----------------------------------------------------------------------------- 
    if (pick_type == 'P') or (pick_type == 'Pgpk') or (pick_type == 'PgCal'):
      pick_path += 'time-corrected-travel-time-picks/'
      
    elif pick_type == 'Pw':
      pick_path += 'time-corrected-waterwave-picks/'
    
    else:
      raise ValueError('Unknown pick_type ' + pick_type)
    
    self.pick_path = Kwarg('pick_path', pick_path, kwargs)    
    
    # ----------------------------------------------------------------------------- 
    # DEFINE INPUT FILES FORMAT (RAW FILES FROM PICKING SOFTWARE)
    # ----------------------------------------------------------------------------- 
    if self.fmt == 'tl': # TomoLab / tlPicker
      self.raw_prefix = 'tlPick_'
      self.raw_suffix = '_' + self.pick_type + '.dat'
    else:
      raise ValueError('Unknown fmt: ' + fmt)

    # DEFINE OUTPUT FORMAT (FILES PROCESSED BY FullwavePy)
    self.ID = 'Pick' + self.pick_type
    self.core = self.ID
    self.ext = 'dat'
    self.suffix = self.core + '.' + self.ext
    self.prefix = self.pname + '-' + self.core
    
    if verbos > 5:
      print('Initialized file: ' + self.name)
      print('Prefix: ' + self.prefix)
      print('Core: ' + self.core)      
      print('Suffix: ' + self.suffix)
      print('Ext: ' + self.ext)

  # -----------------------------------------------------------------------------  
  
  def Duplicate(self, **kwargs):
    """
    Copy pick-files to a project directory.
    
    Notes
    -----
    NOTE: copy only those OBSes that are 
    contained in the project box.
    
    Concatenate all tlPicker files?
    
    """
    verbos = Kwarg('verbos', 1, kwargs)
    
    path = self.pick_path
    
    if verbos > 0:
      print(('Searching for files *' + 
            self.suffix + ' in ' + path))
    self.fnames = []
    for obs in self.proj.obs_list:
      pattern = '*' + obs + self.raw_suffix
      
      fnames = Get_Files(path, pattern, **kwargs)
      if len(fnames) == 1:
        fname = fnames[0]
      elif len(fnames) == 0:
        raise IOError('In ' + path + 'found no files matching pattern: ' + pattern)
      else:
        raise IOError('In ' + path + 'found > 1 file matching pattern: ' + pattern)
      
      #NOTE: COPY WITH A NEW NAME!
      nfname = self.prefix + '_' + obs + '.' + self.ext
      o, e = Bash2('cp ' + fname + ' ' + nfname)
      if verbos > 0:
        print('Copied ' + fname + ' to ' + nfname)
      
      self.fnames.append(nfname)
    
  # -----------------------------------------------------------------------------  
  
  def Read(self, fname, **kwargs):
    """
    Read a file for each station
    (before splitting).
    
    """
    c = Read_File(fname) 
    if self.fmt == 'tl':
      sid_list = [line[2] for line in c]
      t_list = [line[4] for line in c]
    else:
      raise ValueError('Unknown fmt: ' + fmt)
    
    return sid_list, t_list

  # -----------------------------------------------------------------------------    
  
  def Files(self, **kwargs): #NOTE: PSEUDO-INIT. (SEE Notes)
    """
    
    Notes
    -----
    NOTE: it can't be in __init__
    because it's initalized in 
    Proj_Input.__init__ => can't 
    refer to self.proj.inp... which
    is not yet created completely.
    """
    
    # GET PROJECT'S SHOTS
    self.proj.inp.sr.Read(dx=1) #NOTE: dx WRONG <= DUMMY HERE
    sources = self.proj.inp.sr.s
    
    prefix = self.prefix
    
    # DICT OF DICTS (SEE INNER LOOP)
    self.split = {}
    for sid in sources:
      sid = sid[1: ] # NOTE: WE DON'T WANT TO HAVE COMPONENT NUMBER HERE
      #COMPONENT NO. ARE IN HEADERS AS 4140 ETC. AND ARE READ BY FULLWAVE
      # AND WRITTEN TO .geo FILES ETC.
      pattern = prefix + '_' + sid + '*.dat'
      fnames = Get_Files(self.path, pattern, **kwargs)
      self.split[sid] = {}
      for fname in fnames:
        # EXTRACT THE LINE NUMBER NOTE: MAKE IT A FUNCTION?
        fname = Path_Leave(fname, **kwargs)
        core = Strip(fname, **kwargs)
        line_no = core.split('_')[-1]
        self.split[sid][line_no] = Pick_File(fname, self.proj, **kwargs)   
 
  # -----------------------------------------------------------------------------    
   
  def Split(self, **kwargs):
    """
    shot_ranges : dict
      Dictionary where key is a shot line no.
      and the value is the list:
      [first_shot_id, last_shot_id].  
    **kwargs : keyword arguments, optional
    
    Returns
    -------
    picks : dict 
      Dictionary where the key is a shot line no.
      and the value is a list [X, T] where:
      - X - list of shot shot_IDs 
      - T - list of picked times corresponding to X.
    
    Notes
    -----
    Structure of the file:
    | OBS_no | component | shot ID | phase (P/S etc.) | picked time ...
    | uncertainty? | ? | ? | ? | ? | ? | author | ? | ? |
    
    So the most crucial information is in the col. 3 and 5.
    
    The file does NOT distinguish between shot lines (all of them 
    are merged together).
    
    It gets lines directly from pick files!   
    """
    verbos = Kwarg('verbos', 1, kwargs)
    
    picks_path = Kwarg('picks_path', self.proj.path, kwargs)
    shot_file = Kwarg('shot_file', self.shot_file, kwargs)
    try:
      self.fnames
    except AttributeError:
      eprint('No self.fnames defined. Duplicating pick files...\n')
      self.Duplicate(**kwargs)
    
    if verbos > 0:
      print('Using shot file', shot_file, 'to split files:', self.fnames)
    
    for fname in self.fnames:
      rid = Strip(fname).split('_')[-1]
      sid_list, t_list = self.Read(fname)
      
      st_picks = Station_Picks(self.proj, **kwargs)
      split_picks = st_picks.Split(sid_list, t_list, **kwargs)
      for line in split_picks:
        X, T = split_picks[line]
        
        # SAVE TO FILE
        nfname = Strip(fname) + '_' + line + '.' + Ext(fname)
        f = open(nfname, 'w')
        for x, t in zip(X, T):
          f.write(str(x) + ' ' + str(t) + '\n')
        f.close()
        
    self.Files(**kwargs)
  
  # -----------------------------------------------------------------------------  
  

# -------------------------------------------------------------------------------
# GRIDDED DATA
# -------------------------------------------------------------------------------


class Model_File(Proj_File):
  """
  File containing a property Not_defined 
  on a whole model grid.
  
  
  Notes
  -----
  Quite a lot of methods...
  
  """
  this_class = 'Model_File.'
  
  # -----------------------------------------------------------------------------
  # GET THE FILE
  # -----------------------------------------------------------------------------

  def Create(self, **kwargs):
    from lib_fwi_model import Model_Create
    
    ext = Ext(self.name) 
    if ext == 'sgy':
      fname_vtr = self.name[ :-len('sgy')] + 'vtr'
    else:
      raise IOError('Unsupported extension: ' + ext)
        #Preparing models for I/O other than 'fw3d' not implemented yet. Try to duplicate other model in your format or choose fw3d I/O.")
    
    
    #if Ext(self.name) != 'vtr':
      #raise IOError("Preparing models for I/O other than 'fw3d' not implemented yet. Try to duplicate other model in your format or choose fw3d I/O.")
    
    Model_Create(self.fname, dims=self.proj.dims, **kwargs)
    #Save_vtr(self.path + core, model.shape[0], model.shape[1], model.shape[2], model, **kwargs)
    o, e = Bash2('convert_vtr2sgy.sh ' + fname_vtr)

  # -----------------------------------------------------------------------------
  
  def Duplicate(self, file2dupl, **kwargs): 
    """
    WE ASSUME COPYING WITHIN ONE FORMAT
    
    """    
    self.z_origin = Kwarg('z_origin', 'none', kwargs)
    
    super(Model_File, self).Duplicate(file2dupl, **kwargs)
    
    ext = Ext(file2dupl)
    if ext != 'sgy':
      try:
        self.x_origin = kwargs['x_origin']
        self.y_origin = kwargs['y_origin']
      except KeyError:
        raise IOError('You need to provide x_origin & y_origin.')
      
      try:
        dx1, dx2, dx3 = kwargs['dxs']
      except KeyError:
        raise IOError('You need to provide dxs=[dx1, dx2, dx3] for non-SEGY files.')        
      
      self.dxs = [dx1, dx2, dx3] # DELIBERATELY OUTSIDE try
      
    
    if self.z_origin == 'none': 
      #raise IOError('You need to provide z_origin')   
      eprint('Unknown z_origin - run self.Find_Z_Origin\n') # Will try to derive it.\n') 
      #self.z_origin = self.Find_Z_Origin(**kwargs) # FIXME: ADAPT TO VTR
    
    # SAVE ALSO AS vtr FOR PLOTTING
    eprint('Duplicate: Conversion sgy2vtr switched off\n')
    #if ext == 'sgy':
    #  dims = Kwarg('dims', None, kwargs)
    #  if not dims:
    #    nx, ny, nz = self.Find_Nnodes(**kwargs)
    #  else:
    #    nx, ny, nz = dims
    #  o, e = Bash2('convert_sgy2vtr_3D.sh ' + str(int(nx)) + ' ' + self.name)

  # ---------------------------------------------------------------------------

  def Convert(self, **kwargs):
    """
    
    
    """
    this_func = self.this_class + 'Convert: '
    
    verbos = Kwarg('verbos', 1, kwargs)
    nx = Kwarg('nx', self.proj.nx, kwargs)
    force = Kwarg('force', False, kwargs)
    ext = Ext(self.name)
    
    fname_vtr = Strip(self.fname) + '.vtr'
    
    if Exists(fname_vtr):
      if verbos > 0:
        print((this_func + 'File: ' + self.name + 
              ' has already a vtr version'))
      if force:
        if verbos > 0:
          print((this_func + 'Forcing the conversion of ' + self.name))
      else:
        if verbos > 0:
          print((this_func + 'skipping conversion (force=0)'))
        return
    
    if ext == 'vtr':
      if verbos > 0:
        print((this_func + 'File: ' + self.name + ' is in vtr'))
      return
    
    elif ext == 'sgy':
      if nx:
        if verbos > 0:
          print((this_func + 'Converting file ' + self.fname))
        o, e = Bash2('convert_sgy2vtr_3D.sh ' + str(int(nx)) + ' ' + self.fname) 
      else:
        eprint(this_func + 'No nx provided. Skipping conversion.\n')
    else:
      raise ValueError('Unknown ext: ' + ext)
    
  # -----------------------------------------------------------------------------    
      
  # -----------------------------------------------------------------------------  
  # EXTRACT CERTAIN PROPERTIES
  # -----------------------------------------------------------------------------    
  
  def Find_Z_Origin(self, **kwargs):
    """
    WE ASSUME COPYING WITHIN ONE FORMAT
    # FIXME: ADAPT TO VTR
    FIXME: IT DOESN'T WORK FOR CP FILES WHERE dz IS DIFFERENT
    
    """      
    eprint('\nWarning.  IT DOESNT WORK FOR CP FILES WHERE dz IS DIFFERENT\n')
    # FIND Z-SAMPLING OF THE FILE
    dz, e = Bash2('su_range.sh ' + self.name + " | grep dt | awk '{print $2}'", **kwargs)
    print('File ' + self.name + ' has dz = ' + dz)
    self.dz = int(dz)
    
    # FIND SEA LEVEL
    sea_node = self.Find_Sea_Level(**kwargs)
    
    # DETERMINE Z-COORDINATE OF THE UPPERMOST NODES
    z_origin = -sea_node * self.dz
    return z_origin
  
  # -----------------------------------------------------------------------------
  
  def Find_Nnodes(self, **kwargs):
    """
    NOTE: Very slow for some reason.
    
    """    
    from lib_io_fullwave import Read_vtr
    
    ext = Ext(self.name)
    if ext == 'sgy':
      sx_min, e = Bash2('su_range.sh ' + self.name + " | grep sx | awk '{print $2}'", **kwargs)
      sx_max, e = Bash2('su_range.sh ' + self.name + " | grep sx | awk '{print $3}'", **kwargs)
      sy_min, e = Bash2('su_range.sh ' + self.name + " | grep sy | awk '{print $2}'", **kwargs)
      sy_max, e = Bash2('su_range.sh ' + self.name + " | grep sy | awk '{print $3}'", **kwargs)
      
      dx, e = Bash2('su_range.sh ' + self.name + " | grep dt | awk '{print $2}'", **kwargs)
      dx = float(dx)
      
      nx = int((float(sx_max) - float(sx_min)) / dx + 1)
      ny = int((float(sy_max) - float(sy_min)) / dx + 1)
      nz, e = Bash2('su_range.sh ' + self.name + " | grep ns | awk '{print $2}'", **kwargs)
      nz = int(nz)
    
    elif ext == 'vtr':
      nx, ny, nz, data = Read_vtr(self.fname, **kwargs)
    
    else:
      raise ValueError('Extension: ' + ext)
    
    print('Found nodes: ', nx, ny, nz)
    
    return nx, ny, nz
   
  # -----------------------------------------------------------------------------
  
  def Find_Sea_Level(self, **kwargs):
    """
    WE ASSUME COPYING WITHIN ONE FORMAT
    
    """        
    from lib_fwi_model import Model_Check_Sea_Level
    
    nx, ny, nz = self.Find_Nnodes(**kwargs)
    o, e = Bash2('convert_sgy2vtr_3D.sh ' + str(int(nx)) + ' ' + self.name)   
    
    n1, n2, n3, model = Read_vtr(self.fname[ :-len('sgy')]+'vtr', **kwargs)
    v_air = 0
    v_h2o = 1500
    v_h2o_uncertainty = 50  
    kwargs['v_h2o_uncertainty'] = v_h2o_uncertainty
    xy = Model_Check_Sea_Level(model, v_air, v_h2o, **kwargs)
    #print 'min, max', np.min(xy), np.max(xy)
    # FIXME
    sea_node = np.max(xy) # FIXME: ADD/SUBTRACT 0.5?
    eprint('Warning sea_node set as np.max(all found nodes) = ' + str(sea_node) + '\n')
    
    return sea_node
  
  # -----------------------------------------------------------------------------  
  
  def Find_Seabed_Prop(self, **kwargs):
    """
    Find a model value at the seabed.
    
    """
    from lib_generic_PLOTT import Plot
    
    #nx, ny, nz = self.Find_Nnodes(**kwargs)
    #o, e = Bash2('convert_sgy2vtr_3D.sh ' + str(int(nx)) + ' ' + self.name)   
    
    kwargs['alpha'] = Kwarg('alpha', 1, kwargs)
    
    nx, ny, nz, model = Read_vtr(self.fname[ :-len('sgy')]+'vtr', **kwargs)

    v_h2o = 1500
    v_h2o_uncertainty = 100     
    
    v1 = v_h2o - v_h2o_uncertainty
    v2 = v_h2o + v_h2o_uncertainty    
    
    #value_top_00 = model[0, 0, 0]
    seabed = np.zeros((nx, ny, 1)) + (-666.66)
    for x in range(nx):
      for y in range(ny):
        for z in range(nz):
          v_mod = model[x, y, z]
          if (v_mod > v2):
            seabed[x, y, 0] = v_mod
            #print 'Seabed detected with value:', v_mod
            break #NOTE: CRUCIAL
    
    Plot(seabed, data_type='FS', plot_type='map', **kwargs)
          
  # ----------------------------------------------------------------------------- 
  # QC
  # -----------------------------------------------------------------------------   
  
  def Check(self, **kwargs):
    from lib_fwi_project import Project_Input_Check_Dims
    
    fname = Kwarg('fname', self.name, kwargs)
    print(fname)
    
    ext = Ext(fname)
    if ext == 'sgy':
      # FIRST CONVERT!
      nx, ny, nz = self.Find_Nnodes(**kwargs)
      o, e = Bash2('convert_sgy2vtr_3D.sh ' + str(int(nx)) + ' ' + fname)     
    elif ext == 'vtr':
      pass
    
    else:
      raise NotImplementedError('Unsupported extension: ' + ext) 
      
    # CHECK MIN/MAX 
    n1, n2, n3, m = Read_vtr(self.pname + '-TrueVp.vtr')
    print('Header:', n1, n2, n3)
    print('Shape of the array', m.shape)
    print('min/max value of the model: ', np.min(m), np.max(m))
    
    # CHECK DIMS
    Project_Input_Check_Dims(self.pname + '-FreeSurf.vtr', self.pname + '-TrueVp.vtr', **kwargs)

  # ----------------------------------------------------------------------------- 
  # PROCESS
  # ----------------------------------------------------------------------------- 
  
  def Resize(self, **kwargs): 
    """
    CONFORMS TO THE CURRENT proj.box
    
    First update the box!
    
    Parameters
    ----------
    
    Returns
    -------
    
    
    Notes
    -----
    In sgy:
    - it needs z_origin
    - # NOTE: dx OF FILE ASSUMED SAME AS PROJECT'S
    
    In vtr:
    - it needs (x, y, z)_origin
    
    
    FIXME: ADD self.dims = ... 
    => NO NEED TO Find_Nnodes EVERY TIME
    WE WANT TO CONVERT TO vtr!
    
    
    """ 
    this_func = self.this_class + 'Resize: '
    verbos = Kwarg('verbos', 1, kwargs)
    
    if verbos > 0:
      print((this_func + 'Resizing file: ' + self.fname + 
            ' to conform to the project box...'))
    
    ext = Ext(self.name)
    
    if ext == 'sgy':
      key_x = Kwarg('key_x', 'sx', kwargs)
      key_y = Kwarg('key_y', 'sy', kwargs)
      
      x1, x2, y1, y2, z1, z2 = self.proj.box
      
      # CONVERT TO SLICING-NODES (HERE SAMPLES)
      z1n = int((z1 - self.z_origin) / self.proj.dx)
      z2n = int((z2 - self.z_origin) / self.proj.dx)
      
      args =  key_x + ' ' + str(x1) + ' ' + str(x2) + ' '
      args += key_y + ' ' + str(y1) + ' ' + str(y2) + ' ' 
      args += str(z1n) + ' ' + str(z2n) + ' ' + self.name
      
      o, e = Bash2('su_wind_xyz.sh ' + args, path=self.path, **kwargs)
      
      # CONVERT TO vtr
      #eprint('Switched off convert_sgy2vtr_3D.sh\n')
      nx = (x2 - x1) / self.proj.dx + 1 # NOTE: dx OF FILE ASSUMED SAME AS PROJECT'S
      ##nx, ny, nz = self.Find_Nnodes(**kwargs)
      o, e = Bash2('convert_sgy2vtr_3D.sh ' + str(int(nx)) + ' ' + self.name)
    
    
    elif ext == 'vtr':
      from lib_io_fullwave import Read_vtr, Save_vtr
      x1, x2, y1, y2, z1, z2 = self.proj.box
      print('Project box [m]:')
      print('x1, x2, y1, y2, z1, z2', x1, x2, y1, y2, z1, z2)
      #dx = float(self.proj.dx) # WE ASSUME dx IS THE SAME IN FILE AND IN PROJECT
      
      n1, n2, n3, model = Read_vtr(self.fname, **kwargs)
      x1f = self.x_origin 
      x2f = self.x_origin + (n1 - 1) * self.dxs[0]
      y1f = self.y_origin 
      y2f = self.y_origin + (n2 - 1) * self.dxs[1]
      z1f = self.z_origin 
      z2f = self.z_origin + (n3 - 1) * self.dxs[2]   
      print('File box [m]:')
      print('x1f, x2f, y1f, y2f, z1f, z2f', x1f, x2f, y1f, y2f, z1f, z2f)
      
      # CONVERT METRES -> NODES
      x1n = int((x1 - x1f) / self.dxs[0])
      x2n = int((x2 - x1f) / self.dxs[0])
      y1n = int((y1 - y1f) / self.dxs[1])
      y2n = int((y2 - y1f) / self.dxs[1])
      z1n = int((z1 - z1f) / self.dxs[2])
      z2n = int((z2 - z1f) / self.dxs[2])
      print('File box [nodes]:')
      print('x1n, x2n, y1n, y2n, z1n, z2n', x1n, x2n, y1n, y2n, z1n, z2n)
      #print 'x1n, x2n, y1n, y2n', x1n, x2n, y1n, y2n
      #quit()
      
      # SLICE THE ARRAY
      model = np.array(model[x1n : x2n + 1, y1n : y2n + 1, z1n : z2n + 1]) # NOTE: DOUBLE-CHECK + 1!
      print('model.shape after slicing', model.shape)
      
      Save_vtr(self.fname[ :-len('.vtr')], model.shape[0], model.shape[1], model.shape[2], model, **kwargs)
      
      # UPDATE FILE INFO TO PREVENT FURTHER RESIZING WHEN CALLED AGAIN
      self.x_origin = x1
      self.y_origin = y1    
      self.z_origin = z1
      
      o, e = Bash2('convert_vtr2sgy.sh ' + self.name)
      
    else:
      raise NotImplementedError('Unsupported extension: ' + ext)
    
  # -----------------------------------------------------------------------------
  
  def Resample_Py(self, **kwargs):
    """
    Unfortunately it is too memory-intensive.
    Need to use Fortan instead.
    ONLY FOR TINY MODELS (SLOOOW)
    
    Notes
    -----
    For Z-axis / time we could use SU's suresamp. 
    
    """
    ext = Ext(self.name)
    
    if ext == 'sgy':
      nx, ny, nz = self.Find_Nnodes(**kwargs)
      o, e = Bash2('convert_sgy2vtr_3D.sh ' + str(int(nx)) + ' ' + self.name)
      fname_vtr = self.name[ :-len('sgy')] + 'vtr'
    
    elif ext == 'vtr':
      fname_vtr = self.name
    
    else:
      raise NotImplementedError('Unsupported extension: ' + ext) 
    
    
    n1, n2, n3, model = Read_vtr(fname_vtr, **kwargs)
    n1, n2, n3 = model.shape
    
    dx_old = self.dx
    dx_new = self.proj.dx 
    
    if dx_new == dx_old:
      eprint('Model file ' + self.name + ' has already dx the same as defined for this project.\n')
      
    elif dx_new > dx_old:
      eprint('Down-sampling from ' + str(dx_old) + ' to ' + str(dx_new) + '. Assuming it is anti-alias filtered!\n')
      
    else:
      print(('Up-sampling from ' + str(dx_old) + ' to ' + str(dx_new) + '.\n'))
      
    print('Will resample from ', dx_old, 'm to', dx_new, 'm in all 3 dimensions')
    
    points = []
    values = []
    for x in range(model.shape[0]):
      for y in range(model.shape[1]):
        for z in range(model.shape[2]):
          points.append([x, y, z])
          values.append(model[x, y, z])
    
    points = np.array(points)
    values = np.array(values)
    
    dx_ratio = float(dx_old) / dx_new
    
    # NO. OF SAMPLES
    ns1 = complex(0, dx_ratio * n1)
    ns2 = complex(0, dx_ratio * n2)
    ns3 = complex(0, dx_ratio * n3)
    xgrid, ygrid, zgrid = np.mgrid[0:n1:ns1, 0:n2:ns2, 0:n3:ns3] 
    #print xgrid.shape, ygrid.shape
    
    
    #print points.shape

  # -----------------------------------------------------------------------------
 
  def Resample(self, **kwargs):
    """
    Use Mike Warner's Fortran utility ModPrep
    otherwise it is too slow for 3D models.
    
    Notes
    -----
    It requires sgy files.?
    
    """
    quantity = Kwarg('quantity', 'slowness', kwargs)
    
    dxs_new = Kwarg('dxs_new', None, kwargs)
    if not dxs_new: # MIGHT BE PROVIDED E.G. IN Put_Seabed
      dxs_new = self.proj.dxs
    
    dxs_old = self.dxs
    
    
    ext = Ext(self.name)
    
    if ext == 'sgy':
      fname_sgy = self.name
      fname_vtr = self.name[ :-len('sgy')] + 'vtr'
    
    elif ext == 'vtr':
      fname_sgy = self.name[ :-len('vtr')] + 'sgy'
      fname_vtr = self.name
      o, e = Bash2('convert_vtr2sgy.sh ' + self.name)
      
    else:
      raise NotImplementedError('Unsupported extension: ' + ext) 
    
    n1, n2, n3, model = Read_vtr(fname_vtr, **kwargs)
    dims_old = model.shape
    
    # THAT'S WHAT ModPrep DOES FOR SURE:
    n1_new = (n1 - 1) * dxs_old[0] / dxs_new[0] + 1 
    print('n1_new', n1_new)
    # BUT IT'S NOT NECESSARILY WHAT IS CONSISTENT
    
    
    self.proj.inp.mp = Modprep_File('ModPrep', self.proj, **kwargs)
    self.proj.inp.mp.Create(fname_sgy, dims_old, dxs_old, fname_sgy, dxs_new, quantity=quantity)
    self.proj.inp.mp.Run()
    
    #nx, ny, nz = self.Find_Nnodes(**kwargs)
    o, e = Bash2('convert_sgy2vtr_3D.sh ' + str(n1_new) + ' ' + fname_sgy)
    
    self.dxs = dxs_new
    
  # -----------------------------------------------------------------------------
  
  def Pad_Z(self, **kwargs):
    """
    
    """
    from math import ceil
    
    # MUST BE >> IF HIGHER-ORDER SPLINE INTERPOLATION IS USED IN SHIFT (BUT I AM TURNING IT OFF ANYWAY)
    pad_extra = Kwarg('pad_extra', 200, kwargs)
    pad_slowness = Kwarg('pad_slowness', 1. / 1500, kwargs) # s / m 
    
    # READ SLOWNESS MODEL
    mod = self.Read()
    
    # READ ELEVATION
    e = self.proj.inp.fs.Read(units='m', **kwargs)
    print('Extreme elevations as [m a.s.l.]', np.min(e) , np.max(e))
    
    top = np.max(e) 
    bot = abs(np.min(e)) 
    ## h = (np.max(e) - np.min(e)) 
    ## print 'Elevation difference between deepest trench and highest peak [m]: ', h
    #
    dz = self.dxs[-1]
    print('dz', dz)
    
    nnodes_top =  int(ceil(top / dz)) + pad_extra # ADD n MORE NODE TO BE SAFE
    nnodes_bot =  int(ceil(bot / dz)) + pad_extra # ADD n MORE NODE TO BE SAFE
    print('Number of nodes to add (top, bottom): ', nnodes_top, nnodes_bot)
    #
    print('mod.shape', mod.shape)
    new_shape = mod.shape + np.array([0, 0, nnodes_top + nnodes_bot])
    print('new_shape', new_shape)
    mod_padded = np.ones(new_shape) * pad_slowness
    print('mod_padded.shape', mod_padded.shape)
    #
    mod_padded[:, :, nnodes_top:-nnodes_bot] = mod # DOUBLE-CHECK SLICING
    
    model = np.array(mod_padded)
    Save_vtr(self.fname[ :-len('.vtr')], model.shape[0], model.shape[1], model.shape[2], model, **kwargs)
    
    #
    #Plot(s_2m_ext, plot_type='slices_xyz', layout='2+1', alpha=1)
  
  # -----------------------------------------------------------------------------  
  
  def Shear_Z(self, **kwargs):
    """
    """
    from scipy.ndimage.interpolation import shift
    
    # READ SLOWNESS MODEL
    mod = self.Read()
    dz_mod = float(self.dxs[-1])
    print('dz_mod', dz_mod)
    print('mod.shape', mod.shape)
    
    # READ ELEVATION
    e = self.proj.inp.fs.Read(units='m', **kwargs)
    print('Extreme elevations as [m a.s.l.]', np.min(e) , np.max(e))
    print('e.shape', e.shape)
    
    mod_shifted = np.array(mod)
    e_dz2m = np.array(e) / dz_mod # ELEVATION IN NODES (USING dz=2m)
    
    nx, ny, nz =  mod.shape
    for x in range(nx):
      for y in range(ny):
        # THIS COMBINATION OF PARAMS IMPOSES LINEAR INTERPOLATION
        # NOTE: MINUS
            mod_shifted[x, y] = shift(mod[x, y], -e_dz2m[x, y], mode='wrap', prefilter=False, order=1)
    
    return mod_shifted
    
  # -----------------------------------------------------------------------------
  
  def Filter_Z(self, **kwargs):
    """
    Use SU's Butterworth filter.
    
    a.k.a. smooth.
    
    """
    pass  
  
  # -----------------------------------------------------------------------------  
  # ADD SEABED
  # -----------------------------------------------------------------------------
  
  def Put_Seabed(self, **kwargs):
    """
    
    """
    dxs_old = self.dxs
    dz_new = dxs_old[2] / 20.
    dxs_new = [dxs_old[0], dxs_old[1], dz_new]
    self.Resample(dxs_new=dxs_new)
    #self.Pad_Z()
    #self.Shear_Z()
    #self.Filter_Z()
    #self.Resample(dxs_new=dxs_new)

  # -----------------------------------------------------------------------------  
  # ADD ANOMALY
  # -----------------------------------------------------------------------------
  
  def Create_Anomaly(self, **kwargs):
    from lib_fwi_model import Checker_Create, Spheres_Create
    from lib_io_fullwave import Save_vtr
    
    anom_type = Kwarg('anom_type', 'checker', kwargs)
    core = Strip(self.anom.fname) # INIT. IN Init_Input
    
    if anom_type == 'checker':
      anom = Checker_Create(self.proj, **kwargs)
    
    elif anom_type == 'spheres':
      anom = Spheres_Create(self.proj, **kwargs)
    
    else:
      raise NotImplementedError('Unsupported anom_type: ' + anom_type)
      
    # SAVE AND ADD AS AN ATTRIBUTE
    Save_vtr(core, anom.shape[0], anom.shape[1], anom.shape[2], anom, **kwargs)
  
  # -----------------------------------------------------------------------------  
  
  def Add_Anomaly(self, **kwargs):
    """
    Model before adding anomaly 
    should already be copied as StartVp
    or backuped another way.
    
    """
    verbos = Kwarg('verbos', 1, kwargs)
    
    from lib_fwi_model import Model_Add_Anomaly

    #nx = Kwarg('nx', 'none', kwargs)
    #core = self.fname[ :-len('.sgy')]
    #fname_vtr = core + '.vtr'     
    
    #prefix, suffix = Split(self.fname, 'True')
    #start_name = prefix + 'Start' + suffix
    
    
    # COPYING
    #if Exists(start_name): # PREVENT RE-ADDING ANOMALY
      #eprint('Add_Anomaly: model ' + start_name + ' exists => copying it to ' + 
             #self.fname + '\n')
      #o, e = Bash2('cp ' + start_name + ' ' + self.fname)
    
    #else: # SET THE BACKGROUND MODEL AS A START MOD
      #print('Add_Anomaly: before adding anomaly, copying ' + 
             #self.fname + ' to ' + start_name + '\n')
      #o, e = Bash2('cp ' + self.fname + ' ' + start_name)
    
    #if nx != 'none':    
      #o, e = Bash2('convert_sgy2vtr_3D.sh ' + str(int(nx)) + ' ' + self.fname)   
    #else:
      #eprint('Add_Anomaly: Warning. Using vtr file ' + fname_vtr +
             #'. Provide nx to convert the sgy version.\n')
    
    bckp = self.fname + '_bckp'
    core = Strip(self.fname)
    fname_vtr = core + '.vtr'
    
    
    if Exists(bckp):
      print(('Add_Anomaly: copying ' + bckp + ' to ' + self.fname))
      o, e = Bash2('cp ' + bckp + ' ' + self.fname, **kwargs)
    else:
      print(('Add_Anomaly: copying ' + self.fname + ' to ' + bckp))
      o, e = Bash2('cp ' + self.fname + ' ' + bckp, **kwargs)
    
    
    
    kwargs['force'] = 1
    self.Convert(**kwargs)
    n1, n2, n3, model = Read_vtr(fname_vtr, **kwargs)
    
    n1, n2, n3, anom = Read_vtr(self.anom.fname, **kwargs) #NOTE: ALREADY IN .vtr
    
    model = Model_Add_Anomaly(model, anom, **kwargs)
    Save_vtr(core, model.shape[0], model.shape[1], model.shape[2], model, **kwargs)
    
    if verbos > 0:
      eprint('Add_Anomaly: Converting to sgy with convert_vtr2sgy.sh...\n')
    o, e = Bash2('convert_vtr2sgy.sh ' + fname_vtr)
  
  # -----------------------------------------------------------------------------
  # PLOT 1D PROFILES
  # -----------------------------------------------------------------------------  
  
  def Plot_Depth_Profile(self, **kwargs):
    """
    
    """
    from lib_generic_PLOTT import Plot
    
    xy_list = Kwarg('xy_list', [], kwargs)
    mean = Kwarg('mean', False, kwargs)
    xlim = Kwarg('xlim', None, kwargs)
    ylim = Kwarg('ylim', None, kwargs)
    
    nx, ny, nz, model = Read_vtr(self.fname[ :-len('sgy')]+'vtr', **kwargs)
    
    profs = {}
    
    for xy in xy_list:
      x, y = xy
      prof = model[x, y, :]
      key = '(x,y)=(' + str(x) + ',' + str(y) + ')'
      profs[key] = prof
    
    if mean:
      prof_mean = np.zeros(nz)
      n = 0
      for x in range(nx):
        for y in range(ny):
          prof_mean += model[x, y, :]   
          n += 1
      prof_mean /= n
      profs['mean'] = prof_mean
      
    depth = list(range(nz))
    for prof in profs:
      plt.plot(profs[prof], depth, label=prof)


    #plt.grid()
    if xlim:
      plt.xlim(xlim)
    if ylim:
      plt.ylim(ylim)
    
    plt.gca().invert_yaxis()
    
    plt.gca().set_xlabel('velocity [m/s]')
    plt.gca().set_ylabel('depth [nodes]')    
    plt.legend()
  
  # -----------------------------------------------------------------------------  
  # PLOT DIFFERENT SLICES
  # ----------------------------------------------------------------------------- 
  
  def Plot_SR_Slice(self, sr_type, sr_id, **kwargs):
    suffix = Kwarg('suffix', '', kwargs)    
    kwargs['multi'] = 1
    kwargs['plot_sr'] = 1
    kwargs['sr_type'] = sr_type
    kwargs['sr_id'] = sr_id
    kwargs['plot_fw'] = 0
    kwargs['suffix'] = suffix + sr_id
    self.Plot(**kwargs)
  
  # -----------------------------------------------------------------------------  
  
  def Plot_All_SR_Slices(self, **kwargs):
    sr_type = Kwarg('sr_type', 's', kwargs)
    
    self.proj.inp.sr.Read(dx=self.proj.dx)
    
    if sr_type == 's':
      sids = self.proj.inp.sr.s
    elif sr_type == 'r':
      sids = self.proj.inp.sr.r
    else:
      raise ValueError('Wrong sr_type: ' + sr_type)
    
    for sid in sids:
      self.Plot_SR_Slice(sr_type, sid, **kwargs)

  # -----------------------------------------------------------------------------  
  
  def Plot_Middle_Slice(self, **kwargs):
    suffix = Kwarg('suffix', '', kwargs)
    multi = Kwarg('multi', False, kwargs)
     
    xyz = np.array(self.proj.dims) / 2
    #if multi:
      #xyz[0] += self.proj.eleft
      #xyz[1] += self.proj.efront
      #xyz[2] += self.proj.etop
    
    kwargs['xyz'] = [int(i) for i in xyz]
    kwargs['plot_sr'] = 0    
    kwargs['plot_fw'] = 0   
    kwargs['suffix'] = suffix + 'mid'
    self.Plot(**kwargs)

  # -----------------------------------------------------------------------------  
  
  def Plot_Edge_Slices(self, **kwargs):
    """
    
    """
    multi = Kwarg('multi', False, kwargs)
    suffix = Kwarg('suffix', '', kwargs)
    kwargs['plot_sr'] = 0    
    kwargs['plot_fw'] = 0

    #if multi: # SLIGHTLY DIFFERENT FROM Middle_Slice FOR 'DRY'
      #xyz = [self.proj.eleft,
             #self.proj.efront,
             #self.proj.etop]
    #else:
      #xyz = [0,0,0]
    
    kwargs['xyz'] = np.array([0,0,0]) #+ np.array(xyz)
    
    kwargs['suffix'] = suffix + 'edge0'
    self.Plot(**kwargs)
    
    kwargs['xyz'] = np.array(self.proj.dims) - 1 #+ np.array(xyz)
    kwargs['suffix'] = suffix + 'edge1'
    self.Plot(**kwargs)      

  # -----------------------------------------------------------------------------  

  def Plot_XYZ_Slices(self, **kwargs):
    """
    
    """
    xyz_list = self.proj.xyz_list + Kwarg('xyz_list', [], kwargs)
    suffix = Kwarg('suffix', '', kwargs)
    
    for xyz in xyz_list:
      kwargs['xyz'] = xyz
      kwargs['suffix'] = suffix + 'xyz' + str(xyz[0]) + '_' + str(xyz[1]) + '_' + str(xyz[2])
      self.Plot(**kwargs)    

  # ----------------------------------------------------------------------------- 
  
  def Plot_Various_Slices(self, **kwargs):
    plot_sr = Kwarg('plot_sr', 1, kwargs)
    plot_mid = Kwarg('plot_mid', 1, kwargs)
    plot_edge = Kwarg('plot_edge', 1, kwargs)
    plot_xyz = Kwarg('plot_xyz', 1, kwargs)
    
    if plot_sr:
      self.Plot_All_SR_Slices(**kwargs)
    
    if plot_mid:
      self.Plot_Middle_Slice(**kwargs)
    
    if plot_edge:
      self.Plot_Edge_Slices(**kwargs)
    
    if plot_xyz:
      self.Plot_XYZ_Slices(**kwargs)

  # -----------------------------------------------------------------------------
  # GENERIC PLOTTERS
  # -----------------------------------------------------------------------------
  
  def Plot(self, **kwargs):
    from lib_generic_PLOTT import Plot, Save
    
    multi = Kwarg('multi', False, kwargs)
    kwargs['alpha'] = Kwarg('alpha', 1., kwargs)
    save = Kwarg('save', False, kwargs)
    

    if multi:
      self.Plot_All(**kwargs)
    else:
      Plot(self.fname, plot_type='slices_xyz', **kwargs)
    
    if save:
      fname = Save(self.fname, **kwargs)
    
    return 0      

  # -----------------------------------------------------------------------------   
  
  def Plot_All(self, **kwargs):
    """
    
    Notes
    -----
    Saving is done in Plot() that calls 
    this method!
    
    """
    from lib_fwi_generic_PLOTT import Plot_All
    from lib_generic_PLOTT import Save
    
    plot_mod = Kwarg('plot_mod', True, kwargs)
    plot_fs = Kwarg('plot_fs', False, kwargs)
    plot_sr = Kwarg('plot_sr', False, kwargs)
    plot_fw = Kwarg('plot_fw', False, kwargs)
    plot_bw = Kwarg('plot_bw', False, kwargs)
    #save = Kwarg('save', False, kwargs)
    
    if plot_mod:
      if not 'mod' in kwargs:
        kwargs['mod'] = self.fname
    else:
      kwargs['mod'] = None
    
    if plot_fs:
      kwargs['fs'] = self.proj.fs.fname
    else:
      kwargs['fs'] = None    
    
    if plot_sr:
      self.proj.inp.sr.Read(dx=self.proj.dx)
      try:
        sr_type = kwargs['sr_type']
        sr_id = kwargs['sr_id']
      except KeyError:
        raise IOError('For plot_sr=1, you need to provide sr_type, sr_id')

      if sr_type == 's':
        kwargs['sr'] = self.proj.inp.sr.s[sr_id]
      elif sr_type == 'r':
        kwargs['sr'] = self.proj.inp.sr.r[sr_id]    
    else:
      kwargs['sr'] = None    
    
    if plot_fw:
      try:
        i_snap = kwargs['i_snap']
        sr_id = kwargs['sr_id']
      except KeyError:
        raise IOError('For plot_fw=1, you need to provide sr_id and i_snap')        
      
      try:
        fw = self.proj.out.fw.snap
      except AttributeError:
        self.proj.out.fw.Files(**kwargs)

      kwargs['fw'] = self.proj.out.fw.snap[sr_id][i_snap].fname
    else:
      kwargs['fw'] = None    
    
    # PLOT
    Plot_All(self.proj, **kwargs)
    
    return 0        
    
  # -----------------------------------------------------------------------------  
  
  def Plot_Recovered_Anom(self, **kwargs):
    """
    Subtract starting model.
    
    """
    Z_list = []
    for fname in [self.proj.vp.start.fname, self.fname]:
      fname_vtr = Strip(fname) + '.vtr'
      n1, n2, n3, Z = Read_vtr(fname_vtr, **kwargs)
      Z_list.append(Z)
    
    Z = Z_list[1] - Z_list[0]
    kwargs['mod'] = Z #NOTE
    kwargs['multi'] = 1 # OTHERWISE WOULD PLOT THE FILE
    kwargs['cmap'] = 'seismic'
    kwargs['suffix'] = 'Rec_' #NOTE
    
    self.Plot_Various_Slices(**kwargs)

  # -----------------------------------------------------------------------------

  def Animate(self, **kwargs): # GREAT
    """
    1 of 2 types of animations:
     read the file once, plot many (e.g. all) slices.
    
    """
    prefix = Strip(self.name)
    anim_name = prefix + '_anim.gif'
    
    Plot(self.fname, plot_type='slice_anim', 
         prefix=prefix, anim_name=anim_name, **kwargs)

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


class Model_File_CP(Model_File):
  """
  
  """
  
  # ----------------------------------------------------------------------------- 
  
  def __init__(self, fname, proj, path, **kwargs):
    """
    Main difference from Model_File 
    
    Parameters
    ----------
    
    Returns
    -------
    0
    
    Notes
    -----
    proj_name AND io WILL BE INHERITED FROM Proj
    
    """    
    from lib_fwi_project import Project_Suffices
    
    verbos = Kwarg('verbos', 1, kwargs)
    self.pname = proj.name
    self.path = path
    self.io = proj.io
    self.proj = proj
    
    # FILE NAME
    self.name = fname
    self.fname = self.path + self.name

  # ----------------------------------------------------------------------------- 


# -------------------------------------------------------------------------------


class Grad_File(Model_File_CP):
  """
  Gradient.
  
  """
  
  # -----------------------------------------------------------------------------   
  
  def Plot(self, **kwargs):
    kwargs['multi'] = True
    kwargs['cmap'] = 'seismic'
    kwargs['mod_type'] = 'grad' # CAN BE ANYTHING BUT 'model'
    super(Grad_File, self).Plot(**kwargs)

  # -----------------------------------------------------------------------------   
  
  def Plot_Various_Slices(self, **kwargs):
    kwargs['multi'] = True
    kwargs['cmap'] = 'seismic'
    kwargs['mod_type'] = 'grad' # CAN BE ANYTHING BUT 'model'
    super(Grad_File, self).Plot_Various_Slices(**kwargs)
  
  # ----------------------------------------------------------------------------- 


# -------------------------------------------------------------------------------


class Prec_File(Model_File_CP):
  """
  Preconditioner.
  
  """
  
  # -----------------------------------------------------------------------------   
  
  def Plot(self, **kwargs):
    kwargs['multi'] = True
    kwargs['cmap'] = 'hot'    
    kwargs['mod_type'] = 'prec' # CAN BE ANYTHING BUT 'model'
    super(Prec_File, self).Plot(**kwargs)

  # -----------------------------------------------------------------------------   
  
  def Plot_Various_Slices(self, **kwargs):
    kwargs['multi'] = True
    kwargs['cmap'] = 'hot'    
    kwargs['mod_type'] = 'prec' # CAN BE ANYTHING BUT 'model'
    super(Prec_File, self).Plot_Various_Slices(**kwargs)
  
  # ----------------------------------------------------------------------------- 
  

# -------------------------------------------------------------------------------


class Anomaly_File(Model_File):
  """
  vtr 
  
  """

  # -----------------------------------------------------------------------------  
  
  def __init__(self, fname, proj, path, **kwargs):
    """
    Main difference from Model_File 
    Parameters
    ----------
    
    Returns
    -------
    0
    
    Notes
    -----
    proj_name AND io WILL BE INHERITED FROM Proj
    
    """    
    from lib_fwi_project import Project_Suffices
    
    self.pname = proj.name
    self.path = path
    self.io = proj.io
    self.proj = proj
    
    # FILE NAME
    self.name = fname 
    self.fname = self.path + self.name

  # -----------------------------------------------------------------------------  
  
  def Plot(self, **kwargs):
    kwargs['cmap'] = 'seismic'
    super(Anomaly_File, self).Plot(**kwargs)
  
  # ----------------------------------------------------------------------------- 
  

# -------------------------------------------------------------------------------


class Snapshot_File(Model_File):
  """
  
  """
  
  # ----------------------------------------------------------------------------- 
  
  def __init__(self, fname, sid, snap, i_snap, proj, path, **kwargs):
    """
    Main difference from Model_File 
    
    Parameters
    ----------
    
    Returns
    -------
    0
    
    Notes
    -----
    proj_name AND io WILL BE INHERITED FROM Proj
    
    """    
    from lib_fwi_project import Project_Suffices
    
    verbos = Kwarg('verbos', 1, kwargs)
    #super(Snapshot_File, self).__init__('dummy_id', proj, path, **kwargs)

    self.sid = sid 
    self.snap = snap
    self.i_snap = i_snap
    
    self.pname = proj.name
    self.path = path
    
    self.io = proj.io
    self.proj = proj

    pattern = fname # THIS IS BECAUSE WE DON'T KNOW taskID A PRIORI
    fnames = Get_Files(self.path, pattern, **kwargs)
    if len(fnames) > 1:
      raise IOError('Found more than 1 file matching ' + 
                    pattern + ' in ' + self.path)
    self.name = Path_Leave(fnames[0])
    self.fname = self.path + self.name # FULL NAME    
    
  # -----------------------------------------------------------------------------

  def Plot(self, **kwargs):
    #kwargs['cmap'] = 'seismic'
    #kwargs['alpha'] = .8
    kwargs['multi'] = 1
    kwargs['plot_fw'] = 1
    kwargs['plot_sr'] = 1 # UNLESS VARIOUS SLICES
    kwargs['sr_type'] = 's'
    kwargs['sr_id'] = self.sid
    kwargs['i_snap'] = self.i_snap
    
    for mod in self.proj.mods_true:
      kwargs['mod'] = mod.fname
      super(Snapshot_File, self).Plot(**kwargs)

  # -----------------------------------------------------------------------------

  
# ------------------------------------------------------------------------------- 


class FS_File(Proj_File): #FIXME: Model_File?
  """
  
  """ 

  # -----------------------------------------------------------------------------  
  
  def Create(self, **kwargs):
    from lib_fwi_fs import FS_Create
    FS_Create(self.proj_name, **kwargs) # NOTE proj_name

  # ----------------------------------------------------------------------------- 
  
  def Duplicate(self, file2dupl, x_origin, y_origin, **kwargs):  # FIXME: MERGE WITH MODEL
    super(FS_File, self).Duplicate(file2dupl, **kwargs)  
    self.x_origin = x_origin
    self.y_origin = y_origin
    try:
      dx1, dx2, dx3 = kwargs['dxs']
    except KeyError:
      raise IOError('You need to provide dxs=[dx1, dx2, dx3] for non-SEGY files.')        
    
    self.dxs = [dx1, dx2, dx3] # DELIBERATELY OUTSIDE try
    #self.z_origin = z_origin
  
  # -----------------------------------------------------------------------------  
  
  def Convert_grd(self, file_grd, **kwargs):
    """
    Common format used e.g. in GMT.
    
    """    
    pass
  
  # -----------------------------------------------------------------------------   
  
  def Convert_xyz(self, file_xyz, dx, **kwargs):
    """
    Conversion grd->xyz can be done with GMT.
    
    Parameters
    ----------
    file_xyz : str 
      File name with path if needed.
    dx : float
      Grid spacing of xyz file [m].
    
    """
    from lib_io_fullwave import Convert_xyz2vtr
                                                
    Convert_xyz2vtr(file_xyz, dx, **kwargs)                                                
    
  # -----------------------------------------------------------------------------  
  
  def Read(self, **kwargs):
    """
    """
    units = Kwarg('units', 'nodes', kwargs)

    Z = super(FS_File, self).Read(**kwargs)
    
    if units == 'm':
      dz = self.dxs[-1]
      Z *= dz
    
    return Z
  
  # -----------------------------------------------------------------------------  

  def Check(self, **kwargs):
    from lib_fwi_project import Project_Input_Check_Dims
    Project_Input_Check_Dims(self.name, self.pname + '-TrueVp.vtr', **kwargs)
    
    n1, n2, n3, fs = Read_vtr(self.fname)
    z_min = np.min(fs)
    if z_min < 1:
      raise ValueError('At least one Z-coordinate of FS is smaller than 1 node. You need to shift it to put it inside the model grid.')
    
  # -----------------------------------------------------------------------------
  
  def Modify(self, **kwargs):  # NOTE: ALL THE BELOW AFFECT MODEL AS WELL!!!
    raise NotImplementedError('No need so far. Maybe filtering in future?')

  # -----------------------------------------------------------------------------
  
  def Resize(self, **kwargs): # FIXME: MERGE WITH MODEL
    """
    First update box!
    
    Parameters
    ----------
    
    Returns
    -------
    
    
    Notes
    -----
    
    """   
    from lib_io_fullwave import Read_vtr, Save_vtr
    
    dx = float(self.proj.dx)
    
    # DIMENSIONS [m] OF THE BOX FROM PROJECT INFO (DESIRED) 
    box_metre = self.proj.box #[ :-2] # EXCLUDE Z-DIMS
    x1, x2, y1, y2, z1, z2 = box_metre
    
    print('x1, x2, y1, y2, z1, z2', x1, x2, y1, y2, z1, z2)
    
    # CONVERT METRES -> NODES
    #box_nodes = np.array(box_metre) / dx
    
    
    
    #for node, m in zip(box_nodes, box_metre):
    #  print 'node, m', node, m
    #  if not node.is_integer():
    #    raise ValueError('Box dimension: ' + str(m) + ' not divisible by dx: ' + str(dx))
    #
    #x1n, x2n, y1n, y2n = [int(i) for i in box_nodes]
    #
    ##if (x2n - x1n) + 1 == 
    #
    
    
    # DIMENSIONS [m] OF THE BOX FROM FILE USING origin AS SET IN Duplicate
    # (TO CHANGE IN ORDER TO MATCH DESIRED PROJECT DIMENSIONS)
    n1, n2, n3, fs = Read_vtr(self.fname, **kwargs)
    x1f = self.x_origin 
    x2f = self.x_origin + (n1 - 1) * self.dxs[0]
    y1f = self.y_origin 
    y2f = self.y_origin + (n2 - 1) * self.dxs[1]
    #z1f = self.z_origin 
    #z2f = self.z_origin + (n3 - 1) * dx
    #print 'x1f, x2f, y1f, y2f, z1f, z2f', x1f, x2f, y1f, y2f, z1f, z2f
    print('File box [m]:')
    print('x1f, x2f, y1f, y2f', x1f, x2f, y1f, y2f)
    
    # CONVERT METRES -> NODESs
    x1n = int((x1 - x1f) / self.dxs[0])
    x2n = int((x2 - x1f) / self.dxs[0])
    y1n = int((y1 - y1f) / self.dxs[1])
    y2n = int((y2 - y1f) / self.dxs[1])
    #z1n = int((z1 - z1f) / dx)
    #z2n = int((z2 - z1f) / dx)
    #print 'x1n, x2n, y1n, y2n, z1n, z2n', x1n, x2n, y1n, y2n, z1n, z2n
    print('File box [nodes]:')
    print('x1n, x2n, y1n, y2n', x1n, x2n, y1n, y2n)
    
    # SLICE THE ARRAY
    fs = np.array(fs[x1n : x2n + 1, y1n : y2n + 1]) # NOTE: DOUBLE-CHECK + 1!
    print('fs.shape after slicing', fs.shape)
    
    nz = 1
    Save_vtr(self.fname[ :-len('.vtr')], fs.shape[0], fs.shape[1], nz, fs, **kwargs)
  
    # UPDATE FILE INFO TO PREVENT FURTHER RESIZING WHEN CALLED AGAIN
    self.x_origin = x1
    self.y_origin = y1
    #self.z_origin = z1  
  
  # -----------------------------------------------------------------------------  

  def Run_Fsprep(self, **kwargs):
    pass
  
  # -----------------------------------------------------------------------------
  
  def Plot(self, **kwargs):
    from lib_generic_PLOTT import Plot
     # REVERSED BECAUSE OUR Z-AXIS POINTS BELOW FS
    kwargs['cmap'] = 'cmo.topo_r'
    Plot(self.fname, data_type='FS', plot_type='slices_xyz', **kwargs)

  # -----------------------------------------------------------------------------


# ------------------------------------------------------------------------------- 
# WIGGLY DATA
# ------------------------------------------------------------------------------- 


class Wavelet_File(Proj_File):
  """
  
  """  

  # -----------------------------------------------------------------------------  

  def Create(self, **kwargs):
    """
    Synthetic wavelet, e.g. Ricker.
    
    """
    from lib_fwi_wavelet import Wavelet_Create
    Wavelet_Create(self.proj.name, dt=self.proj.dt, 
                   ns=self.proj.ns, path=self.path, **kwargs) # NOTE proj_name
    o = Bash('convert_sgy2vtr_2D.sh ' + self.fname)

  # ----------------------------------------------------------------------------- 

  def Duplicate(self, file2dupl, **kwargs):
    super(Wavelet_File, self).Duplicate(file2dupl, **kwargs)
    o, e = Bash2('convert_sgy2vtr_2D.sh ' + self.fname)    
    
    eprint('\n\n CAUTION: remember to re-run SegyPrep!\n\n')

  # -----------------------------------------------------------------------------  

  def Derive(self, **kwargs): #NOTE
    """
    From data.
    
    Can one automate this?
    
    """
    pass
    
  # ----------------------------------------------------------------------------- 

  def Filter(self, freqs, **kwargs):
    """
    SU'S Butterworth filter
    
    Later we can improve(?) filtering 
    by padding (Jack: 75% zero padding), 
    adding a bit of noise (Jack: 0.01), etc.
    but this seems to make little difference.
    
    """
    zerophase = Kwarg('zerophase', False, kwargs)
    filt_type = Kwarg('filt_type', 'bandpass', kwargs)
    
    ext = Ext(self.name, **kwargs)
    core = self.name[ :-(len(ext)+1)]
    #core2 = core + '.' + ext
    
    freqs = [str(i) for i in freqs]


    if zerophase:
      phase_str = '_zerophase'
    else:
      phase_str = '_minphase'
    zerophase = str(zerophase)
    
    
    if filt_type == 'bandpass':
      if len(freqs) != 4:
        raise IOError('You have to provide list of 4 frequencies for a bandpass filter.')
        
      fstoplo, fpasslo, fpasshi, fstophi = freqs
      nfname = core + '_' + filt_type + phase_str + fstoplo + '-' + fpasslo + '-' + fpasshi + '-' + fstophi + '.' + ext
      
      #print nfname
      o, e = Bash2('segyread tape=' + self.fname + ' | ' + 
                   'subfilt' + 
                   ' fstoplo=' + fstoplo +
                   ' fpasslo=' + fpasslo + 
                   ' fpasshi=' + fpasshi +
                   ' fstophi=' + fstophi + 
                   ' zerophase=' + zerophase + ' | ' +
                   'segyhdrs | ' + 
                   'segywrite tape=' + nfname)    
    
    elif filt_type == 'highpass':
      if len(freqs) != 2:
        raise IOError('You have to provide list of 2 frequencies for a highpass filter.')
        
      fstoplo, fpasslo = freqs
      nfname = core + '_' + filt_type + phase_str + fstoplo + '-' + fpasslo + '.' + ext
      
      #print nfname
      o, e = Bash2('segyread tape=' + self.fname + ' | ' + 
                   'subfilt' + 
                   ' locut=1 ' + 
                   ' hicut=0 ' +
                   ' fstoplo=' + fstoplo +
                   ' fpasslo=' + fpasslo + 
                   ' zerophase=' + zerophase + ' | ' +
                   'segyhdrs | ' + 
                   'segywrite tape=' + nfname)         
      
    else:
      raise NotImplementedError('Unsupported filter type: ' + filt_type)
    
    print('Wavelet_File.Filter: backuping ' + self.fname)
    o, e = Bash2('cp ' + self.fname + ' ' + self.fname + '_bckp')
    print('Wavelet_File.Filter: replacing old wavelet ' + self.fname)
    o, e = Bash2('cp ' + nfname + ' ' + self.fname)
    
    
    
  # -----------------------------------------------------------------------------   
  
  def Plot(self, **kwargs):
    from lib_fwi_generic_PLOTT import Plot_Wavelet_Full
    
    save = Kwarg('save', False, kwargs)
    kwargs['nx'] = 1
    Plot_Wavelet_Full(self.fname, self.proj.dt, **kwargs)
    if save:
      fname = Strip(self.fname) + '.png'
      plt.savefig(fname)
      plt.close()
      

  # ----------------------------------------------------------------------------- 


# -------------------------------------------------------------------------------


class Gather_File(Proj_File):
  """
  Split, filtered, or altered in 
  another way versions of files 
  such as:
  - MGL1521_W188_4.sgy
  - Observed.sgy
  - RawSeis.sgy
  - OutSeis.sgy
  - Synthetic.sgy  
  
  Notes
  -----
  HANDY!
  
  Formats other than .sgy not supported
  at the moment.
  
  """

  # -----------------------------------------------------------------------------   
  # BASICS
  # -----------------------------------------------------------------------------   
  
  def __init__(self, fname, proj, path, **kwargs):
    """
    Main difference from Model_File 
    Parameters
    ----------
    
    Returns
    -------
    0
    
    Notes
    -----
    proj_name AND io WILL BE INHERITED FROM Proj
    
    """    
    from lib_fwi_project import Project_Suffices
    
    this_func = 'Gather_File.__init__: '
    
    verbos = Kwarg('verbos', 1, kwargs)
    convert = Kwarg('convert', True, kwargs)
    
    self.pname = proj.name
    self.path = path
    self.io = proj.io
    self.proj = proj
    
    # FILE NAME
    self.name = fname
    self.fname = self.path + self.name
    
    # THIS IS A FIXED FORMAT ASSUMED FOR THESE FILES
    fname_split = Split(fname, '_')
    sid = fname_split[-2] # SHOT ID
    lid = Strip(fname_split[-1]) # LINE ID
    self.sid = sid 
    self.lid = lid
    self.core = fname_split[0][len(self.pname)+1: ]
    
    if convert:
      if verbos > 0:
        print(this_func + 'converting ' + self.name + ' with convert_sgy2vtr_2D.sh')
        print('(set convert=False to search for old vtr files instead of conversion).')
      o, e = Bash2('convert_sgy2vtr_2D.sh ' + self.fname)
    else:
      if verbos > 0:
        eprint(this_func + 'Conversion disabled. Will use .vtr versions if present.\n') 

  # -----------------------------------------------------------------------------   

  def Split(self, **kwargs):
    """
    
    """
    raise ValueError('Gather has already been split!')
    
  # -----------------------------------------------------------------------------
  
  def Gethw(self, **kwargs):
    """
    Get header-word values of the fname SEG-Y file
    and save them into a fname_keyword.txt.
    fldr - SEGY keyword for shot id (channel no.).
    
    Returns
    -------
    keys : list 
      List of keys for each trace.
    
    """
    verbos = Kwarg('verbos', 1, kwargs)
    
    keyword = Kwarg('keyword', 'fldr', kwargs)
    o, e = Bash2('su_gethw.sh ' + keyword + ' ' + self.name)
    if verbos > 3:
      print(o, e)
    
    fname = Strip(self.fname) + '_' + keyword + '.txt'
    
    c = Read_File(fname, **kwargs)
    keys = []
    for i in c:
      keys.append(int(i[0]))
    
    return keys
    
  # -----------------------------------------------------------------------------
  
  def Get_Chnls(self, **kwargs):
    """
    Get channel IDs.
    
    """
    verbos = Kwarg('verbos', 1, kwargs)
    
    chnls = self.Gethw(keyword='fldr', **kwargs)
    self.chnls = chnls
    
  # -----------------------------------------------------------------------------    
  
  def Check_Chnls(self, **kwargs):
    """
    Check continuity.
    
    """    
    try:
      self.chnls
    except AttributeError:   
      self.Get_Chnls(**kwargs)
    
    chnls = self.chnls
      
    
    err = False
    for i in range(len(chnls) - 1): # EXCLUDE THE LAST ONE
      chnl = chnls[i]
      chnl_next = chnls[i+1]
      if chnl_next - chnl != 1:
        err = True
        eprint('File ' + self.name + ' has a gap between channels ' + str(chnl) + ' and ' + str(chnl_next) + '\n')
    
    if err:
      raise IOError('File ' + self.name + ' has discontinuous channels (see warnings).')
    
    print('Check successful. File ' + self.name + ' has continuous channels.')
    print('No. of channels:', len(chnls))
    print('Range: ', chnls[0], '-', chnls[-1])

  # -----------------------------------------------------------------------------  
  # PROCESSING
  # ----------------------------------------------------------------------------- 

  def Filter(self, sid, line, component, freqs, **kwargs):
    """
    SU'S Butterworth filter
    
    Later we can improve filtering 
    by padding (Jack: 75% zero padding), 
    adding a bit of noise (Jack: 0.01), etc.
    but this seems to make little difference.
    
    """
    ext = Ext(self.name, **kwargs)
    core = self.name[ :-(len(ext)+1)]
    core2 = core + '_' + str(component) + str(sid) + '_' + str(line) 
    fname = core2 + '.' + ext
    
    freqs = [str(i) for i in freqs]

    zerophase = Kwarg('zerophase', False, kwargs)
    if zerophase:
      phase_str = '_zerophase'
    else:
      phase_str = '_minphase'
    zerophase = str(zerophase)
    
    filt_type = Kwarg('filt_type', 'bandpass', kwargs)
    
    if filt_type == 'bandpass':
      if len(freqs) != 4:
        raise IOError('You have to provide list of 4 frequencies of a bandpass.')
        
      fstoplo, fpasslo, fpasshi, fstophi = freqs
      nfname = core2 + '_' + filt_type + phase_str + fstoplo + '-' + fpasslo + '-' + fpasshi + '-' + fstophi + '.' + ext
      print(nfname)

      o, e = Bash2('segyread tape=' + fname + ' | subfilt' + 
                   ' fstoplo=' + fstoplo +
                   ' fpasslo=' + fpasslo + 
                   ' fpasshi=' + fpasshi +
                   ' fstophi=' + fstophi + 
                   ' zerophase=' + zerophase + ' | ' +
                   'segyhdrs | ' + 
                   'segywrite tape=' + nfname)    
    
    else:
      raise NotImplementedError('Unsupported filter type: ' + filt_type)

  # -----------------------------------------------------------------------------
  
  def Extract(self, **kwargs): #EMPTY
    """
    Using keyword(s), e.g.
    for some offset range.
    
    """
    pass
    
  # -----------------------------------------------------------------------------

  def Kill_Traces(self, **kwargs): #Add_To_Header #EMPTY
    """
    Actually flag them and 
    let SegyPrep do the rest.
    
    """
    pass

  # -----------------------------------------------------------------------------

  def Add_Picks(self, **kwargs): #Add_To_Header #EMPTY
    pass

  # -----------------------------------------------------------------------------
  
  def Add_To_Header(self, **kwargs): #EMPTY
    pass
  
  # -----------------------------------------------------------------------------  
  
  def Mute(self, **kwargs): #EMPTY
    pass

  # -----------------------------------------------------------------------------

  def Align(self, **kwargs): #EMPTY
    """
    On picks.
    
    """
    pass
  
  # -----------------------------------------------------------------------------

  def Resample(self, **kwargs): #EMPTY
    """
    Or let SegyPrep do it?
    
    """
    pass
  
  # -----------------------------------------------------------------------------

  def Window(self, **kwargs): #EMPTY
    """
    Or let SegyPrep or Fullwave3D do it?
    
    """
    pass

  # -----------------------------------------------------------------------------    
  # PLOT
  # ----------------------------------------------------------------------------- 
  
  def Juxtapose(self, against, **kwargs):
    from lib_generic_PLOTT import Save
    from lib_fwi_generic_PLOTT import Plot_Juxtaposition
    
    save = Kwarg('save', False, kwargs)
    fname1 = self.fname 
    kwargs['suffix'] = Kwarg('suffix', '', kwargs)
    
    if against == 'outseis':
      fname2 = self.proj.inp.outseis.split[self.sid][self.lid].fname
    elif against == 'synth':
      fname2 = self.proj.out.synth.split[self.sid][self.lid].fname
    elif against == 'CP':
      fname2 = self.proj.out.synCP.out.synth.split[self.sid][self.lid].fname
    else:
      fname2 = against
    
    Z_list = []
    for fname in [fname1, fname2]:
      fname = Strip(fname) + '.vtr'
      n1, n2, n3, Z = Read_vtr(fname, **kwargs)
      Z_list.append(Z)
    
    Plot_Juxtaposition(Z_list, **kwargs)

    if save:
      print(kwargs['suffix'], fname1, fname2)
      core1 = Path_Leave(Strip(fname1))
      core2 = Path_Leave(Strip(fname2))
      #print 'cc', core1, core2
      kwargs['save_as'] = Path(fname1) + core1 + '_vs_' + core2 + kwargs['suffix']  + '.png'
      fname = Save(self.fname, **kwargs)    
      #print 'fname', fname
      #print 'eeeee'
      #quit()
    
  # -----------------------------------------------------------------------------    

  def Plot(self, **kwargs):
    """
    This works only if the data is 2D or 1D.
    Otherwise split and plot Gather_File.
    
    """
    from lib_generic_PLOTT import Subplots, Subplot
    #from lib_generic_CONST import figsize_default
    
    xunit = Kwarg('xunit', None, kwargs)
    plot_picks = Kwarg('plot_picks', False, kwargs)
    plot_wiggles = Kwarg('plot_wiggles', False, kwargs)
    
    kwargs['cmap'] = Kwarg('cmap', 'seismic', kwargs)
    kwargs['alpha'] = Kwarg('alpha', 1, kwargs)
    kwargs['yflip'] = Kwarg('yflip', 1, kwargs)
    kwargs['plot_type'] = 'slice'
    kwargs['slice_coord'] = 'y'
    kwargs['coord_value'] = 0
    kwargs['yunit'] = 5
    
    if xunit == 'chnls':
      try:
        self.chnls
      except AttributeError:
        self.Get_Chnls(**kwargs)
      kwargs['xaxis'] = self.chnls
    
    if plot_picks:
      picks = self.proj.inp.picks.split[self.sid[1: ]][self.lid].Read()
      scatt = [[i, 0, j] for i, j in zip(picks[0], picks[1])] #NOTE: DUMMY y TO SLICE
      kwargs['scatts'] = [scatt]
    
    #kwargs['plot_type'] = 'slice' # FIXME
    super(Gather_File, self).Plot(**kwargs)  

    if plot_wiggles:
      plt.figure()
      super(Gather_File, self).Plot(data_type='data', suffix='wigg', **kwargs)

  # -----------------------------------------------------------------------------  


# -------------------------------------------------------------------------------
# RUNFILES & LAUNCHERS
# -------------------------------------------------------------------------------


class Modprep_File(Proj_File):
  """
  
  """  

  # -----------------------------------------------------------------------------  
  
  def Create(self, f_old, dims_old, dxs_old, f_new, dxs_new, quantity='slowness', **kwargs):
    """
    
    """
    
    modprep = {}
    modprep['old file'] =  f_old
    modprep['new file'] =  f_new
    n1, n2, n3 = dims_old
    modprep['old nx1'] =  n1
    modprep['old nx2'] =  n2
    modprep['old nx3'] =  n3
    n1, n2, n3 = dxs_old
    modprep['old dx1'] =  n1
    modprep['old dx2'] =  n2
    modprep['old dx3'] =  n3
    n1, n2, n3 = dxs_new
    modprep['new dx1'] =  n1
    modprep['new dx2'] =  n2
    modprep['new dx3'] =  n3
    
    modprep['type'] =  quantity
    
    Save_Dict(self.fname, modprep)  
    
  # -----------------------------------------------------------------------------  
  
  def Run(self, **kwargs):
    o, e = Bash2('run_modprep.sh ' + self.name)  
    
  # -----------------------------------------------------------------------------

  def Log(self, **kwargs):
    o, e = Bash2('cat ' + self.name[ :-len('key')] + 'log')
    print(o, e)

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------   


class Segyprep_File(Proj_File):
  """
  
  """  

  # -----------------------------------------------------------------------------  
  
  def Create(self, **kwargs):
    from lib_fwi_generic import Segyprep_Create
    segyprep = Segyprep_Create(self.proj, **kwargs)
    
    #if segyprep['outseis'] == 'yes':
      #self.proj.inp.outseis = Data_File('OutSeis', self.proj, **kwargs)
    
  # -----------------------------------------------------------------------------
  
  def Modify(self, **kwargs):
    raise NotImplementedError('No need so far.')
  
  # -----------------------------------------------------------------------------  
  
  def Run(self, **kwargs): # NOTE CHECK IF RawSeis.txt EXISTS FOR sgy GEOMETRY
    #from lib_fwi_generic import Segyprep_Run_n_Shift_Z
    #Segyprep_Run_n_Shift_Z(self.proj.name, self.proj.z_sea, self.proj.dx, **kwargs)
    
    print('Running SegyPrep...')
    print('To see the live-progress, run SegyPrep directly with fwi_run_segyprep.sh proj.name')
    
    o, e  = Bash2('fwi_run_segyprep.sh ' + self.proj.name, path=self.path, **kwargs)
    print(o)
    eprint('\n')
    eprint('\n')
    eprint('shifting receivers for FS needs FIX!!! Now switched off!!! \n')
    eprint('\n')
    eprint('\n')
    
    #quit()
    
    
    #try:  
    #  SR_Modify_Shift_Z_geo(self.proj.name, z_sea, dx, proj_path=proj_path)  
    #except:
    #  eprint(this_func + 'Could not shift .geo files. Probably due to fw3d  format.\n')    
    
  # -----------------------------------------------------------------------------

  def Log(self, **kwargs):
    o, e = Bash2('cat ' + self.fname[ :-len('key')] + 'log')
    print(o, e)

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


#class Fsprep_File(Proj_File): # JUST A LAUNCHER?


# -------------------------------------------------------------------------------


class Runfile(Proj_File):
  """
  
  """  

  ## -----------------------------------------------------------------------------  
  
  #def __init__(self, ID, proj, path, **kwargs):
    #"""
    
    #"""
    #self.names = []
    #super(Runfile, self).__init__('runfile_synth', proj, path, **kwargs)
    #self.names.append(self.fname)
    #super(Runfile, self).__init__('runfile_inv', proj, path, **kwargs)
    #self.names.append(self.fname)  

  # -----------------------------------------------------------------------------  
  
  def Create(self, **kwargs):
    """
    
    """
    from lib_fwi_generic import Runfile_Create, Runfile_Set_Iterations
    
    Runfile_Create(self.proj, **kwargs)
    if self.proj.problem != 'synthetic':
      Runfile_Set_Iterations(self.fname, self.proj.blocks, **kwargs)

  # -----------------------------------------------------------------------------

  def Modify(self, **kwargs):
    """
    """
    from lib_fwi_generic import Runfile_Update
    from lib_fwi_generic import Runfile_Set_Iterations
    
    files = Kwarg('files', 'both', kwargs)
    
    if files == 'both':
      fnames = self.names
    elif files == 'synth':
      fnames = self.names[0]
    elif files == 'inv':
      fnames = self.names[1]
      
    for fname in fnames:    
      Runfile_Update(fname, **kwargs)
      Runfile_Set_Iterations(fname, self.proj.blocks)

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


class Runfile_CP(Runfile):
  """
  NOTE: CP-Runfiles CONTAIN step-length VALUES
  
  Mike: values $1$-$3$ are small but not tiny. 
  If they are tiny it means 'fwi is trying to do the same thing
  over and over again and it fails'
  
  """
  
  pass


# -------------------------------------------------------------------------------


class PBS_File(Proj_File):
  """
  
  """  
  
  # -----------------------------------------------------------------------------
  
  def Create(self, **kwargs):
    from lib_fwi_generic import PBS_File_Create
    PBS_File_Create(self.proj, **kwargs)

  # -----------------------------------------------------------------------------

  def Modify(self, **kwargs):
    raise NotImplementedError('Not yet')

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------
# OTHER (TEXT) FILES
# -------------------------------------------------------------------------------


class SR_Files(Proj_File):
  """
  Complications due to 2 not 1 file.
  
  Notes
  -----
  It is kept as a single hook for 2 files
  to plot both at the same type.

  """  
  
  # -----------------------------------------------------------------------------

  def __init__(self, SID, RID, proj, path, **kwargs):
    self.names = []
    self.fnames = []
    super(SR_Files, self).__init__(SID, proj, path, **kwargs)
    self.names.append(self.name)
    self.fnames.append(self.fname)
    super(SR_Files, self).__init__(RID, proj, path, **kwargs)
    self.names.append(self.name)
    self.fnames.append(self.fname)
    #self.func_prefix = 'SR_'
    self.Read(dx=self.proj.dx, **kwargs)
    
  # -----------------------------------------------------------------------------   
  
  def Create(self, **kwargs): # FIXME: IT WILL BE CALLED TWICE
    from lib_fwi_generic import SR_Prepare
    SR_Prepare(self.proj_name, **kwargs)
    
  # -----------------------------------------------------------------------------   
  
  def Check(self, **kwargs):
    self.Check_FS(**kwargs)
    #self.Check...
  
  # ----------------------------------------------------------------------------- 
  
  def Check_FS(self, **kwargs):
    from lib_fwi_generic import SR_Check_FS_Position
    SR_Check_FS_Position(self.proj.name, adjust=False)

  # -----------------------------------------------------------------------------   
  
  def Check_Model(self, **kwargs):
    pass
    #from   
  
  # ----------------------------------------------------------------------------- 
  
  def Check_Sparsity(self, **kwargs):
    pass
  
  # -----------------------------------------------------------------------------   

  def Read(self, **kwargs):
    from lib_fwi_generic import SR_Read
    dims, sources, receivers = SR_Read(self.proj.name, adjust=False, 
                                       proj_path=self.path, **kwargs)
    self.s = sources
    self.r = receivers

  # -----------------------------------------------------------------------------   
  
  def Cat(self, **kwargs):
    for fname in self.fnames:
      o, e = Bash2('cat ' + fname)
      print('Content of ' + fname + ':')
      print(o)
  
  # -----------------------------------------------------------------------------  
  
  def Show_Offsets(self, **kwargs):
    """
    Plot a histogram.
    
    """
    from lib_fwi_generic import SR_Read
    from lib_math_generic import Dist
    
    kwargs['dx'] = self.proj.dx
    kwargs['proj_path'] = self.path
    dims, sources, receivers = SR_Read(self.proj.name, adjust=False, 
                                       convert=False, **kwargs)
    
    self.s = sources
    self.r = receivers
    
    offsets = {}
    for s in sources:
      sxyz = sources[s] # RECIPROCITY => THESE ARE OBSes
      offsets[s] = []
      for r in receivers:
        rxyz = receivers[r]
        offsets[s].append(Dist(sxyz, rxyz))
        
    for s in sources:
      #plt.plot(offsets[s])
      num_bins = 20
      plt.hist(offsets[s], num_bins, alpha=0.5, label=s) 
      #plt.hist(offsets[s], num_bins, alpha=0.5) 
    
    ax = plt.gca()
    plt.title('Offsets distribution')
    ax.set_xlabel('offset [m]')
    ax.set_ylabel('counts')
    #plt.legend()
    plt.show()
      #break
    
  # -----------------------------------------------------------------------------
  
  def Plot_Pairs(self, **kwargs):
    from lib_fwi_generic_PLOTT import Plot_SR_Pairs
    
  # -----------------------------------------------------------------------------     
  
  def Plot(self, **kwargs):
    """
    
    """
    from lib_fwi_generic_PLOTT import Plot_SR
    from lib_generic_PLOTT import Save
    
    save = Kwarg('save', False, kwargs)
    
    Ss = self.proj.inp.sr.s
    Rs = self.proj.inp.sr.r

    kwargs['dims'] = self.proj.dims
    #kwargs['pairs'] = 1
    #kwargs['max_offset'] = 50
    Plot_SR(Ss, Rs, **kwargs)

    
    if save:
      fname = self.path + self.proj.name + '-SR.geo'
      fname = Save(fname, **kwargs)
    
    return 0
  
  # -----------------------------------------------------------------------------   
  
  #def Plot_Retained
  #def Plot_Extra
    
  # ----------------------------------------------------------------------------- 


# -------------------------------------------------------------------------------  


class Pick_File(Proj_File):
  """
  
  """
 
  # -----------------------------------------------------------------------------   
 
  def __init__(self, fname, proj, path, **kwargs):
    """
    Main difference from Model_File 
    
    Parameters
    ----------
    
    Returns
    -------
    0
    
    Notes
    -----
    proj_name AND io WILL BE INHERITED FROM Proj
    
    """    
    verbos = Kwarg('verbos', 1, kwargs)

    self.pname = proj.name
    self.path = path
    self.io = proj.io
    self.proj = proj
    
    # FILE NAME
    self.name = fname
    self.fname = self.path + self.name
    
    self.sid, self.lid = self.Parse_Fname(self.name)

  # -----------------------------------------------------------------------------  

  def Parse_Fname(self, fname, **kwargs):
    """
    """
    fname = Path_Leave(fname)
    fname = Strip(fname)
    #print fname
    sid = fname.split('_')[1] # FIXME: MAKE IT GENERIC
    lid = fname.split('_')[-1]    
    return sid, lid

  # -----------------------------------------------------------------------------  

  def Read(self, **kwargs):
    """
    
    """    
    c = Read_File(self.fname, **kwargs)
    chnls, times = [int(i[0]) for i in c], [float(i[1]) for i in c]
    return chnls, times

  # -----------------------------------------------------------------------------  
  
  def Export_To_Reveal(self, **kwargs):
    """
    
    Parameters
    ----------
    picks : list
      [[x, t], ...]
    
    Notes
    -----
    Uses shot IDs.
    
    """
    path = self.proj.path + '/reveal/tables/'
    fname = path + Strip(self.name) + '.tbl'
    
    
    chnls, times = self.Read(**kwargs)
    
    f = open(fname, 'w')
    
    # WRITE HEADER
    f.write('{\n')
    f.write('"version"        : 2,\n')
    f.write('"interptype"     : "f_x",\n')
    f.write('"xkey"           : "SRC_ID",\n')
    f.write('"values"         : "TIME",\n')
    f.write('"options"        : {\n')
    f.write('        "adjust_applied_stat":false\n')
    f.write('    },\n')
    f.write('"f_of_x_picks"   :\n')
    
    f.write('  [\n')
    i = 0
    for chnl, time in zip(chnls, times):
      time = int(time * 1000) # CONVERT TO WHOLE MILISECS
      f.write('   [' + str(chnl) + '       , ' + 
              str(time) + '        ]')
      if i < len(chnls) - 1:
        f.write(',\n')
      else:
        f.write('\n')
      i += 1
    
    f.write('  ]\n')
    
    # DON'T FORGET ABOUT TRAILING CURLY BRACKET
    f.write('}\n')
    f.close()

  # -----------------------------------------------------------------------------     
  
  def Create_Ascii(self, component, **kwargs): #FIXME
    """
    Ascii files are almost ready to put
    into SEGY headers.
    
    """    
    #fname = self.name(sid, line, component, **kwargs)
    #print fname
    #self.proj.name + '-' + Path_Leave(fname)
    quit()
    
    for sgy_file in sgy_files:
      fname = Path_Leave(sgy_file)
      
      
      quit() #FIXME
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
    
  # ----------------------------------------------------------------------------- 
  
  def Check(self, **kwargs):
    """
    Check continuity.
    
    """    
    chnls, times = self.Read(**kwargs)
    
    err = False
    for i in range(len(chnls) - 1): # EXCLUDE THE LAST ONE
      chnl = chnls[i]
      chnl_next = chnls[i+1]
      if chnl_next - chnl != 1:
        err = True
        eprint('File ' + self.name + ' has a gap between channels ' + str(chnl) + ' and ' + str(chnl_next) + '\n')
    
    if err:
      raise IOError('File ' + self.name + ' has discontinuous picks (see warnings).')
    
    print('Check successful. File ' + self.name + ' has continuously picked channels.')
    print('No. of channels:', len(chnls))

  # ----------------------------------------------------------------------------- 

  def Plot(self, **kwargs): # OR Plot_Raw?
    """
    sid : str 
      4140
    line : int
      27
    Notes
    -----
    Shows continuity.
    
    """
    plot_data = Kwarg('plot_data', False, kwargs)
    
    chnls, times = self.Read(**kwargs)
    #print fname, chnls, times
    
    if plot_data:
      sid, lid = self.sid, self.lid
      
      try:
        self.proj.inp.outseis.split
      except AttributeError:
        self.proj.inp.outseis.Files(**kwargs)
      
      self.proj.inp.outseis.split[sid][lid].Plot(**kwargs)
    else:
      plt.plot(chnls, times, 'o', c='lightgreen')
      plt.gca().invert_yaxis()
      plt.gca().set_xlabel('channel ID')
      plt.gca().set_ylabel('time [s]')

  # ----------------------------------------------------------------------------- 

  def Add_To_SEGY(self, **kwargs):
    pass
  
  # ----------------------------------------------------------------------------- 


# -------------------------------------------------------------------------------


class RawSeis_File(Proj_File):
  """
  This is a .txt file!
  proj-RawSeis.txt
  
  Notes
  -----
  Better to have relative paths.
  
  """ 
  
  # -----------------------------------------------------------------------------  
  
  def Create(self, **kwargs):
    """
    Only OBSes within the model - speed gain!
    
    """
    from lib_fwi_project import Project_Suffices
    
    data_path = Kwarg('data_path', 
                      '/media/kmc3817/DATADRIVE1/heavy_PhD/DATA/Santorini_2015/seismic/OBS/segy_local_coords/',
                      kwargs)
    component = Kwarg('component', 4, kwargs) 
    comp = str(component)
    prefix = Kwarg('prefix', 'MGL1521', kwargs)
    suffix = Kwarg('suffix', comp+'.sgy', kwargs)

    #NOTE: DELETE THE OLD FILE, SINCE WE'LL USE '>>' REDIRECTION 
    o, e = Bash2('rm ' + self.fname)
    
    # FILL IN WITH OBS FILES CONTAINED IN THE MODEL GRID
    obs_list = self.proj.obs_list # IT WAS SET IN Set_Geometry
    for obs in obs_list:
      pattern = prefix + '*' + obs + '*' + suffix
      fnames = Get_Files(data_path, pattern, **kwargs)
      if len(fnames) == 1:
        fname = fnames[0]
      else:
        raise IOError('In ' + path + 'found more than 1 file matching pattern: ' + pattern)
      
      o, e = Bash2('ls ' + fname + ' >> ' + self.fname)
      
      
  # -----------------------------------------------------------------------------
  
  def Preserve_Lines(self, i_min, i_max, **kwargs):
    """
      
    """
    
    c = Read_File_Not_Split(self.fname, **kwargs)
    c = c[i_min-1 : i_max] 
    f = open(self.fname, 'w')
    for line in c:
      f.write(line)
    f.close()
    
  # -----------------------------------------------------------------------------
  

# -------------------------------------------------------------------------------


class Log_File(Proj_File):
  pass # IT STILL INHERIT EVERYTHING!
  

# -------------------------------------------------------------------------------
# OTHER CLASSES  (USUALLY ATTRIBUTES OF PROJECT FILES)
# -------------------------------------------------------------------------------


class Proj_Def(object):
  """
  Move to Runfile?
  
  """
  
  # -----------------------------------------------------------------------------
  
  def __init__(self, proj, **kwargs):
    """
    """
    from lib_fwi_project_CONST import z_sea_default
    
    self.proj = proj  
    
    # SECTION A OF THE RUNFILE
    proj.problem = Kwarg('problem', 'synthetic', kwargs)
    proj.domain = Kwarg('domain', 'time', kwargs)
    proj.dim = Kwarg('dim', '3D', kwargs)
    proj.equation = Kwarg('equation', 'acoustic', kwargs)
    proj.anisotropy = Kwarg('anisotropy', 'none', kwargs)
    proj.kernel = Kwarg('kernel', 'low', kwargs)
    
    # THESE ARE SWITCHED ON IF FILES EXIST
    proj.qp = Kwarg('qp', False, kwargs) # P-WAVE QUALITY FACTOR
    proj.qs = Kwarg('qs', False, kwargs) # S-WAVE ...    
    
    # SECTION C OF THE RUNFILE
    proj.io = Kwarg('io', 'sgy', kwargs) 
    kwargs = Del('io', kwargs)
    
    #FIXME 
    proj.z_sea = Kwarg('z_sea', z_sea_default, kwargs)  

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


class Geometry(object):
  """
  """

  # -----------------------------------------------------------------------------  
  
  def __init__(self, proj, **kwargs):  
    """
    """
    verbos = Kwarg('verbos', 1, kwargs)
    timespace = Kwarg('timespace', None, kwargs)
    plot = Kwarg('plot', True, kwargs)
    
    self.proj = proj
    
    if timespace:
      box, dx, ns, dt = timespace
      kwargs['box'] = box
      kwargs['dx'] = dx
      kwargs['ns'] = ns
      kwargs['dt'] = dt
    else:
      eprint('Timespace data not provided ' + 
             '=> setting default box, dx, ns, dt\n')
      
    self.Set(**kwargs)
    
    if plot > 0:
      self.Plot(**kwargs)     

  # ----------------------------------------------------------------------------- 
  
  def Set(self, **kwargs):
    """
    
    """  
    from lib_fwi_generic import Set_Extra_Nodes
    
    dx = Kwarg('dx', 50, kwargs)
    self.proj.dx = dx # WE'LL RUN FULLWAVE ON dx=dy=dz GRIDS SO DEFINE dx FOR CONVENIENCE
    self.proj.dxs = [dx, dx, dx] # FOR RESAMPLING ETC.
    self.proj.dt = Kwarg('dt', 0.001, kwargs)
    self.proj.ns = Kwarg('ns', 1000, kwargs)
    self.proj.ttime = self.proj.ns * self.proj.dt
    
    self.proj.box = Kwarg('box', [0, 1, 0, 1, 0, 1], kwargs)
    x1, x2, y1, y2, z1, z2 = self.proj.box
    nx1 = int((x2 - x1) / self.proj.dx) + 1 
    nx2 = int((y2 - y1) / self.proj.dx) + 1  
    nx3 = int((z2 - z1) / self.proj.dx) + 1 
    self.proj.nx = nx1
    self.proj.ny = nx2
    self.proj.nz = nx3
    dims = [nx1, nx2, nx3]
    self.proj.dims = dims
    
    try:
      fname = self.proj.path + '/inp/' + self.proj.name + '-Runfile.key'
      etop, eleft, efront, ebot, eright, eback = Set_Extra_Nodes(fname, **kwargs)
      self.proj.etop = etop 
      self.proj.eleft = eleft
      self.proj.efront = efront
      self.proj.ebot = ebot
      self.proj.eright = eright
      self.proj.eback = eback
    except IOError:
      eprint('Warning ' + fname + ' not found. Cannot set etop, eleft, etc.\n')
    
    self.proj.xyz_list = Kwarg('xyz_list', [], kwargs) # LIST OF KEY SLICES
    
    self.Find_OBSes(**kwargs)  
    self.Set_Shot_Ranges(**kwargs)
    self.Set_SEGY_Mapping(**kwargs)
    
  # ----------------------------------------------------------------------------- 

  def Find_OBSes(self, **kwargs):
    """
    Determine which OBSes are contained in 
    the model.
    
    Notes
    -----
    It makes RawSeis.txt contain as few data 
    files as possible. This is essential 
    for SegyPrep runs not to be 
    ridiculuously slow (as they are 
    for ~100 huge receiver files)
    
    It takes advantage of pre-prepared 
    lists of all sources and receivers 
    (within the box used by Ben as of 15.05.2019)
    
    """
    
    # READ IN ALL SR (BIG BOX - BEN'S INVERSIONS)
    
    from lib_fwi_generic import SR_Read_All_PROTEUS   
    
    verbos = Kwarg('verbos', 1, kwargs)
    path = '/home/kmc3817/heavy_PhD/meta_data/'

    sources, receivers = SR_Read_All_PROTEUS(**kwargs)
    
    # FIND ONLY A SUBSET (FOR CURRENT PROJECT'S BOX)
    obs_list = []
    for key in receivers:
      x = receivers[key][0]
      y = receivers[key][1]
      if (x > self.proj.box[0]) and (x < self.proj.box[1]):
        if (y > self.proj.box[2]) and (y < self.proj.box[3]):
          key = key[1: ] # EXCLUDE CHANNEL NUMBER PRESENT IN SEGY HEADERS => 
          obs_list.append(key)
        
    if len(obs_list) == 0:
      raise ValueError('No OBS receivers contained in the model!')
    
    if verbos > 0:
      print('No. of OBSes within the box: ' + str(len(obs_list)))
      print('List of OBSes within the box:', obs_list)
      
    self.proj.obs_list = obs_list
    
  # -----------------------------------------------------------------------------
  
  def Set_Shot_Ranges(self, **kwargs):
    self.proj.shot_file = Kwarg('shot_file', 
                  "/home/kmc3817/heavy_PhD/meta_data/shot_lines_ranges.txt", 
                  kwargs)    
    shot_lines = Shot_Lines(self.proj.shot_file, **kwargs)
    self.proj.shot_ranges = shot_lines.Read(**kwargs)
      
  # ----------------------------------------------------------------------------- 
  
  def Set_SEGY_Mapping(self, **kwargs):
    """
    Determine mapping of SEG-Y header
    words used in the project.
    
    """
    segy = {}
    segy['station'] = 'tracf'
    segy['line'] = 'ep'
    self.proj.SEGY = segy

  # -----------------------------------------------------------------------------
  
  def Plot(self, **kwargs):
    """
    
    Notes
    -----
    FIXME: Plot different receiver types differently.
    
    FIXME: Make geometry separate class?
    proj.geom.Plot()
    
    Choosing file is at the moment
    specific to Santorini.
    
    IT'S CLUTTER WHEN YOU ANNOTATE SHOTS TOO
     plt.annotate(key, (x + shift_labels, y + shift_labels), clip_on=True) 
    clip_on IS REQUIRED, OTHERWISE BUGGY (https://github.com/matplotlib/matplotlib/issues/7895)    
    
    NOTE: .geo FILES CONTAIN COORDINATES 
    IN NODES MULTPLIED BY GRID SPACING, NOT 'REAL' (EXPERIMENT) COORDS!!!  
    
    """
    from lib_fwi_generic import SR_Read_All_PROTEUS
    from lib_generic_PLOTT import Plot_Square, Plot_Slices_XYZ
    from lib_generic_CONST import figsize_default
    
    full_map = Kwarg('full_map', False, kwargs)
    xlim = Kwarg('xlim', [-65000, 65000], kwargs)
    ylim = Kwarg('ylim', [-20000, 35000], kwargs)
    lims_pad = Kwarg('lims_pad', 5000, kwargs)
    plot_fs = Kwarg('plot_fs', 0, kwargs) # NOTE: Plot_2D
    
    if lims_pad:
      xlim = [self.proj.box[0] - lims_pad, self.proj.box[1] + lims_pad]
      ylim = [self.proj.box[2] - lims_pad, self.proj.box[3] + lims_pad]
    
    figsize = Kwarg('figsize', figsize_default, kwargs)
    s_factor = Kwarg('s_factor', 1, kwargs) # MAGNIFY MARKERS
    shift_labels = Kwarg('shift_labels', 100, kwargs) # SHIFT ANNOTATIONS
    

    # -------------------------------------------------------------------------
    # SELECT THE FILE WITH BATHYMETRY - THE SMALLEST POSSIBLE
    # EMBRACING THE PROJECT (FOR SPEED)
    # NOTE: WE ASSUME x1<x2, y1<y2
    # -------------------------------------------------------------------------
    path = '/home/kmc3817/heavy_PhD/PROJECTS/topo/'
    x1, x2, y1, y2 = self.proj.box[ :-2]
    
    # CHRISTIANA BASIN
    if x1 >= -6e4 and x2 <= 0 and y1 >= -5e3 and y2 <=15e3:
      fname = path + 'bathy_x_-6e4_0_y_-5e3_15e3_cell_50.vtr'
      extent = [-6e4, 0, -5e3, 15e3]
    
    # CALDERA
    elif x1 >= -3e4 and x2 <= 3e4 and y1 >= -5e3 and y2 <=15e3:
      fname = path + 'bathy_x_-3e4_3e4_y_-5e3_15e3_cell_50.vtr'
      extent = [-3e4, 3e4, -5e3, 15e3]
    
    # ANYDROS BASIN, ANYDROS-AMORGOS FAULT ZONE & ANAFI BASIN
    elif x1 >= 0 and x2 <= 6e4 and y1 >= -5e3 and y2 <=15e3:
      fname = path + 'bathy_x_0_6e4_y_-5e3_15e3_cell_50.vtr'
      extent = [0, 6e4, -5e3, 15e3]    
    
    # MAX EXTENT I HAVE
    else:
      fname = path + 'santorini_merged_local_xrange_-80km_80km_yrange_-40km_40km_grid_50m.vtr'
      extent = [-8e4, 8e4, -4e4, 4e4] 
      
    if full_map: # PLOT FULL BATHYMETRY REGARDLESS OF THE PROJECT BOX 
      fname = path + 'santorini_merged_local_xrange_-80km_80km_yrange_-40km_40km_grid_50m.vtr'
      extent = [-8e4, 8e4, -4e4, 4e4]       
      
    #NOTE FLIP IT FOR Plot_Map
    # IT IS ALSO DONE IN Plot_Slices_XYZ # FIXME: MERGE IT SOMEHOW
    extent[-1], extent[-2] = extent[-2], extent[-1]
    
    # -------------------------------------------------------------------------
    # PLOT
    # -------------------------------------------------------------------------    
    
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    

    # TOPOGRAPHY/BATHYMETRY
    Plot(fname, data_type='FS', plot_type='map', 
         cmap='cmo.topo_r', alpha=.8, 
         extent=extent, units='m', **kwargs)
    
    
    # PROJECT BOX
    Plot_Square(self.proj.box[0], self.proj.box[1], 
                self.proj.box[2], self.proj.box[3], 
                c='r', lw='4')
    
    
    # SOURCES AND RECEIVERS NOTE: LOOP OVER S, R
    path = '/home/kmc3817/heavy_PhD/meta_data/'
    sources, receivers = SR_Read_All_PROTEUS(**kwargs)    
    
    sx, sy = [], []
    for key in sources:

      x = sources[key][0]
      y = sources[key][1]
      sx.append(x)
      sy.append(y)   

    plt.scatter(sx, sy, s=1*s_factor, c='lightgray')
    
    sx, sy = [], []
    for key in receivers:
      x = receivers[key][0]
      y = receivers[key][1]
      
      sx.append(x)
      sy.append(y)
      plt.annotate(key, (x + shift_labels, y + shift_labels), 
                   clip_on=True) # NOTE: clip_on IS REQUIRED
    plt.scatter(sx, sy, s=10*s_factor, c='orange')    
    
    plt.xlim(xlim)
    plt.ylim(ylim)
      
  # ----------------------------------------------------------------------------- 


# -------------------------------------------------------------------------------


class Env_Vars(object):
  """
  env = Env_Vars(self)
  print env.var['dump_RawGrad']
  env.Read()
  env.Set('RawPrec', False) # REDUNDANT
  OR proj.env.var['dump_RawPrec'] = False
  env.Write()
  
  
  
  
  
  """
  
  # ----------------------------------------------------------------------------- 
  
  def __init__(self, **kwargs):
    """
    Defines the default values.
    
    """
    var = {}
    
    var['dump_fw'] = Kwarg('dump_fw', -500, kwargs) # HOW TO DISABLE? FIXME
    
    var['sched_timestamp'] = Kwarg('sched_timestamp', True, kwargs)
    var['slave_timestamp'] = Kwarg('slave_timestamp', True, kwargs)
    
    var['dump_dat'] = Kwarg('dump_dat', True, kwargs)
    var['dump_compare'] = Kwarg('dump_compare', False, kwargs)
    
    var['dump_grad'] = Kwarg('dump_Grad', True, kwargs)
    var['dump_rawgrad'] = Kwarg('dump_RawGrad', True, kwargs)
    var['dump_prec'] = Kwarg('dump_Prec', True, kwargs)
    var['dump_rawprec'] = Kwarg('dump_RawPrec', True, kwargs)
    
    self.var = var
    
  # -------------------------------------------------------------------------------
    
  def Write(self, **kwargs): # ENCODES ACTUAL NAMES FOR BASH!
    pass

  # ------------------------------------------------------------------------------


# -------------------------------------------------------------------------------


class Paths(object):
  """
  """

  # ----------------------------------------------------------------------------- 
  
  def __init__(self, proj, **kwargs):
    """
    
    """
    from lib_fwi_project import Project_Suffices
    
    proj.suffices = Project_Suffices(proj.io)
    proj.machine = Kwarg('machine', 'kmc', kwargs)
    path = Paths_Prepare(proj.machine)
    
    proj.paths = {}
    
    proj.paths['segyprep'] = path['phd_light'] + 'fullwave3D/segyprep_v3.16/bin/segyprep_v3.16'
    proj.paths['modprep'] = path['phd_light'] + 'fullwave3D/modprep/modprep.exe'
    proj.paths['fsprep'] = path['phd_light'] + 'fsprep/fsprep'
    proj.paths['fullwave'] = path['phd_light'] + 'fullwave3D/rev690/bin/fullwave3D.exe' 

  # -----------------------------------------------------------------------------           


# -------------------------------------------------------------------------------


class Station_Picks(object):
  """
  
  """
  
  # ----------------------------------------------------------------------------- 
  
  def __init__(self, proj, **kwargs):
    """
    """
    self.proj = proj
  
  # ----------------------------------------------------------------------------- 
  
  def Split(self, sid_list, t_list, **kwargs):
    """
    
    """
    shot_ranges = Kwarg('shot_ranges', self.proj.shot_ranges, kwargs)
    
    picks = {}
    for line_no in shot_ranges:
      shot_min = int(shot_ranges[line_no][0])
      shot_max = int(shot_ranges[line_no][1])
    
      X, T = [], [] 
      for i, sid in enumerate(sid_list):
        sid = float(sid)
        t = t_list[i]
        if (sid >= shot_min) and (sid <= shot_max):
          X.append(int(sid))
          T.append(float(t))
      
      # DON'T OUTPUT A LINE IF IT WASN'T PICKED FOR THIS STATION
      if len(X) == 0:
        continue
      
      picks[line_no] = [X, T] 
    
    self.picks = picks
    
    return picks

  # -----------------------------------------------------------------------------           

           
# -------------------------------------------------------------------------------


class Shot_Lines(object):
  """
  NOTE: It assumes no gaps between the first and last shot. 

  """

  # -----------------------------------------------------------------------------    
  
  def __init__(self, shot_file, **kwargs):
    self.shot_file = shot_file
    #shot_ranges = self.Read(shot_file=shot_file, **kwargs)
    #return shot_ranges

  # -----------------------------------------------------------------------------  
    
  def Read(self, **kwargs):
    """
    Read shot ranges (id of first shot, id of the last shot)
    of all shot lines
    
    Parameters
    ----------
    file_name : string
      File name. It should include extension. 
      It can include path if needed.
    **kwargs : keyword arguments, optional
    
    Returns
    -------
    shot_ranges : dict
      Dictionary where key is a shot line no.s
      and the value is the list:
      [first_shot_id, last_shot_id].
    
    Notes
    -----
    NOTE: It assumes no gaps between the first and last shot    
    
    File structure:
    
    shot_line_no first_shot_id, last_shot_id
    shot_line_no first_shot_id, last_shot_id
    ...
    
    """
    this_func = this_lib + 'Read_Shot_Lines_Ranges: '
    verbos = Kwarg('verbos', 1, kwargs)
    
    shot_file = Kwarg('shot_file', self.shot_file, kwargs)
    
    c = Read_File(shot_file)
    
    shot_ranges = {}
    for l in c:
      shot_line_no = l[0]
      first_shot_id = l[1]
      last_shot_id = l[2] 
      shot_ranges[shot_line_no] = [int(first_shot_id), int(last_shot_id)]
    
    return shot_ranges

  # -----------------------------------------------------------------------------  


# -------------------------------------------------------------------------------
  

class Functional(object):
  """
  
  """
  
  def __init__(self, proj, **kwargs):
    """
    
    """
    self.proj = proj
    self.path = self.proj.out.path
    self.name = self.proj.name + '-Functional.png' # FIXME: DUMP ASCII INSTEAD?
    self.fname = self.path + self.proj.name + '-Functional.png'
  
  # -----------------------------------------------------------------------------    
  
  def Read(self, **kwargs):
    """
    Get value of fit for all sources and for 
    all iterations.
    
    Parameters
    ----------   
    
    Returns
    -------
    functional : dict 
      Dictionary with source-IDs as keys, 
      every key contains a list of fit (%) for 
      all iterations.
    
    Notes
    -----
    PLOT FIT AS ~ SIZE, OR COLOR (TIM LIN, PHASE PLOTS)? 
    
    """
     # TRUE IF EVERY LINE OF OUTPUT LOG STARTS WITH TIME (SLAVES_SHOWTIMESTAMP="yes")
    
    
    timestamp = Kwarg('timestamp', self.proj.env.var['slave_timestamp'],
                      kwargs) # NOT SCHEDULER? FIXME
    
    if timestamp:
      add_to_index = 1
    else: 
      add_to_index = 0
    
    functional = {}
    fname = self.proj.out.path + self.proj.out.invout.name  # InvOut.log
    c = Read_File_Not_Split(fname, **kwargs)
    for line in c:
      split = line.split(None)
      # PARSE LINES CONTAINING FIT INFO
      if len(split) > 5 and split[1 + add_to_index] == 'calcResidsInfo:':
        #print split
        
        # SLAVE NO.
        first = split[0 + add_to_index]
        slave_no = first.split('(Slave')[1]
        slave_no = slave_no.split(';')[0]
        #print 'slave_no', slave_no 
        
        # SOURCE ID
        sid = first.split('CSRef')[1]
        sid = sid.split(';')[0]
        #print 'sid', sid
        
        # FIT [%] #NOTE: D-C THAT'S THE CORRECT VALUE
        percent = split[4 + add_to_index][:-1] # 4TH WORD OF LINE WITHOUT % SIGN
        
        # INITIALIZE A LIST FOR ALL ITERATIONS IF NOT DONE IT YET
        if not sid in functional:
          functional[sid] = []
        
        # APPEND THIS ITERATIONS (ACTUALLY TWICE EVERY ITERATION, SEE BELOW)
        functional[sid].append(float(percent))
    
    #NOTE: ESSENTIAL: Out.log CONTAINS 2 PIECES OF FIT INFO FOR EACH ITERATION
    # MOST PROBABLY THE SECOND ONE IS AFTER BACKPROPAGATION?! ANYWAY, WE GONNA GET 
    # RID OF IT FOR NOW BY DECIMATING LIST OF FITS (TAKING EVERY SECOND) SO THAT 
    # WE HAVE 1 NUMBER PER SOURCE PER ITERATIONS
    for key in functional:
      functional[key] = functional[key][::2]
      #print key, functional[key]
    
    # ADD TO BLOCK INFORMATION (FIXME: WE WOULD NEED TO PARSE MORE INFO, TO 
    # DISTINGUISH BETWEEN DIFFERENT BLOCKS, NO NEED FOR THIS, MAYBE IN THE FUTURE
    #for block in self.proj.inp.runfile.blocks:
    #  for iteration in range(block.niters):
    #    print 'a'
    #quit()  
    
    self.shot = functional 

  # ----------------------------------------------------------------------------- 
  
  def Plot(self, **kwargs):
    """
    
    """
    from lib_generic_PLOTT import Save
    save = Kwarg('save', False, kwargs)
    
    self.Read(**kwargs)

    plt.figure()
    
    shot = self.shot
    for shot_id in shot:
      fit_evo = shot[shot_id]
      plt.plot(list(range(1, len(fit_evo) + 1)), fit_evo, '.-', label=shot_id)
    
    
    plt.legend(prop={'size': 6})
    plt.gca().set_xlabel('Iteration')
    plt.gca().set_ylabel('Trace fit [%]')
    
    if save:
      fname = Save(self.fname, **kwargs)
    
  # -----------------------------------------------------------------------------    

    
# -------------------------------------------------------------------------------


class Iteration(object): # DEL OR HANDY REDUNDANCY?
  """  
  """
  def __init__(self, **kwargs):  
    pass


# -------------------------------------------------------------------------------


class Iteration_Block(object):
  """
  FIXME ADD ALL OTHER PARAMS FROM BLOCKS F, G, H & M!!!!! 
  
  """  
  def __init__(self, freq=3, niters=5, **kwargs):
    self.freq = freq # LOWPASS FREQUENCY
    self.niters = niters # NO. OF ITERATIONS WITHIN THE BLOCK
    self.minoff = Kwarg('minoff', None, kwargs)
    self.maxoff = Kwarg('maxoff', None, kwargs)
    
    
# -------------------------------------------------------------------------------
# SOURCES / RECEIVERS
# -------------------------------------------------------------------------------


class SR(object): #NOTE: DEL?
  """
  Source or a receiver.
  
  """
  
  # -----------------------------------------------------------------------------
  
  def __init__(self, ID, **kwargs):
    self.ID = ID

  # -----------------------------------------------------------------------------
    
  def Get_Coords(self, **kwargs):
    pass 
  
  # -----------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


class OBS(SR): #NOTE: DEL?
  '''
  NOTE: MAYBE JUST MIMIC SEG-Y?!
  
  Example of initilization:
    id = '33'
    obs = OBS(id)
    print obs.quality
  
  '''
  def Foo(self):
    pass
  #def __init__(self, ID, sensor=err_value, shot_lines=[err_value], picked_lines=[err_value], quality='good'): 
    #self.id = ID # CAPITALS CAUSE id IS RESERVED
    #self.sensor = sensor # WHOI / SIO
    #self.sensid = sensor + ID
    #self.x = err_value
    #self.y = err_value
    #self.shot_lines = shot_lines # ALONG WHICH THE OBS LIES, SEE Shot_line CLASS FOR POSSIBLE ATRIBUTES ('reshot' etc.)
    ##self.shot_line_2nd_side = shot_line_2nd_side # PARALLEL TO shot_line BUT ON THE OTHER SIDE OF THE CALDERA
    #self.picked_lines = picked_lines
    #self.quality = quality # FIXME: TAKE INTO ACCOUNT ALL COMPONENTS


# -------------------------------------------------------------------------------


class Shot(SR):
  '''
  Read from SEG-Y headers
  
  '''
  def __init__(self, sid, **kwargs):
    self.sid = sid 
    self.reciprocity = Kwarg('reciprocity', 0, kwargs)
    self.fit
    

# -------------------------------------------------------------------------------
