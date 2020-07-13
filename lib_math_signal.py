"""
Author: Kajetan Chrapkiewicz, 2018. 
All rights reserved. Ask for permision writing to K.Chrapkiewicz17@imperial.ac.uk.

This library provides framework procedures for...
At the moment fully implemented is:

1. ...

"""


## MODULES 
import sys
import numpy as np

from lib_generic import *
from lib_generic_CONST import verbos_func

## CONSTS
this_lib = 'lib_math_signal.py/'

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
# Multidimensional signal processing
# -------------------------------------------------------------------------------


def Data_Modify(g, **kwargs):
  """
  3D.
  
  Parameters
  ----------
  g : array 
    Gather in vtr format.
  func : function
    Trace-modifier.
  
  **kwargs : keyword arguments (optional)
    Current capabilities:
    trace_math : list
      List of functions to apply subsequently 
      on each trace. The order is following:
      [func1, func2, ...] 
      first func1 will be applied and so on.
      Note that order of the elements is 
      opposite to the composite function's 
      notation:
      ...(func2(func1(trace))
  
  Returns
  -------
  0
  
  Notes
  -----
  Example of usage:
  Z = Data_Modify(Toy_Data_3D(), trace_math=[lambda x: x**2], **kwargs)
  
  trace
  
  """
  this_func = this_lib + 'Data_Modify: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_math_generic import RMS
  
  trace_math = Kwarg('trace_math', [lambda x: x], kwargs)
  normalize = Kwarg('normalize', None, kwargs)
  
  # DEFINE 2 MOST COMMON TRACE OPERATIONS FOR CONVENIENCE
  if normalize == 'max':
    trace_math = [lambda x: x / np.max(x)]
  elif normalize == 'rms':
    trace_math = [lambda x: x / RMS(x)]
  
  nx, ny, nsamp = g.shape
  for x in range(nx):
    for y in range(ny):
      for func in trace_math:
        g[x][y] = func(g[x][y])
  
  if verbos > 5:
    print(this_func, 'Min, max after modifications:', np.min(g), np.max(g))
  
  #print this_func + 'END'
  return g


# -------------------------------------------------------------------------------
# WAVELETS AND IMPULSES
# -------------------------------------------------------------------------------


def Sinc(x, **kwargs):
  """
  Calculate value of the sinc function defined as:
         sinc(x) = sin(pi * x) / (pi * x).
  
  Parameters
  ----------
  x : float
    Argument of the sinc(x); Any real number.
  
  Returns
  -------
  value : float
    Value of the sinc(x).
    
  Notes
  -----
  Definition containing pi number is favorable 
  because in this case Sinc has its roots at integer 
  numbers (finite-difference nodes), not pi-multiples.
  
  Numerical stability is addressed.
  
  """
  from lib_generic_CONST import epsi
  
  pi = 3.141592654
   
  if abs(x) < epsi:
    value = 1.0
  
  else:
    value = np.sin(pi * x) / (pi * x)
  
  return value


# -------------------------------------------------------------------------------


def Dirac_Comb_1D(l, dx, pads, ampl, **kwargs):
  """
  Create a l-length series of 1D Dirac deltas
  of a given amplitude separated by dx. 

  Parameters
  ----------
  l : int 
    Lengh
  std : float 
    As above.
    
  Returns
  -------
  value : float 
    Value of the function.
  
  Notes
  -----
  Note that spikes lie in the middle of real checkers
  which are in turn obtained through convolution.
  So if spike is at z = 0, first checker will start
  at z = -0.5 * dx and end at z = 0.5 * dx.
  Nothing wrong about it.
  
  """  
  this_func = this_lib + 'Dirac_Comb_1D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  code = Kwarg('code', 'fw3d', kwargs) # FIXME: SHOULD BE UNIFIED
  
  # PREPARE 1D GRID...
  X = np.arange(1, l+1) 
  
  # ...AND INITILIAZE VALUES ON IT
  A = np.zeros(len(X))
  
  try:
    xmin = pads[0] + 1
    xmax = l - pads[1]
  except TypeError:
    eprint(this_func + 'Error. Padding must be a list of 2 values.\n')
    quit()
  
  dx = int(dx) # FIXME
  
  if verbos > 4:
    print(this_func, 'X.shape', X.shape)
    print(this_func, 'xmin, xmax, dx', xmin, xmax, dx)
  #quit()
  
  if code == 'fw3d':
    spikes = X[xmin : xmax : dx]
    
    sign = 1
    for i, x in enumerate(X):
      if x in spikes:
        add = sign * ampl
        A[i] += add
        sign *= -1
        
  elif code == 'cps':
    first = X[0] + pads[0] - 2 #FIXME?
    step = dx
    last = X[-1] - pads[1]
    
    #print 'first, last, step', first, last, step
    
    spikes = X[first : last : step]
    #print 'spikes', spikes
    sign = 1
    #no_of_spikes = pads[1] # FIXME
    
    #j = 1
    for i in range(len(X)):
      if X[i] in spikes:
        add = sign * ampl
        A[i] += add
        sign *= -1
        #if j > no_of_spikes:
          #break
        #j += 1    
  
  else:
    eprint(this_func + 'Error! Unknown code: ' + code + '\n')
    quit()
  
  #print this_func + 'END'
  return X, A


# -------------------------------------------------------------------------------


def Impulse_Rect_1D(width, length, **kwargs):
  """
  Create a rectangle OF HEIGHT=1 CENTERED IN THE SERIES   

  Parameters
  ----------
  std : float 
    As above.
    
  Returns
  -------
  value : float 
    Value of the function.
  
  Notes
  -----
  Note that spikes lie in the middle of real checkers
  which are in turn obtained through convolution.
  So if spike is at z = 0, first checker will start
  at z = -0.5 * dx and end at z = 0.5 * dx.
  Nothing wrong about it.
  
  """  
  this_func = this_lib + 'Impulse_Rect_1D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')    
  t = np.zeros(length)
  i1, i2 = Impulse_Rect_1D_Centre_Window(int(width), int(length)) # FIXME ADD INT-CHECK
  #print 'i1, i2', i1, i2
  t[i1 : i2] = 1
  return t


# -------------------------------------------------------------------------------


def Impulse_Rect_1D_Centre_Window(width, length, **kwargs):
  '''
  Return corners of a 1D rectangle pulse
  Both width and length should be integers.
  NOTE: Additionaly width should be an even number!
  
  '''
  
  this_func = this_lib + 'Centre_Window: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  if width > length:
    width = length
    eprint(this_func + 'Width of the window > length of series. Putting width = length: ' + str(length) + '\n')
  
  mode = 'new'
  if mode == 'new':
    hw = int(width / 2)
    #print this_func, 'half-width of the window: ', hw
    
    div = length % 2
    if div == 1: # ODD LENGTH
      c = int(length / 2) + 1
      #print 'centre of the odd-length list', c
      i1 = c - hw
      i2 = c + hw
    elif div == 0: # EVEN LENGTH
      c1 = int(length / 2)
      c2 = int(length / 2) + 1
      #print 'centres of the even-length list', c1, c2
      i1 = c1 - (hw - 1)
      i2 = c2 + (hw - 1)
    else:
      eprint(this_func + 'Error! Weird division by 2: ' + str(div) + '\n')
      quit()
      
    # BECAUSE PYTHON COUNTS FROM 0
    i1 -= 1
    i2 -= 1
    ##
    
  elif mode == 'old':
    hw = width / 2   
    if length % 2 == 0:
      c = length / 2 - 1
    else:
      c = length / 2
    
    print('centre of the window', c)
    print('hw', hw)
    i1 = c - hw 
    i2 = c + hw + 1 # +1 IS BECAUSE PYTHON SLICES WEIRDLY
    
    # FIXME (MAKE SHIFT CONSISTENT)
    
    if width % 2 == 0 and length % 2 != 0:
      print('case1')
      i1 = c - hw 
      i2 = c + hw
    elif width % 2 == 0 and length % 2 == 0:
      print('case2')
      i1 = c - hw + 1
      i2 = c + hw 
      if hw == 1:
        i1 = c
        i2 = c + 2
    else:
      print('another case')
      i1 = c - hw 
      i2 = c + hw + 1 # +1 IS BECAUSE PYTHON SLICES WEIRDLY
  else:
    eprint(this_func + 'Error! Unknown mode: ' + mode + '\n')
    quit()
  
  #print this_func + 'END'
  return i1, i2


# -------------------------------------------------------------------------------
# WINDOWING
# -------------------------------------------------------------------------------


def Windowed(x, f, W, **kwargs):
  """
  
  
  Parameters
  ----------
  f : function 
    To window.
  W : function 
    Windowing function.
    
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Window: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  value = f(x, **kwargs) * W(x, **kwargs)

  #print this_func + 'END'
  return value


# -------------------------------------------------------------------------------


def Window_Kaiser(x, **kwargs):
  """
  Value of the Kaiser-windowing function at point x.
  Bessel function can computed using different 
  implementations.
  
  Parameters
  ----------
  b_kaiser : float
    Shape of window
  r_kaiser : float
    Radius of window
  bessel : str 
    Implementation to choose: fw3D or python built-in.
  x : float
    Distance from the centre of sinc. 
    
  Returns
  -------
  value : float
    Value of the Kaiser-windowed sinc function at point x.
  
  Notes 
  -----
  See Hicks 2002 for details.
  
  It can't be vectorized.
  
  """
  this_func = this_lib + 'Kaiser_Window: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  b_kaiser = Kwarg('b_kaiser', 4.14, kwargs)
  r_kaiser = Kwarg('r_kaiser', 3, kwargs)
  bessel = Kwarg('bessel', 'py', kwargs)
  
  from scipy.special import i0
  from lib_generic_CONST import epsi
  
  #print this_func, 'b', b_kaiser
  
  if bessel == 'py':
    bess = lambda x : i0(x)
    
  elif bessel == 'fw3d':
    bess = lambda x : Bessel_m01_FW3D(x)
  
  else:
    raise ValueError('Unknown implementation of the Bessel function: ' + bessel + '\n')
    quit()
  
  if r_kaiser < epsi:
    eprint(this_func + 'Error. r_kaiser too small: ' + str(r_kaiser) + '\n')
    quit()
  
  frac = (x / r_kaiser) ** 2
  
  if frac > 1:
    value = 0. # NOTE: IT EXACTLY FOLLOWS THE HICKS' DEFINITION
  
  else:
    value = bess(b_kaiser * np.sqrt(1 - frac)) / bess(b_kaiser)
  
  #print this_func + 'END'
  return value


# -------------------------------------------------------------------------------
# FILTERING
# -------------------------------------------------------------------------------


def Filter_Wiener(dinp, dout, s1, **kwargs):
  from lib_io_fullwave import Save_vtr, Read_vtr
  
  for core, trace in zip(['dinp', 'dout', 's1'], [dinp, dout, s1]):
    Z = Array1D_Convert_To_3D(trace, **kwargs)
    Save_vtr(core, Z.shape[0], Z.shape[1], Z.shape[2], Z, **kwargs)
    o, e = Bash2('convert_vtr2sgy.sh ' + core + '.vtr')
    o, e = Bash2('segyread tape=' + core + '.sgy > ' + core + '.su')
  
  o, e = Bash2('sushape < s1.su wfile=dinp.su dfile=dout.su dt=1 > s2.su')
  o, e = Bash2('segyhdrs < s2.su | segywrite tape=s2.sgy')
  o, e = Bash2('convert_sgy2vtr_2D.sh s2.sgy')
  
  n1, n2, n3, s2 = Read_vtr('s2.vtr', **kwargs)
  
  s2 = s2[0][0]
  
  return s2


# ------------------------------------------------------------------------------


def Filter1D_Rect_Unnorm(data, width, **kwargs):
  # TO POLISH (EXPLAIN DOUBLING THE WIDHT BY np.convolve)
  import matplotlib.pyplot as plt

  y = Impulse_Rect_1D(width, len(data))

  #Mode same returns output of length max(M, N). Boundary effects are still visible.
  conv_array = np.convolve(data, y, mode = 'same')   
  return conv_array


# -------------------------------------------------------------------------------


def Filter1D_Gauss_Unnorm(data, sigma, **kwargs):
  # TO POLISH (EXPLAIN DOUBLING THE WIDHT BY np.convolve)
  import matplotlib.pyplot as plt
  x = np.arange(-50, 50, 0.5)
  #print 'sigma', sigma, 'fwhm', std2fwhm(sigma)
  
  # WE DIVIDE SIGMA, BECAUSE CONVOLUTION WIDENS GAUSS TWICE AS MUCH AS IT SHOULD DO!!!!!!!!!
  y = Gaussian_1D_Unnorm(x, -0.5, sigma / 2) # SMALL SHIFT NEEDED!
  
  if False: # TESTS
    print('tests')
    #data = np.zeros(1000)
    #data[500] = 1
    #from scipy import signal
    #y = signal.Gaussian_1D(200, std=sigma)
    #l = 50  
      #w = 15
    #y = np.zeros((l))
    #y[l / 2 - w - 1 : l / 2 + w] = 1
    #plt.plot(y)
    #plt.plot(data)
    #plt.plot(x, y)
    #plt.show()
    quit()
  
  #Mode same returns output of length max(M, N). Boundary effects are still visible.
  conv_array = np.convolve(data, y, mode = 'same')   
  return conv_array


# -------------------------------------------------------------------------------


def Filter1D_Gauss_Norm(data, sigma, **kwargs):
  # FILTER data WITH GAUSSIAN FILTER OF VARIANCE sigma
  from scipy import ndimage as nd # for Gaussian_1D filter
  
  return nd.filters.Gaussian_1D_filter1d(data, sigma) * sqrt(2 * pi) # last term necessary provided Python uses unitary FT


# -------------------------------------------------------------------------------
# FOURIER ANALYSIS
# -------------------------------------------------------------------------------


def DFT(**kwargs):
  """
  Discrete Fourier transform.
  
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
  this_func = this_lib + 'DFT: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')



  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------
# SPECIAL FUNCIONS
# -------------------------------------------------------------------------------


def Bessel_m01_FW3D(x, **kwargs):
  """
  Fullwave3D's approximation of the modified zero-order Bessel's function of the first kind.
  
  Parameters
  ----------
  x : float
    Argument of the function.
    
  Returns
  -------
  s : float
    Value of the function. Named identical to fullwave3D.
  
  Notes
  -----
  From the original source code: 
  'accuracy is <1/4%, which is quite adequate for distributed source'.
  There are discrepancies compared to Python's built-in function
  but these are for x < 0 which does not matter in Kaiser window
  (we look only at x > 0). FIXME: d-c.
  
  """
  this_func = this_lib + 'Bessel_m01_FW3D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  v = x * 0.5
  a = 1.0
  s = 1.0
  i = 0
  while a > 0.03:
    i = i + 1
    a = a * (v / i)
    s = s + a ** 2

  #print this_func + 'END'
  return s


# -------------------------------------------------------------------------------
# GAUSSIAN DISTRIBUTION
# -------------------------------------------------------------------------------


def Gaussian_2D_Unnorm(x, y, mu_x=0, mu_y=0, sig_x=1, sig_y=1, **kwargs):
  """
  Value of a unnormalized Gaussian 2D distribution 
  in a point. 
  
  Parameters
  ----------
  x, y : float
    Point to estimate value in.
  mu_x, mu_y : float 
    Mean of the distribution.
    Default: 0.
  sig_x, sig_y : float 
    Standard deviaton of the distribution.
    Default: 1.
    
  Returns
  -------
  value : float 
    Value of the function.
  
  Notes
  -----
  It can take values 0-1 only.
  
  """  
  this_func = this_lib + 'Gaussian_2D_Unnorm: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')    
  
  value = Gaussian_1D_Unnorm(x, mu_x, sig_x) * Gaussian_1D_Unnorm(y, mu_y, sig_y)
  
  #print this_func + 'END'
  return value


# -------------------------------------------------------------------------------


def Gaussian_1D(x, mu=0, sigma=1, **kwargs):
  """
  Value of a normalized Gaussian 1D distribution 
  in a point. 
  
  Parameters
  ----------
  x : float
    Point to estimate value in.
  mu : float 
    Mean of the distribution.
    Default: 0.
  sigma : float 
    Standard deviaton of the distribution.
    Default: 1.
    
  Returns
  -------
  value : float 
    Value of the function.
  
  Notes
  -----
  If only 2 args are provided 
  they are assumed to be (x, mu).
  
  """
  this_func = this_lib + 'Gaussian_1D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  value = 1 / (sigma * np.sqrt(2 * np.pi)) * Gaussian_1D_Unnorm(x, mu, sigma)
  
  #print this_func + 'END'
  return value


# -------------------------------------------------------------------------------


def Gaussian_1D_Unnorm(x, mu=0, sigma=1, **kwargs):
  """
  Value of a unnormalized Gaussian 1D distribution 
  in a point. 
  
  Parameters
  ----------
  x : float
    Point to estimate value in.
  mu : float 
    Mean of the distribution.
    Default: 0.
  sigma : float 
    Standard deviaton of the distribution.
    Default: 1.
    
  Returns
  -------
  value : float 
    Value of the function.
  
  Notes
  -----
  It can take values 0-1 only.
  
  """  
  this_func = this_lib + 'Gaussian_1D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  mu = float(mu) # OTHERWISE BLOTCHY
  sigma = float(sigma) ##
  
  value = np.exp(-((x - mu)**2) / (2 * sigma**2)) 
  
  #print this_func + 'END'
  return value  
 
 
# -------------------------------------------------------------------------------


def fwhm2std(fwhm, **kwargs):
  """
  Convert full width at half maximum (FWHM) 
  to a standard deviation of the Gaussian.
  
  Parameters
  ----------
  fwhm : float 
    As above.
    
  Returns
  -------
  value : float 
    Value of the function.
  
  Notes
  -----
  
  """  
  this_func = this_lib + 'fwhm2std: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  value = fwhm / (2. * np.sqrt(2 * np.log(2)))
  
  #print this_func + 'END'
  return value


# -------------------------------------------------------------------------------


def std2fwhm(std, **kwargs):
  """
  Convert a standard deviation of the Gaussian (std)
  to a full width at half maximum (FWHM).

  Parameters
  ----------
  std : float 
    As above.
    
  Returns
  -------
  value : float 
    Value of the function.
  
  Notes
  -----
  
  """  
  this_func = this_lib + 'fwhm2std: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  value = std * (2. * np.sqrt(2 * np.log(2)))
  
  #print this_func + 'END'
  return value  


# -------------------------------------------------------------------------------
# VARIA
# -------------------------------------------------------------------------------


def Signal_Modify_Simply(values, **kwargs):
  """
  
  
  Parameters
  ----------
  
  
  
  Returns
  -------
  
  
  Notes
  -----
  Order of scale and shift!!!!
  
  """
  this_func = this_lib + 'Signal_Modify_Simply: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  try:
    flagged = kwargs['flagged']
  except KeyError:
    flagged = []

  try:
    scale = kwargs['scale']
  except KeyError:
    scale = 1
    
  try:
    shift = kwargs['shift']
  except KeyError:
    shift = 0    
  
  nvalues = []
  
  for val in values:
    if (len(flagged) != 0) and (val == flagged[0]):
      nvalues.append(flagged[1])
      continue

    val = int(val * scale) # NOTE: IMPLICIT ROUNDING UP!!!  
    val += shift 

    nvalues.append(val)
    
  #print this_func + 'END'
  return nvalues


# -------------------------------------------------------------------------------

