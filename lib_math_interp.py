"""
Author: Kajetan Chrapkiewicz, 2018. 
Copywright: Ask for permision writing to K.Chrapkiewicz17@imperial.ac.uk.

This library provides framework procedures for interpolation.
At the moment fully implemented is:

1. 2D polynomial interpolation in Lagrange's, Newton's and Neville's form.
2. 2D / 3D sinc interpolation.

"""

## MODULES
import numpy as np
import matplotlib.pyplot as plt

## MY MODULES
from lib_generic import eprint
from lib_generic_CONST import epsi
##


## CONSTS
this_lib = 'lib_math_interp.py/'

## FUNCTIONS


#-------------------------------------------------------------------------------------


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


#-------------------------------------------------------------------------------------


def Interp_RBF_2d(x0, y0, X, Y, ZZ, **kwargs):
  """
  Find a value in f(x0, y0) given the array
  ZZ of discretized values of f on a 
  meshgrid(X, Y).
  
  Parameters
  ----------
  x0, y0 : float 
    Coordinates of a point to interpolate in.
  ZZ : array 
    2D array of function values evaluated on
    an original grid.
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Interp_RBF_2d: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  wx, wy = Interp_RBF_2d_Weights(x0, y0, X, Y)
  z0 = Interp_2d(wy, ZZ, wx)

  #print this_func + 'END'
  return z0


#-------------------------------------------------------------------------------------


def Interp_RBF_2d_Weights(x0, y0, X, Y, **kwargs):
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
  this_func = this_lib + 'Interp_Linear_2d_Weights: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  p = 2
  
  weights_x = np.array([1. / (abs(x - x0) + epsi)**2 for x in X])
  weights_y = np.array([1. / (abs(y - y0) + epsi)**2 for y in Y])

  wx = weights_x / sum(weights_x)
  wy = weights_y / sum(weights_y)

  #print this_func + 'END'
  return wx, wy


#-------------------------------------------------------------------------------------


def Interp_2d(wy, Z, wx, **kwargs):
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
  this_func = this_lib + 'Interp2d: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
    
  # MATRIX ELEMENT {wy|Z|wx}
  #wyZwx = np.dot(wy, np.dot(Z, wx))
  
  # SHOULD BE EQUIVALENT TO THE ABOVE (TESTING)
  wyZwx = 0.
  Zwx = np.zeros(len(Z))
  for i in range(len(Z)):
    for j in range(len(Z[0])):
      Zwx[i] += Z[i][j] * wx[j]
    wyZwx += wy[i] * Zwx[i]
      
  #print this_func + 'END'
  return wyZwx


#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------


def Interpolate(weights, data_values):
  """
  Generic interpolator applying a weighted-average interpolation
  given weights and values at data points.
  
  Parameters
  ----------
  weights : list / array
    List of weights that can as well be looked at as 
    interpolation coefficients.
  data_values : list / array
    List of values at data points. Each value must 
    correspond to the respective element of 'weights' list.
    In particular they need to be of the same length.

  Returns
  -------
  value : float 
    Value of the function in the point in which 
    the interpolation is performed.
  
  Notes
  -----
  It might as well be named weighted average.
  Note that coordinates do not to be provided 
  explicitly. They are included implicitly in weights.
  
  """
  this_func = this_lib + 'Interpolate: '
  #print this_func, 'START'
  
  value = 0.
  
  for i in range(len(weights)):
    value += weights[i] * data_values[i]
  
  #print this_func, 'END'
  return value


#-------------------------------------------------------------------------------------


def Interpolate_Bilinear(x, y, Q_corners):
  """
  Bilinear interpolation of the z-coordinate at (x,y)-point lying between Q points 
  which are corners of a square (not rectangular) area.
  
  Parameters
  ----------
  x : float
    X-coordinate of the point to interpolate in.
  y : float
    Y-coordinate of the point to interpolate in.
  Q_corners : list 
    List of 4 distinct 3D vectors (arrays / lists / tuples)
    being corners of a square (!) when projected on a XY-plane.
    They must be ordered in the following way: 
    first bottom-left and then clockwise, i.e. Q11, Q12, Q22, Q21.
  
  Returns 
  -------
  z : float
    Interpolated Z-coordinate of the point.
    
  Notes
  -----
  Formula is taken from wikipedia and it is tested 
  using surface plots. It resembles Python's cubic interpolation.
  
  Corners must be ordered in 
  FIXME: CHECK ORDER OF CORNERS COMPARED TO Neares_Bilinear_Corners
  FIXME: does it work outside the tile? (it is linear after all)
  FIXME: does it need to be square? It does, at least to have 
  quite simple analytical formulas.
  
  """
  this_func = this_lib + 'Interpolate_Bilinear: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  Q11, Q12, Q22, Q21 = Q_corners
  x1, y1, fQ11 = Q11
  x1, y2, fQ12 = Q12  
  x2, y2, fQ22 = Q22
  x2, y1, fQ21 = Q21
  
  if ((abs(x2 - x1) < epsi) or (abs(y2 - y1) < epsi)):
    eprint(this_func + 'Error. Coordinates of the corners must differ.\n')
    eprint(this_func + 'x1, x2, y1, y2: ' + str(x1) + ', ' + str(x2) + ', ' + str(y1) + ', ' + str(y2) + '\n')
    quit()
    
  #if ((x < x1) or (x > x2) or (y < y1) or (y > y2)):
  #  eprint(this_func + 'Point (' + str(x) + ',' + str(y) + ') outside the interpolation corners\n')
  
  z = (1. / ((x2 - x1) * (y2 - y1))) * ((y2 - y) * ((x2 - x) * fQ11 + (x - x1) * fQ21) + (y - y1) * ((x2 - x) * fQ12 + (x - x1) * fQ22))
  
  #print this_func, 'Interpolated value: ', z
  #print this_func + 'END'
  return z


#-------------------------------------------------------------------------------------


def Unpack_Args_Bilinear_Interp(m, m_ref, Q_corners):
  """
  Prepare parameters for exact computation of gradient and 
  Hessian for a bilinear-surface problem.
  
  Parameters
  ----------
  m : array 
    Point in the model space for which to calcalute the Hessian.
  m_ref : array
    Point in the model space of the same dimension as m. It is a reference
    from which the distance to the surface is measured.
  Q_corners : list
    List of points, each of the same dimensions as m. They are corners of 
    a square inside which bilinear interpolation is performed.
  
  Returns
  -------
  x, y, z : float
    Coordinates of the model.
  dx, dy : float 
    Intervals between bilinear-interpolation data-points (Q_corners).
    Not to confuse with dx, dy used in FD calculations.
  dx_G, dy_G, dz_G : float
    Coordinate-distances between model m and reference model m_ref (aka G).
  x1, x2 : float
    X-coordinates of Q_corners.
  y1, y2 : float
    Y-coordinates of Q_corners.
  fQ11, fQ12, fQ22, fQ21 : float
    Z-coordinates of Q_corners.
  
  Notes
  -----
  Qij = (xi, yj, fQij)
  
  """
  this_func = this_lib + 'Unpack_Args_Bilinear_Interp: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  x, y, z = m
  x_ref, y_ref, z_ref = m_ref
  Q11, Q12, Q22, Q21 = Q_corners

  x1, y1, fQ11 = Q11
  x1, y2, fQ12 = Q12
  x2, y2, fQ22 = Q22
  x2, y1, fQ21 = Q21
  
  dx, dy = [x2 - x1, y2 - y1]
  dx_G, dy_G, dz_G = [x - x_ref, y - y_ref, z - z_ref]
  
  #print this_func, 'x1, x2, y1, y2', x1, x2, y1, y2
  #print this_func, 'fQ11, fQ12, fQ22, fQ21', fQ11, fQ12, fQ22, fQ21
  #print this_func, 'x, y, z of the m0: ', x, y, z
  #print this_func, 'dx, dy', dx, dy
  #print this_func, 'dx_G, dy_G, dz_G', dx_G, dy_G, dz_G
  
  #print this_func + 'END'
  return x, y, z, dx, dy, dx_G, dy_G, dz_G, x1, x2, y1, y2, fQ11, fQ12, fQ22, fQ21


#-------------------------------------------------------------------------------------


def Model_Is_Inside_Bilinear_Tile(m, x1, x2, y1, y2):
  """
  Check if the model is outside the area of the bilinear interpolation.
  
  Parameters 
  ----------
  m : array
    Model which is supposed to be a 3D vector.
  x1, x2, y1, y2 : float
    Coordinates of respective corners of the (square) bilinear interpolation area.
    
  Returns 
  -------
  answer : bool
    True if projection of m onto XY plane lies inside the square defined by x1, x2, y1, y2.
    False otherwise.
    
  Notes
  -----
  NOTE: It allows for edges of the square.
  FIXME: Double-check it works in interpol. and optim.
  FIXME: Make it general as is_inside_cube???
  
  """
  this_func = this_lib + 'Model_Is_Inside_Bilinear_Tile: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  xm, ym = m[0], m[1]
  
  if ((xm < x1 or xm > x2) or (ym < y1) or (ym > y2)):
    answer = False
    #eprint(this_func + 'Model is outside the bilinear tile\n')
  
  else:
    answer = True
  
  #print this_func + 'END'
  return answer


#-------------------------------------------------------------------------------------


def Kaiser_Sinc_Weights(b_kaiser, r_kaiser, bessel, data_x, x0):
  """
  Calculate weights of data points for 
  the Kaiser-windowed sinc interpolation.
  
  Parameters
  ----------
  b_kaiser : float
    Shape of window
  r_kaiser : float
    Radius of window
  bessel : str 
    Implementation to choose: fw3D or python built-in.
  data_x : list
    List of X-coordinates of data points measured in arbitrary
    frame of reference.
  x0 : float
    Coordinate of a point to interpolate in the same frame
    of reference as data_x.
  
  Returns
  -------
  weights : list 
    List of weights of data points for 
    the Kaiser-windowed sinc interpolation in 
    the same order as in data_x list.
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Kaiser_Sinc_Weights: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  from lib_math_signal import Kaiser_Sinc_1d
  
  weights = []
  for i in range(len(data_x)):
    x = data_x[i]
    xdist = abs(x - x0)
    w = Kaiser_Sinc_1d(b_kaiser, r_kaiser, bessel, xdist)
    weights.append(w)
  
  #print this_func + 'END'
  return weights


#-------------------------------------------------------------------------------------


def Lagrange_Polyn_Weights(data_x, x0):
  """
  Calculate weights for 1D polynomial interpolation in Lagrange's form.
  
  Parameters
  ----------
  data_x : list
    List of X-coordinates of data points measured in arbitrary
    frame of reference.
  x0 : float
    Coordinate of a point to interpolate in the same frame
    of reference as data_x.
  
  Returns
  -------
  weights : list 
    List of weights of data points for 
    the interpolation.  
  
  Notes
  -----
  It allows x0 to be one of the data points xm; in this case lj(x0) = kronecker(j, m).  
    
  """
  this_func = this_lib + 'Lagrange_Polyn_Weights: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  weights = [] # COLLECTION OF BASE POLYNOMIALS lj(x) EVALUATED AT POINT TO INTERPOLATE IN: lj(x0)
  for j in range(len(data_x)):
    xj = data_x[j]
    lj = 1. # INITIAL VALUE J-TH BASE POLYNOMIAL lj(x0) 
    for m in range(len(data_x)):
      xm = data_x[m]
      if j != m:
        lj =  lj * (x0 - xm) / (xj - xm) # NOTE: WE CAN'T USE PYTHONIC *= HERE (> 1 TERM)
    weights.append(lj)
  
  #print this_func + 'END'
  return weights


#-------------------------------------------------------------------------------------


def Inv_Dist_Weights(mode, p, data_x, x0):
  """
  Calculate weights for 1D inverse-distance-weighting (IDW) interpolation.
  
  Parameters
  ----------
  mode : TEMPORARY FIXME
    
  p : float 
    Power to which the distance is raised (recip = 1 / d^p).
  data_x : list
    List of X-coordinates of data points measured in arbitrary
    frame of reference.
  x0 : float
    Coordinate of a point to interpolate in the same frame
    of reference as data_x.
  
  Returns
  -------
  weights : list 
    List of weights of data points for 
    the interpolation.  
  
  Notes
  -----
  It is one of the RBF interpolation methods.
  For details see wiki/Inverse_distance_weighting. 
  Here p=1 is prefered.
    
  """
  this_func = this_lib + 'Inv_Dist_Weights: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from numpy.linalg import norm
  
  reciprocals = []
  reciprocals_sum = 0.
  for i in range(len(data_x)):
    x = data_x[i]
    if mode == 'scalar':
      xdist = abs(x - x0)
    
    elif mode == 'vector': # HERE WE ASSUME x, x0 TO BE 2D np.array 
      xdist = norm(x - x0)
    
    #if xdist > 0:
    recip = 1. / ((xdist + epsi) ** p)
    #else:
      #recip = 0.
      #print 'recip=0 for x=', x
    
    reciprocals.append(recip)
    reciprocals_sum += recip
  
  weights = [] 
  for recip in reciprocals:
    w = recip / reciprocals_sum
    weights.append(w)
  
  #print 'weights[0]', weights[0]
  #print 'weights[-1]', weights[-1]
  
  #print this_func + 'END'
  return weights


#-------------------------------------------------------------------------------------












def Interpolate_Neville(datax, datay, x0):
  """
  Neville's algorithm of polynomial interpolation at a single point.
  
  Parameters
  ----------
  datax : list / array 
    x-coordinates of data points [x, y(x)].
  datay : list / array 
    y-coordinates of data points [x, y(x)]; 'abcissas'?
  x0 : float
    point in which to interpolate.

  Returns
  -------
  y0 : float
    P(x0), i.e. a value of polynomial in x0.
  
  Notes
  -----
  It does NOT output coefficients for data points.
  It is probably the best option when 
  For this see Netwon's algorithm of polynomial interpolation.
  
  """
  this_func = this_lib + 'Interpolate_Neville: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  n = len(datax) - 1
  
  datap = []
  for y in datay:
    datap.append(y)
  
  for k in range(1, n + 1):
    for i in np.linspace(n, k, n-k+1):
      i = int(i)
      datap[i] = datap[i] + (x0 - datax[i]) * (datap[i] - datap[i - 1]) / (datax[i] - datax[i - k])
  
  y0 = datap[n]
  
  #print this_func + 'END'  
  return y0 

def Newton_Weights(datax, datay):
  """
  
  Parameters
  ----------
  datax : list / array 
    x-coordinates of data points [x, y(x)].
  datay : list / array 
    y-coordinates of data points [x, y(x)]; 'abcissas'?

  Returns
  -------
  coeffs : list #FIXME
    Interpolation weights of data points.
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Newton_Coeffs: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  m = len(datax)
  n = m - 1
  
  
  coeffs = []
  for y in datay:
    coeffs.append(y)
  
  for k in range(1, m):
    for i in range(k, m): #np.linspace(k, n + 1, n + 1 - k):
      coeffs[i]  = (coeffs[i] - coeffs[k - 1]) / (datax[i] - datax[k - 1])
  
  #print this_func + 'END'    
  return coeffs  

def Interpolate_Newton(coeffs, datax, x0):
  """
  
  Parameters
  ----------
  datax : list / array 
    x-coordinates of data points [x, y(x)].

  Returns
  -------
  y0 : float
    y0 := P(x0)
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Interpolate_Newton: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  n = len(datax) - 1
  y0 = coeffs[-1]
  for k in range(1, n + 1):
    #print 'k', k
    y0 = coeffs[n - k] + y0 * (x0 - datax[n-k])
    #print 'y0', y0
    
  #print this_func + 'END'
  return y0

def Sinc_Rect_Coeffs(x0, dx, data):
    """
    We assume that original function was discretized
    on a grid dx, and the data is a fragment of this 
    grid windowed around x0 (FIXME?)
    CAUTION: data points must be spaced dx, not more sparsely!!!
    
    """
    coeffs = []
    for n in range(len(data)):
      x = data[n]
      c = np.sinc((x0 - x) / dx)
      coeffs.append(c)
    return coeffs





def Interpolate_Kaiser_Sinc_2D(b_kaiser, r_kaiser, bessel, data, x0, y0):
  z0 = 0
  for d in data:
    x = d[0]
    y = d[1]
    z = d[2]
    z0 += z * Kaiser_Sinc_2D(b_kaiser, r_kaiser, bessel, abs(x - x0), abs(y - y0))

  #print this_func + 'END'
  return z0



