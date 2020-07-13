"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permision writing to k.chrapkiewicz17@imperial.ac.uk.

At the moment the module contains implementations of:

  1. simple mathematical operations such as Euclidan distance, 
     radian-degree conversion etc.

"""


## MODULES
import numpy as np
from numpy.linalg import norm
from numpy import pi as PI

# MY OWN ONES
from lib_generic import eprint, Kwarg
from lib_generic_CONST import err_value, epsi, verbos_func


## CONSTANTS
this_lib = 'lib_math_generic.py/'


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
  return 


# -------------------------------------------------------------------------------
# DIFFERENTIAL GEOMETRY
# -------------------------------------------------------------------------------


def Polar2Cartes(r, phi, **kwargs):
  """
  Convert 2D polar coordinates 
  into Cartesian ones.
  
  Parameters
  ----------
  r : float 
    Radius.
  phi : float
    Angle (rad).
  
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  x, y : float 
    Cartesian.
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Polar2Cartes: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  x = r * np.cos(phi)
  y = r * np.sin(phi)
  
  #print this_func + 'END'
  return x, y


# -------------------------------------------------------------------------------
# MEANS
# -------------------------------------------------------------------------------


def RMS(vector, **kwargs):
  """
  Root mean square.
  
  Parameters
  ----------
  vector : list / array
    List of numbers.
  
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  rms : float 
    RMS of the vector.
  
  Notes
  -----
  rms(y) = (y_i * y_i / n)^(1/2)
  
  """
  this_func = this_lib + 'RMS: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  n = len(vector)
  rms = np.sqrt(sum(np.array([i**2 for i in vector])) / n)

  #print this_func + 'END'
  return rms


# -------------------------------------------------------------------------------
# COMMON FUNCTIONS
# -------------------------------------------------------------------------------


def Sine_2D(x, y, **kwargs):
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
  this_func = this_lib + 'Sine_2D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  k_xy = Kwarg('k_xy', [1, 1], kwargs)  
  phi_xy = Kwarg('phi_xy', [0, 0], kwargs)
  #zmin = Kwarg('zmin', None, kwargs)
  #zmax = Kwarg('zmax', None, kwargs)
  z0 = Kwarg('z0', 0, kwargs)
  A = Kwarg('A', 1, kwargs)
  
  kx, ky = k_xy
  phi_x, phi_y = phi_xy
  
  zx = Sine(x, k=kx, phi=phi_x, y0=0, A=1)
  zy = Sine(y, k=ky, phi=phi_y, y0=0, A=1) # FIXME: CLEARLY WRONG
  
  
  print(this_func, 'zx, zy', zx, zy)
  
  z = z0 + A * zx * zy
  
  #print this_func + 'END'
  return z


# -------------------------------------------------------------------------------


def Sine(x, **kwargs):
  """
  Sine function.
  
  
  Parameters
  ----------
  x : float 
    Time (not angle!)
   
  **kwargs : keyword arguments (optional)
    Current capabilities:
    k : float 
      Wave number.
    phi : float 
      Phase.
    
   
  
  Returns
  -------
  y : float
    Value of the function evaluated at x.
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Sine: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  k = Kwarg('k', 1, kwargs)
  phi = Kwarg('phi', 0, kwargs)
  ymin = Kwarg('ymin', None, kwargs)
  ymax = Kwarg('ymax', None, kwargs)
  y0 = Kwarg('y0', 'None', kwargs)
  A = Kwarg('A', None, kwargs)
  
  y = np.sin((2*np.pi*x / k) + phi) 
  
  if y0 == 'None': # BECAUSE WE MAY WANT y0=0 WHICH IS TREATED AS NONE
    if ymin and ymax:
      y0 = (ymin + ymax) / 2.
    else:
      raise ValueError('You need to provide either (ymin, ymax) or (y0, A)')
    
  
  if not A:
    if ymin and ymax:
      A = (ymax - ymin) / 2.
    else:
      raise ValueError('You need to provide either (ymin, ymax) or (y0, A)')
  
  y = y0 + A * y
  
  #print this_func + 'END'
  return y


# -------------------------------------------------------------------------------


def Sine_2D_Old(x, y, lamb_x, lamb_y, phi_x, phi_y, z_min, z_max, **kwargs):
  """
  
  Parameters
  ----------
  x : float
    Coordinate along X-axis.
  y : float
    Coordinate along Y-axis.
  lamb_x
    Wavelength
    
  phi_x 
    Phase shift.
    
  z_min 
    Min value.
  
  
  
  Returns
  -------
  value : float 
    Value of the function.
    
  Notes
  -----
  Set phi = np.pi / 2 to get a cosine.
  
  Amplitude (A) is measured from z=0. Peak-to-peak amplitude = 2A.
  
  """
  this_func = this_lib + 'Sine_2D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  z0 = (z_min + z_max) / 2.
  A = (z_max - z_min) / 2.
  
  value_x = np.sin(phi_x + 2 * PI * x / lamb_x) 
  value_y = np.sin(phi_y + 2 * PI * y / lamb_y)
  
  value = z0 + A * value_x * value_y
  
  #print this_func + 'END'
  return value 
  

# -------------------------------------------------------------------------------


def Find_Linear_Function(p1, p2):
  """
  Find a and b coefficients of 
  y = ax + b given to points (x,y).
  
  Parameters
  ----------
  p1, p2 : array 
    2 arrays (x,y).
    
  Returns
  -------
  a, b : float 
    Coefficents of y = ax + b
    equation.
    
  Notes
  -----
  
  """
  this_func = this_lib + 'Find_Linear_Function: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  x1, y1 = p1 
  x2, y2 = p2
  dy = y2 - y1
  dx = x2 - x1
  a = dy / float(dx)
  b = y1 - x1 * dy / float(dx)
  
  #print this_func + 'END'
  return a, b


# -------------------------------------------------------------------------------


def Straight_Line(X, alpha, y0):
  """
  Find y-values of a linear function y = ax + b
  given X coordinates, slope angle and
  y(X[0]).
  
  Parameters
  ----------
  X : list
    List of X coordinates.
  alpha : float 
    Slope (angle) of the function
    in degrees.
  y0 : float 
    Value at X[0].
    
  Returns
  -------
  Y : list 
    List of y-values.
    
  Notes
  -----
  
  """
  this_func = this_lib + 'Straight_Line: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')    
  
  x0 = X[0]
  
  a = np.tan(deg2rad(alpha))
  b = y0 - a * x0
  
  Y = [a * x + b for x in X]
  
  #print this_func + 'END'
  return Y


# -------------------------------------------------------------------------------
# TENSOR OPERATIONS
# -------------------------------------------------------------------------------


def Invert_Matrix_2x2(M):
  """
  Invert 2x2 matrix.
  
  Parameters
  ----------
  M : array
    2x2 matrix.
    
  Returns
  -------
  M_inv : array
    Inverse of M.
    
  Notes
  -----
  It uses a well-known, simple formula:
  A^(-1) = ((a, b), (c, d))^(-1) = ((d, -b), (-c, a)) / det(A)
  Determinant is also computed analytically taking into account
  values close to zero.
  
  """
  this_func = this_lib + 'Invert_Matrix_2x2: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  a = M[0][0]
  b = M[0][1]
  c = M[1][0]
  d = M[1][1]
  
  det = float(a * d - b *c)
  
  if det < epsi:
    eprint(this_func + 'Error. Determinant smaller than epsi, det = ' + str(det) + ' No inverse can be calculated\n')
    quit()
  
  M_inv_raw = np.array(((d, -b), (-c, a)))
  M_inv = M_inv_raw / det 
   
  #print this_func + 'END'
  return M_inv


# -------------------------------------------------------------------------------


def Normed(v):
  """
  Normalize a vector.
  
  Parameters
  ----------
  v : array
    Vector to normalize.
    
  Returns
  -------
  nv : array
    Normalized vector v.
    
  Notes
  -----
  Raised error if norm -> 0.
  
  NOTE: IT CORRESPONDS TO RMS 
  NORMALIZATION, DOESN'T IT?
  
  """
  this_func = this_lib + 'Normed: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  v = np.array(v)
  
  if norm(v) < epsi:
    eprint(this_func + 'Error. Norm approx. 0\n')
    quit()
  
  else:
    nv = v / norm(v)
    
  #print this_func + 'END'
  return nv


# -------------------------------------------------------------------------------
# GEOMETRY
# -------------------------------------------------------------------------------


def Rotate_2D(origin, point, angle):
  """
  Rotate a point clockwise around the origin.
  
  Parameters
  ----------
  origin : array 
    2D axis of rotation.
  point : array
    Point to rotate.
  angle : float 
    Angle of rotation (clockwise)
    in degrees.
    
  Returns
  -------
  qx 
  qy
  qz
    
  Notes
  -----
  3D WITH ROT. AXIS ALONG Y DIRECTION
  
  """
  this_func = this_lib + 'Rotate_2D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')   
  
  ox, oz = origin
  px, py, pz = point
  angle = deg2rad(angle)
  
  qx = ox + np.cos(angle) * (px - ox) - np.sin(angle) * (pz - oz)
  qy = py # FIXME
  qz = oz + np.sin(angle) * (px - ox) + np.cos(angle) * (pz - oz)
  
  #print this_func + 'END'  
  return qx, qy, qz


# -------------------------------------------------------------------------------


def Dist(point1, point2, **kwargs):
  """
  Calculate Euclidean (L2) 
  distance between 2 points 
  of a dimension N.
  
  Parameters
  ----------
  point1, point2 : array 
    2 arrays of coordinates.
    
  Returns
  -------
  dist : float 
    Distance between the points.
    
  Notes
  -----
  
  """
  this_func = this_lib + 'Dist: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START') 
  
  if len(point1) != len(point2):
    eprint(this_func + 'Error. Points have different dimensions. \n')
    quit()
  else:
    ndims = len(point1)
  
  squares_sum = 0
  for i in range(ndims):
    squares_sum += (point1[i] - point2[i]) ** 2
    
  dist = np.sqrt(squares_sum)
  
  #print this_func + 'END'
  return dist


# -------------------------------------------------------------------------------


def deg2rad(angle):
  """
  Convert degrees to radians.
  
  Parameters
  ----------
  angle : float 
    Angle in degrees.
    
  Returns
  -------
  angle_rad : float 
    Angle in radians.
    
  Notes
  -----
  
  """
  this_func = this_lib + 'deg2rad: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  angle_rad = angle * np.pi / 180
  
  
  #print this_func + 'END'
  return angle_rad


# -------------------------------------------------------------------------------


def rad2deg(angle):
  """
  Convert radians to degrees.
  
  Parameters
  ----------
  angle : float 
    Angle in radians.
    
  Returns
  -------
  angle_deg : float 
    Angle in degrees.
    
  Notes
  -----
  
  """
  this_func = this_lib + 'rad2deg: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  angle_deg = angle * 180 / np.pi 
  
  
  #print this_func + 'END'
  return angle_deg


# -------------------------------------------------------------------------------
# SURFACES
# -------------------------------------------------------------------------------


def Closest_Surface_Point(G, **kwargs):
  """
  Find a point on a surface that minimizes 
  distance from a reference point (not on the surface).
  
  Parameters
  ----------
  G : array
    Reference point.
  **kwargs : keyword arguments, optional
      Just passing it down.
  
  Returns
  -------
  point : array 
    Nearest point.
  dist : float 
    Distance between G and point.
  
  Notes
  -----
  # FIXME: TO FILL
  
  """
  this_func = this_lib + 'Closest_Surface_Point: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  


  #print this_func + 'END'
  return


# -------------------------------------------------------------------------------


def Nearest_Bilinear_Corners(G, radius, overhang): 
  """
  Framework function to find nearest bilinear tiles
  (described as their corners) of the point G 
  within the radius.
  
  Parameters
  ----------
  G : array
    Point in the 'middle' of the bilin. tile.
  radius : int 
    Radius of the frame to find corners in.
  overhang : bool 
    Is overhanging allowed?
  
  Returns
  -------
  Qs : list 
    List of tiles (lists of corners).
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Nearest_Bilinear_Corners: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  if not overhang:
    Qs = Nearest_Bilinear_Corners_No_Overhang(G, radius)#, surface)
  
  else:
    eprint(this_func + 'Error. Overhangs not yet implemented\n')
    quit()

  #print this_func + 'END'
  return Qs


# -------------------------------------------------------------------------------


def Nearest_Bilinear_Corners_No_Overhang(G, radius):
  """
  Find nearest bilinear tiles
  (described as their corners) of the point G 
  within the radius.
  
  Parameters
  ----------
  G : array
    Point in the 'middle' of the bilin. tile.
  radius : int 
    Radius of the frame to find corners in.
  
  Returns
  -------
  Qs : list 
    List of tiles (lists of corners).
  
  Notes
  -----
  FIXME: Deal with model edges!
  
  """  
  this_func = this_lib + 'Nearest_Bilinear_Corners_No_Overhang: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  x_G, y_G = G[0], G[1]
  
  Qs = []
  midpoints = Find_Grid_Midpoints([x_G, y_G], radius)
  Qs = Find_Corners(midpoints)
  
  #print this_func + 'END'
  return Qs


# -------------------------------------------------------------------------------


def Find_Grid_Midpoints(frame_centre, radius):
  """
  Within a frame  
  find grid midpoints which are average
  of the grid-cell corners. 
  
  Parameters
  ----------
  frame_centre : array 
    Coordiantes of the frame centre [x,y].
  radius : int 
    Radius of the frame.
  
  Returns
  -------
  midpoints : list 
    List of midpoints.
  
  Notes
  -----
  Assumes dx = 1 grid nodes (no sub- nor super-sampling)
  
  """
  this_func = this_lib + 'Find_Grid_Midpoints: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  dx = 1

  x0 = frame_centre[0]
  y0 = frame_centre[1]
  
  #xrang = range(-radius, radius + 1)
  
  midpoints = []
  xy_range = np.arange(-radius + 0.5 * dx, radius + 0.5 *dx)
  
  for dx in xy_range:
    for dy in xy_range:
      #if dx == 0 or dy == 0:
        #continue 
      midpoints.append([x0 + dx, y0 + dy])
  
  #print this_func + 'END'
  return midpoints


# -------------------------------------------------------------------------------


def Find_Corners(midpoints):
  """
  Find corners surrounding 
  every midpoint of the grid.
  
  Parameters
  ----------
  midpoints : list 
    List of midpoints.
  
  Returns
  -------
  corners : list 
    List of corners
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Find_Corners: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  corners = []
  for mp in midpoints:
    mp_corners = []
    x = mp[0]
    y = mp[1]
    
    for dx in [-0.5, 0.5]:
      for dy in [-0.5, 0.5]:
        mp_corners.append([x+dx, y+dy])
    
    corners.append(mp_corners)

  #print this_func + 'END'
  return corners


# -------------------------------------------------------------------------------

