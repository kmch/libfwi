"""
Author: Kajetan Chrapkiewicz, 2018. 
Copywright: Ask for permision writing to K.Chrapkiewicz17@imperial.ac.uk.

This library provides framework procedures for a general, non-linear optimization.
At the moment fully implemented is:

1. (quasi)-Newton optimization of a distance from a (bi)linear surface in 3D.

"""


## MODULES
import numpy as np
import matplotlib.pyplot as plt

## MY MODULES
from lib_generic import eprint
#from constants import *
##


## CONSTS
this_lib = 'lib_math_optim.py/'

## FUNCTIONS


#-----------------------------------------------------------------------------------------------------------


def Optimize_Newton(problem, approx, threshold, niters, m0, args):
  """
  Perform an iterative (quasi-)Newton optimization of a general non-linear problem.
  
  Parameters
  ----------
  problem : str 
    Type of the problem to solve.
    Current capabilities:
      - 'lin2d' (minimum distance from a 2D surface interpolated bilinearly).
  approx: str
    Approximation of the Hessian matrix used in the optimization.
    Current capabilities:
      - 'full' - full Hessian.
      - 'diag' - diagonal Hessian.
      - 'alpha' - steepest descent (alpha stands for a step-length).
      - 'alpha_calc' - steepest descent with line-search? #FIXME
  threshold : float
    Stopping criterion expressed by: norm(current_model - previous_model).
  niters : int 
    Maximum no. of iterations to perform - additional stopping criterion
    if the algorithm does not converge.     
  m0 : array 
    Starting model for the optimization expressed as a vector of parameters.
  args : list
    Additional parameters problem-specific, passed on to inner functions. 
  
  Returns
  -------
  grads : list
    List of gradients for each iteration.
  Hs : list
      List of Hessians for each iteration.
  models : list
      List of current models for each iteration.
  
  Notes
  -----
  It is written in a low-level (non-Pythonic) style
  to retain portability to Fortran etc.
  
  """
  this_func = this_lib + 'Optimize_Newton: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  from numpy.linalg import norm
  
  # INITIALIZE
  models, grads, Hs = [], [], []
  #print this_func, 'Starting model: ', m0
  m = m0 # NOTE
  iter_no = 1

  # MAIN LOOP
  while True:
    #print this_func, 'Iteration ', iter_no, ' starts...\n'
    
    # CHECK CONSTRAINTS
    if Model_Satisfies_Constraints(problem, m, args):
      models.append(m)
    else:
      break
    
    # GRADIENT CALCULATION
    grad = Gradient(problem, m, args)
    grads.append(grad) 

    # HESSIAN CALCULATION
    H = Hessian(problem, approx, m, args)
    Hs.append(H)
    
    # MODEL UPDATE
    m_prev = m
    m = Update_Model(problem, m, grad, H, args)
    #print this_func, 'Previous model: ', m_prev
    #print this_func, 'Current model: ', m
    
    # STOPPING CRTIERION
    mnorm = norm(m - m_prev)
    #print this_func, 'Norm of the update: ', mnorm
    if (mnorm < threshold) or (iter_no > (niters - 1)): # -1 SINCE STARTED AT 0
      break
    
    # ITERATION DONE
    #print this_func, '...iteration ', iter_no, ' ends.\n'
    iter_no += 1 
  
  # FINAL REPORT
  if len(models) == 0: 
    eprint(this_func + 'Error. No models found because the starting model was already outside interpolation domain\n')
    quit()
  
  elif len(models) == 1:
    eprint(this_func + 'Final model is the starting model\n')
  
  #else: 
  #  print this_func, 'Final model: ', models[-1]
  #  print this_func, 'No. of iterations performed: ', iter_no # FIXME: or i - 1?
  
  #print this_func + 'END'
  return grads, Hs, models


#-----------------------------------------------------------------------------------------------------------


def Model_Satisfies_Constraints(problem, m, args):
  """
  Check if model satisfies the constraints defined for 
  a given optimization problem. 
  
  Parameters
  ----------
  problem : str 
    Type of the problem to solve.
    Current capabilities:
      - 'lin2d' : minimum distance from a 2D surface interpolated bilinearly.
  m : array 
    Point in the model space for which to calcalute the gradient.
  args : list
    Additional parameters. They are problem-specific and
    passed on to inner functions.   
  
  Returns
  -------
  answer : bool
    True if model satisfies the constraints. 
    False otherwise.
  
  Notes
  -----
  Example of the constraints: model must lie inside the square 
  defined for a bilinear interpolation.
  
  """
  this_func = this_lib + 'Model_Satisfies_Constraints: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_math_interp import Unpack_Args_Bilinear_Interp, Model_Is_Inside_Bilinear_Tile
  
  if problem == 'lin2d':
    m_ref = args[0]
    Q_corners = args[1: ]
    x, y, z, dx, dy, dx_G, dy_G, dz_G, x1, x2, y1, y2, fQ11, fQ12, fQ22, fQ21 = Unpack_Args_Bilinear_Interp(m, m_ref, Q_corners)
    answer = Model_Is_Inside_Bilinear_Tile(m, x1, x2, y1, y2)
    
  else: 
    eprint(this_func + 'Error. Unknown problem: ' + problem + '\n')
    quit()

  if answer == False:
    eprint(this_func + 'Current model does NOT satisfy the constraints of the optimization problem\n')
  
  #print this_func + 'END'
  return answer


# -----------------------------------------------------------------------------------------------------------


def Gradient(problem, m, args): 
  """
  Calculate a gradient of the error function for a given model.
  
  Parameters
  ----------
  problem : str 
    Type of the problem to solve.
    Current capabilities:
      - 'lin2d' : minimum distance from a 2D surface interpolated bilinearly.
  m : array 
    Point in the model space for which to calcalute the gradient.
  args : list
    Additional parameters. They are problem-specific and
    passed on to inner functions. 
  
  Returns
  -------
  grad : array
    Gradient vector. If any errors are raised, empty [] list is returned.
  
  Notes
  -----
  The returned gradient is NOT normed NOR multiplied by (-1).
  
  It is written in a low-level (non-Pythonic) style
  to remain easily portable to Fortran etc.
  
  """
  this_func = this_lib + 'Gradient: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  from numpy.linalg import norm  
  
  if problem == 'lin1d':
    grad = Gradient_Linear(m, args)
  
  elif problem == 'lin2d':
    m_ref = args[0]
    Q_corners = args[1: ]
    grad = Gradient_Bilinear(m, m_ref, Q_corners)
  
  else:
    eprint(this_func + 'Error. Problem: ' + problem + ' not yet implemented.\n')
    quit()
  
  if norm(grad) == 0:
    eprint(this_func + 'Zero norm of the gradient.\n')
  
  #print this_func, 'Gradient: ', grad
  
  #print this_func + 'END'
  return grad


#-----------------------------------------------------------------------------------------------------------


def Gradient_Bilinear(m, m_ref, Q_corners):
  """
  Calculate a gradient of the error function for a given model.
  The error function is a distance between a point on a 2D surface
  immersed in 3D space, and a reference point lying beyond the surface.
  The surface is defined on a grid nodes and interpolated bilinearly
  in between.
  
  Parameters
  ----------
  m : array 
    Point in the model space for which to calcalute the gradient.
  m_ref : array
    Point in the model space of the same dimension as m. It is a reference
    from which the distance to the surface is measured.
  Q_corners : list
    List of 4 distinct 3D vectors (arrays / lists / tuples)
    being corners of a square (!) when projected on a XY-plane.
    They must be ordered in the following way: 
    first bottom-left and then clockwise.
  
  Returns
  -------
  grad : array
    Gradient vector. If any errors are raised, empty [] list is returned.
  
  Notes
  -----
  One can calculate the gradient either using analytical formula
  or a central finite difference approximation.
  
  """
  this_func = this_lib + 'Gradient_Bilinear: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_math_interp import Interpolate_Bilinear, Unpack_Args_Bilinear_Interp
  
  x, y, z, dx, dy, dx_G, dy_G, dz_G, x1, x2, y1, y2, fQ11, fQ12, fQ22, fQ21 = Unpack_Args_Bilinear_Interp(m, m_ref, Q_corners)
  
  mode = 'num'
  if mode == 'num':
    grad = Gradient_Bilinear_Num(m, m_ref, Q_corners)
    
  elif mode == 'exact':
    # FIXME: SHOULDN'T IT BE 2 * dx + 2 * ... ?
    # WE ONLY WANT DIRECTION (WE NORM IT ANYWAY)
  
    df_dx = ((y2 - y) * (fQ21 - fQ11) + (y - y1) * (fQ22 - fQ12)) / float(dx * dy)
    
    df_dy = (x2 - x) * (fQ12 - fQ11) + (x - x1) * (fQ22 - fQ21) / float(dx * dy)
    
    d_phi_d_x = 2 * (dx_G + dz_G * df_dx)
    d_phi_d_y = 2 * (dy_G + dz_G * df_dy)
    
    gradx = d_phi_d_x
    grady = d_phi_d_y
    
    grad = np.array((gradx, grady))
    
    #gradx = dx + (dz / (x2 - x1) / (y2 - y1)) * ((y2 - y) * (fQ21 - fQ11) + (y - y1) * (fQ22 - fQ12)) 
    #grady = dy + (dz / (x2 - x1) / (y2 - y1)) * ((x2 - x) * (fQ12 - fQ11) + (x - x1) * (fQ22 - fQ21))  
  
  else:
    eprint(this_func + 'Error. Wrong mode: ' + mode + '\n')
    quit()
  
  #print this_func + 'END'
  return grad


#-----------------------------------------------------------------------------------------------------------


def Gradient_Bilinear_Num(m, m_ref, Q_corners):
  """
  Calculate a central-difference approximation
  of the Gradient_Bilinear.
  
  Parameters
  ----------
  m : array 
    Point in the model space for which to calcalute the gradient.
  m_ref : array
    Point in the model space of the same dimension as m. It is a reference
    point from which the distance to the surface is measured.
  Q_corners : list
    List of 4 distinct 3D vectors (arrays / lists / tuples)
    being corners of a square (!) when projected on a XY-plane.
    They must be ordered in the following way: 
    first bottom-left and then clockwise. 
  
  Returns
  -------
  grad : array
    Gradient vector. If any errors are raised, empty [] list is returned.
  
  Notes FIXME
  -----
  It is NOT written in a general form, using a loop over coordinates
  to remain explicit.
  
  """
  this_func = this_lib + 'Gradient_Bilinear_Num: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from numpy.linalg import norm
  from lib_math_interp import Interpolate_Bilinear
  from constants import differential
  
  x, y, z = m  
  #print this_func, ' x, y', x, y
  #for i, Q in enumerate(Q_corners):
  #  print this_func, 'Q no.', i + 1, ':', Q
  
  
  # FINITE DIFFERENCES
  dx = differential # FIXME: THE SMALLER THE BETTER?
  dy = dx
  
  
  # X-COORDINATES OF CENTRAL-DIFFERENCE
  x_plus_dx = x + dx
  x_minus_dx = x - dx
  
  # POINTS ON SURFACE AT X-COORDINATES
  z_x_plus_dx = Interpolate_Bilinear(x_plus_dx, y, Q_corners)
  z_x_minus_dx = Interpolate_Bilinear(x_minus_dx, y, Q_corners)
  
  r_x_plus_dx = np.array([x_plus_dx, y, z_x_plus_dx])
  r_x_minus_dx = np.array([x_minus_dx, y, z_x_minus_dx])
  
  # X COMPONENT OF THE FINITE DIFFERENCE OF THE ERROR FUNCTION
  phi_x_plus_dx = norm(m_ref - r_x_plus_dx)
  phi_x_minus_dx = norm(m_ref - r_x_minus_dx)
  
  # FD-APPROX. OF X-PARTIAL DERIVATIVE OF THE ERROR FUNCTION
  d_phi_dx = (phi_x_plus_dx - phi_x_minus_dx) / (2 * dx)
    
  
  # Y-COORDINATES OF CENTRAL-DIFFERENCE
  y_plus_dy = y + dy
  y_minus_dy = y - dy
  
  # POINTS ON SURFACE AT Y-COORDINATES
  z_y_plus_dy = Interpolate_Bilinear(x, y_plus_dy, Q_corners)
  z_y_minus_dy = Interpolate_Bilinear(x, y_minus_dy, Q_corners)
  
  r_y_plus_dy = np.array([x, y_plus_dy, z_y_plus_dy])
  r_y_minus_dy = np.array([x, y_minus_dy, z_y_minus_dy])
  
  # y COMPONENT OF THE FINITE DIFFERENCE OF THE ERROR FUNCTION
  phi_y_plus_dy = norm(m_ref - r_y_plus_dy)
  phi_y_minus_dy = norm(m_ref - r_y_minus_dy)

  # FD-APPROX. OF Y-PARTIAL DERIVATIVE OF THE ERROR FUNCTION
  d_phi_dy = (phi_y_plus_dy - phi_y_minus_dy) / (2 * dy)
   
   
  # GRADIENT 
  grad = np.array((d_phi_dx, d_phi_dy))
  
  #print this_func + 'END'
  return grad


#-----------------------------------------------------------------------------------------------------------


def Hessian(problem, approx, m, args):
  """
  Calculate a Hessian of the error function for a given model.
  
  Parameters
  ----------
  problem : str 
    Type of the problem to solve.
    Current capabilities:
      - 'lin2d' (minimum distance from a 2D surface interpolated bilinearly).
  approx : str 
    Approximation of the Hessian. 
    Current capabilities:
      - 'full' - full Hessian.
      - 'diag' - diagonal Hessian.
      - 'alpha' - steepest descent (alpha stands for a step-length).
      - 'alpha_calc' - steepest descent with line-search? #FIXME
  m : array 
    Point in the model space for which to calcalute the gradient.
  args : list
    Additional parameters; problem-specific; passed on to inner functions. 
  
  Returns
  -------
  H : array
    Hessian matrix. No. of dimensions depends on the problem.
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Hessian: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  if problem == 'lin1d':
    a = args[0]
    H = Hessian_Linear(a)
  
  elif problem == 'lin2d':
    m_ref = args[0]
    Q_corners = args[1: ]
    H = Hessian_Bilinear(approx, m, m_ref, Q_corners)  
  
  else:
    eprint(this_func + 'Error. Problem: ' + problem + ' not yet implemented\n')
    quit()
  
  #print this_func, 'Hessian: ', H
  
  #print this_func + 'END'
  return H


#-----------------------------------------------------------------------------------------------------------


def Hessian_Bilinear(approx, m, m_ref, Q_corners):
  """
  Calculate a Hessian of the error function for a given model.
  The error function is a distance between a point on a 2D surface
  immersed in 3D space, and a reference point lying beyond the surface.
  The surface is defined on a grid nodes and interpolated bilinearly
  in between.
  
  Parameters
  ----------
  approx : str 
    Approximation of the Hessian. 
    Current capabilities:
      - 'full' - full Hessian.
      - 'diag' - diagonal Hessian.
      - 'alpha' - steepest descent (alpha stands for a step-length).
      - 'alpha_calc' - steepest descent with line-search? #FIXME  
  m : array 
    Point in the model space for which to calcalute the Hessian.
  m_ref : array
    Point in the model space of the same dimension as m. It is a reference
    from which the distance to the surface is measured.
  Q_corners : list
    List of 4 distinct 3D vectors (arrays / lists / tuples)
    being corners of a square (!) when projected on a XY-plane.
    They must be ordered in the following way: 
    first bottom-left and then clockwise.
  
  Returns
  -------
  H : array
    Hessian matrix. No. of dimensions depends on the problem.
  
  Notes
  -----
  One can calculate the Hessian either using analytical formula
  or a central finite difference approximation.
  
  """
  this_func = this_lib + 'Hessian_Bilinear: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_math_interp import Unpack_Args_Bilinear_Interp
  
  mode = 'num'
  if mode == 'num':
    if approx != 'diag':
      eprint(this_func + 'Error. approx must be <diag> for a num. Hessian\n')
      quit()
      
    H11, H22 = Hessian_Bilinear_Num(m, m_ref, Q_corners)

  elif mode == 'exact':
    x, y, z, dx, dy, dx_G, dy_G, dz_G, x1, x2, y1, y2, fQ11, fQ12, fQ22, fQ21 = Unpack_Args_Bilinear_Interp(m, m_ref, Q_corners)
  
    #NOTE SHOULD BE THE SAME AS IN GRAD
    df_dx = ((y2 - y) * (fQ21 - fQ11) + (y - y1) * (fQ22 - fQ12)) / (dx * dy)
    
    df_dy = (x2 - x) * (fQ12 - fQ11) + (x - x1) * (fQ22 - fQ21) / (dx * dy)
    ##
    
    # DIAGONAL ELEMENTS
    d2_phi_dx2 = 2 * (1 + df_dx ** 2)
    d2_phi_dy2 = 2 * (1 + df_dy ** 2)
    H11 = d2_phi_dx2
    H22 = d2_phi_dy2
    ##
  
  else:
    eprint(this_func + 'Error. Wrong mode: ' + mode + ' \n')
    quit()


  # OFF-DIAGONAL ELEMENTS
  if approx == 'none':
    #print this_func, 'Optimization method: Netwon without any approximation'
    d2_f_dxdy = (fQ11 + fQ22 - fQ12 - fQ21) / (dx * dy)
    #print 'd2_f_dxdy', d2_f_dxdy
    d2_phi_dx_dy = 2 * (df_dx * df_dy + dz_G * d2_f_dxdy)
    d2_phi_dy_dx = d2_phi_dx_dy # SINCE d2_f_dydx = d2_f_dxdy
    H12 = d2_phi_dx_dy
    H21 = d2_phi_dy_dx
  
  elif approx == 'diag':
    #print this_func, 'Optimization method: Quasi-Netwon with diagonal Hessian'
    H12 = 0.
    H21 = 0.
   
  elif approx == 'alpha':
    alpha =  0.5
    #print this_func, 'Optimization method: gradient descent with arbitrary step-length equal to ', alpha
    H11 = alpha
    H22 = alpha
    H12 = 0.
    H21 = 0.  

  elif approx == 'alpha_calc':
    #print this_func, 'Optimization method: gradient descent with calculated step-length'
    dm =  0.000001 # WE WANT A SMALL CHANGE IN THE DIRECTION OF GRADIENT
    H11 = dm
    H22 = dm
    H12 = 0.
    H21 = 0.    
  
  else:
    eprint(this_func + 'Error. Wrong approx: ' + str(approx) + '\n') 
    quit()
  
  H = np.array(((H11, H12), (H21, H22)))
  
  #print this_func + 'END'
  return H


#-----------------------------------------------------------------------------------------------------------


def Hessian_Bilinear_Num(m, m_ref, Q_corners):
  """
  Calculate a central-difference approximation
  of the Gradient_Bilinear.
  
  Parameters
  ----------
  m : array 
    Point in the model space for which to calcalute the gradient.
  m_ref : array
    Point in the model space of the same dimension as m. It is a reference
    point from which the distance to the surface is measured.
  Q_corners : list
    List of 4 distinct 3D vectors (arrays / lists / tuples)
    being corners of a square (!) when projected on a XY-plane.
    They must be ordered in the following way: 
    first bottom-left and then clockwise.
  
  Returns
  -------
  H11 : float
    First element of the Hessian diagonal
  H22 : float 
    Second element of the Hessian diagonal
  
  Notes
  -----
  It is NOT written in a general form, using a loop over coordinates
  to remain explicit.
  
  """
  this_func = this_lib + 'Hessian_Bilinear_Num: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from numpy.linalg import norm
  from lib_math_interp import Interpolate_Bilinear
  from constants import differential
  
  x, y, z = m  
  
  # FINITE DIFFERENCES
  dx = differential # FIXME: THE SMALLER THE BETTER?
  dy = dx
  
  
  # X-COORDINATES OF CENTRAL-DIFFERENCE
  x_plus_dx = x + dx
  x_minus_dx = x - dx
  
  # POINTS ON SURFACE AT X-COORDINATES
  z_x_plus_dx = Interpolate_Bilinear(x_plus_dx, y, Q_corners)
  z_x_minus_dx = Interpolate_Bilinear(x_minus_dx, y, Q_corners)
  
  r_x_plus_dx = np.array([x_plus_dx, y, z_x_plus_dx])
  r_x_minus_dx = np.array([x_minus_dx, y, z_x_minus_dx])
  
  # X COMPONENT OF THE FINITE DIFFERENCE OF THE ERROR FUNCTION
  phi_x_plus_dx = norm(m_ref - r_x_plus_dx)
  phi_x_minus_dx = norm(m_ref - r_x_minus_dx)  
  phi_x = norm(m_ref - np.array([x, y, z]))
      
  # FD-APPROX. OF 2ND X-PARTIAL DERIVATIVE OF THE ERROR FUNCTION
  d2_phi_dx2 = (phi_x_plus_dx - 2 * phi_x + phi_x_minus_dx) / (dx ** 2)
  H11 = d2_phi_dx2

  
  # Y-COORDINATES OF CENTRAL-DIFFERENCE
  y_plus_dy = y + dy
  y_minus_dy = y - dy
  
  # POINTS ON SURFACE AT Y-COORDINATES
  z_y_plus_dy = Interpolate_Bilinear(x, y_plus_dy, Q_corners)
  z_y_minus_dy = Interpolate_Bilinear(x, y_minus_dy, Q_corners)
  
  r_y_plus_dy = np.array([x, y_plus_dy, z_y_plus_dy])
  r_y_minus_dy = np.array([x, y_minus_dy, z_y_minus_dy])
  
  # Y COMPONENT OF THE FINITE DIFFERENCE OF THE ERROR FUNCTION
  phi_y_plus_dy = norm(m_ref - r_y_plus_dy)
  phi_y_minus_dy = norm(m_ref - r_y_minus_dy)
  phi_y = norm(m_ref - np.array([x, y, z]))
      
  # FD-APPROX. OF 2ND Y-PARTIAL DERIVATIVE OF THE ERROR FUNCTION
  d2_phi_dy2 = (phi_y_plus_dy - 2 * phi_y + phi_y_minus_dy) / (dy ** 2)
  H22 = d2_phi_dy2
   
  #print this_func + 'END'
  return H11, H22  


#-----------------------------------------------------------------------------------------------------------


def Update_Model(problem, m, grad, H, args): 
  """
  Update the iteratively optimized model given the 
  gradient and the Hessian.
  
  Parameters
  ----------
  problem : str 
    Type of the problem to solve. It determines the method 
    of inverting the Hessian.
    Current capabilities:
      - 'lin2d' (minimum distance from a 2D surface interpolated bilinearly).
  m : array 
    Point in the model space for which to calcalute the gradient.
  grad : array
    Gradient of the error function.
  H : array
    Hessian of the error function.  
  args : list
    Additional parameters; problem-specific; passed on to inner functions.   
  
  Returns 
  -------
  m : array
    Updated model.
    
  Notes
  -----
  FIXME: mess with "alpha*" approximations.
  
  """
  this_func = this_lib + 'Update_Model: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  dm = Model_Difference(problem, grad, H)
  m0 = m # FOR STEP-LENGHT CALC
  
  if problem == 'lin1d':
    #FIXME
    x = m[0]
    a, b = args
    x += dm
    y = a * x + b
    m = [x, y]
  
  elif problem == 'lin2d':
    m_ref = args[0]
    Q_corners = args[1: ]
    m = Update_Model_Bilinear(dm, m, m_ref, Q_corners)
  
  else:
    eprint(this_func + 'Error. Problem: ' + problem + ' not yet implemented\n')
    quit()

  m = np.array(m)
  
  #print this_func + 'END'
  return m


#-----------------------------------------------------------------------------------------------------------


def Update_Model_Bilinear(dm, m, m_ref, Q_corners): 
  """
  Apply the precomputed update dm to the current model m
  so that it never goes beyond the bilinear tile.
  
  Parameters
  ----------
  dm : array
    Point in the model space expressing the model update.
  m : array 
    Point in the model space for which to calcalute the gradient.
  m_ref : array
    Point in the model space of the same dimension as m. It is a reference
    point from which the distance to the surface is measured.
  Q_corners : list
    List of 4 distinct 3D vectors (arrays / lists / tuples)
    being corners of a square (!) when projected on a XY-plane.
    They must be ordered in the following way: 
    first bottom-left and then clockwise.
  
  Returns 
  -------
  m : array
    Updated model.
    
  Notes
  -----
  It iteratively decreases the factor multiplying dm
  so that eventually the updated model falls within
  the bilinear tile, provided the current one fell
  within it too. This assumption is valid since 
  other option is not allowed throughout the optimization.
  
  """
  this_func = this_lib + 'Update_Model_Bilinear: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_math_interp import Interpolate_Bilinear, Unpack_Args_Bilinear_Interp, Model_Is_Inside_Bilinear_Tile
  

  x, y, z, dx, dy, dx_G, dy_G, dz_G, x1, x2, y1, y2, fQ11, fQ12, fQ22, fQ21 = Unpack_Args_Bilinear_Interp(m, m_ref, Q_corners)
  
  dm_x = dm[0]
  dm_y = dm[1]
  
  i = 1
  damp = 1.0
  while True:
    #print this_func, 'Iteration no. ', i
    x_test = x + damp * dm_x
    y_test = y + damp * dm_y
    z_test = Interpolate_Bilinear(x_test, y_test, Q_corners)
    m_test = [x_test, y_test, z_test]
  
    if Model_Is_Inside_Bilinear_Tile(m_test, x1, x2, y1, y2):
      break
    
    else:
      damp *= 0.5
    
    i += 1
    
  m = m_test
      
  #print this_func + 'END'
  return m
  
  
#-----------------------------------------------------------------------------------------------------------


def Model_Difference(problem, grad, H):
  """
  Calculate the Newton's optimization update
  using the gradient and the Hessian (derived 
  from the Taylor's expansion).
  
  Parameters
  ----------
  problem : str 
    Type of the problem to solve. It determines the method 
    of inverting the Hessian.
    Current capabilities:
      - 'lin2d' (minimum distance from a 2D surface interpolated bilinearly).
  grad : array
    Gradient of the error function.
  H : array
    Hessian of the error function.
  
  Returns
  -------
  dm : array
    Model update.
    
  Notes
  -----
  It is derived from the Taylor series expansion 
  up to quadratic terms inclusive.
  It is expressed as an antigradient preconditioned 
  by the inverse of the Hessian.
  FIXME: DOES IT NEED SCALING?
  See e.g. Mike Warner's notes on a basic theory of FWI.
  
  """
  this_func = this_lib + 'Model_Difference: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  H_inv = Invert_Hessian(problem, H)
  dm = H_inv.dot(-grad) # Newton's update
  
  #print this_func + 'END'
  return dm


#-----------------------------------------------------------------------------------------------------------


def Invert_Hessian(problem, H):
  """
  Calculate inverse of the Hessian matrix.
  
  Parameters
  ----------
  problem : str 
    Type of the problem to solve. It determines the method 
    of inverting the Hessian.
    Current capabilities:
      - 'lin2d' (minimum distance from a 2D surface interpolated bilinearly).
  H : array / float 
    Hessian matrix. For some problems it might be just a number.
 
  Returns
  -------
  H_inv : array 
  
  Notes
  -----
  FIXME; ADD DIMENSIONALITY CHECKS?
  
  """
  this_func = this_lib + 'Invert_Hessian: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_math import Invert_Matrix_2x2
  
  if problem == 'lin1d':
    inv_H = 1. / H 
  
  elif problem == 'lin2d':
    inv_H = Invert_Matrix_2x2(H)
  
  else:
    eprint(this_func + 'Error. Problem: ' + problem + ' not yet implemented\n')
    quit()    
  
  #print this_func, 'inv_H', inv_H
  #print this_func + 'END'  
  return inv_H


#-----------------------------------------------------------------------------------------------------------
















#
def Calculate_Step_Length(m0, m1, G):
  """
  According to "FWI basics" by MW
  
  """
  
  this_func = this_lib + 'Calculate_Step_Length: '
  print(this_func + 'START')

  d_d0 = dist3d(m0, G)
  d_d1 = dist3d(m1, G)
  
  q = abs(d_d0 - d_d1)
  
  if q > epsi:
    alpha = d_d0 / q
  else:
    eprint('Error. Value of q (see <FWI basics> by MW) close to zero \n')
    quit()

  print(this_func + 'END')
  return alpha


# FIXME: args CHANGED
def Gradient_Linear(m, m_ref, a, b):
  """
  FIXME: ADD Lin1d which used to work like:
  problem = 'lin1d'
  m0 = [1, -1]
  m = [13, 13]
  a = 1
  b = 0
  args = [a, b]
  niters = 3  
  
  """
  
  x = m[0]
  y = m[1]
  x_ref = m_ref[0]
  y_ref = m_ref[1]
  return x * (a**2 + 1) + a * (b - y_ref) - x_ref 


# FIXME: args CHANGED
def Hessian_Linear(a):
  return a**2 + 1


