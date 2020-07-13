"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permision writing to k.chrapkiewicz17@imperial.ac.uk.

This library provides framework procedures for...
At the moment fully implemented is:

1. ...

"""

## MODULES
import matplotlib


## CONSTS
this_lib = 'lib_generic_CONST.py/'


# -------------------------------------------------------------------------------


## FILE NAMING
tmp_prefix = 'tmp_'


## CODE DEV.
verbos_func = -1 # PRINT this_func, ': START' AND NOTHING ELSE


## USEFUL CONSTANTS
epsi = 1e-10 # epsilon (small number)
big = 1e6
very_big = 1e15  
err_value = -666.666



# FIXME: MOVE TO PLT STYLE --------------------------------------------------------


## GLOBAL PLOTTING OPTIONS
fig_resol = 2000 # FIGURE RESOLUTION [dpi]
figsize_default = [24,10] # FOR JUPYTER
alpha_default = 0.5
slice_coord_default = 'y'
# ORDER OF APPEARING ON CANVAS (FRONT/BACK)
zorder_back = 0
zorder_mid = 10
zorder_front = 100

axes_lw = 2
tick_length = 5
tick_width = axes_lw
shift = 0.5 # FOR 2D MODELS: CENTERING PIXELS: X - shift, Z - shift)

font = {'size' : 15} #40
matplotlib.rc('font', **font)
matplotlib.rc('axes', lw = axes_lw)
matplotlib.rc('grid', lw = 2)
square_side = 5. # inches 
longer_side = 2 * square_side #20 
shorter_side = square_side / 2. 
tiny_side = 5 ##
opacity = 0.4 # the less the more transparent
lwidth = 2 #4
legend_font_size = 12
marker_size = 50
