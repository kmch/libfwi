"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permision writing to k.chrapkiewicz17@imperial.ac.uk.

"""


## CONSTS
this_lib = 'lib_fwi_fs_CONST.py/'


margin = 4                # MARGIN OF ZERO GHOST DATA
pad = 1                   # PLOTT MARGINS
zpad_min = 50             
node_rad = 4              # RADIUS WITHIN WHICH CUT FS
psize = 70                # POINT SIZE
nwidth = 5                # WITH OF NORMAL LINE TO FS
nopac = 0.5               # OPACITY OF NORMAL LINE TO FS
ghosts_no = 1             # PER EACH GHOST
intersects_no = 1         # PER EACH GHOST
ficts_no = 4              # PER EACH GHOST !!!!!!!!!!!!!!!!!!
aux_no = 4

# FLAGS FOR FS POINTS
in_flag = 0
acc_flag = -1
ext_flag = -666

# FLAGS FOR GHOST NODES
normal_flag = 0  
fs_node_flag = 1 
margin_flag = 2  


## FINITE-DIFFERENCE METHOD
differential = 1e-6       # FINITE APPROXIMATION OF INFINITESIMAL dx; CAN'T BE < ~1e-7 OTHERWISE CLIPPED TO 0
min_stencil_size = 2      
hicks_radius = 3          # FREE SURFACE