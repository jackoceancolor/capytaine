# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 10:56:50 2020

@author: akeeste
"""

# setup environment
import numpy as np
import os
import sys
sys.path.insert(1,'C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases')

import call_capytaine as cc # call_capytaine.py has some mods from david's original function

import capytaine as cpt
import xarray as xr
import pandas as pd
import logging

logging.basicConfig(level=logging.WARNING)


# # Sphere
# files = (os.getcwd() + os.path.sep + 'Sphere/sphere.dat', )
# cgs = ((0,0,-2.0), )
# names = ('sphere', )

# # Ellipsoid
# files = (os.getcwd() + os.path.sep + 'Ellipsoid/ellipsoid.dat', )
# cgs = ((0,0,0), )
# names = ('ellipsoid', )

# # Cylinder
# files = (os.getcwd() + os.path.sep + 'Cylinder/cylinder.dat', )
# cgs = ((0,0,0), )
# names = ('sphere1', )

# # Cubes
# issues with t_cube
files = (os.getcwd() + os.path.sep + 'Cubes/r_cube.gdf',
              os.getcwd() + os.path.sep + 'Cubes/t_cube.gdf')
cgs = ((0,0,-2.5),
            (100,0,-2.5))
names = ('r_cube_capytaine',
              't_cube_capytaine')

# # Coer-comp
# files = (os.getcwd() + os.path.sep + 'Coer_Comp/coer_comp.dat', )
# cgs = ((0,0,-0.2), )
# names = ('coercomp', )

# # RM3
# files = (os.getcwd() + os.path.sep + 'RM3/float.dat',
#               os.getcwd() + os.path.sep + 'RM3/spar.dat')
# cgs = ((0,0,-0.72),
#         (0,0,-21.29))
# names = ('rm3_float',
#               'rm3_spar')

# # OSWEC
# issue with the mesh?
# files = (os.getcwd() + os.path.sep + 'OSWEC/flap.gdf',
#               os.getcwd() + os.path.sep + 'OSWEC/base.gdf')
# cgs = ((0,0,-3.90),
#         (0,0,-10.90))
# names = ('oswec_flap',
#               'oswec_base')


# Create body in same manner as call_capytaine -------------------------------#
meshFName = files
CoG = cgs
body_name = names

body_dict = {}
for i in np.arange(0, len(meshFName)):
    body_dict["body{0}".format(i)] = cpt.FloatingBody.from_file(meshFName[i])
    body_dict["body{0}".format(i)].center_of_mass = CoG[i]
    body_dict["body{0}".format(i)].keep_immersed_part()
    if body_name != '':
        body_dict["body{0}".format(i)].name = body_name[i]
    body_dict["body{0}".format(i)].add_all_rigid_body_dofs()
bodies = list(body_dict.values())



# Add hydrostatics to each mesh ----------------------------------------------#
for body in bodies:
    body = cc.hydrostatics(body)
    khs = body.hydrostatic_stiffness.values/(1023*9.81)
    tmp = np.zeros([3,3])
    tmp[0,0] = khs[0]
    tmp[0,1] = khs[1]
    tmp[0,2] = khs[2]
    tmp[1,0] = khs[1]
    tmp[1,1] = khs[3]
    tmp[1,2] = khs[4]
    tmp[2,0] = khs[2]
    tmp[2,1] = khs[4]
    tmp[2,2] = khs[5]
    print(f'body name:\n {body.name} \n\n'
          f'disp. vol.: \n {body.displaced_volume.values}\n\n'
          f'cog: \n {body.center_of_mass.values}\n\n'
          f'cob: \n {body.center_of_buoyancy.values}\n\n'
          # f'K_hs: \n {body.hydrostatic_stiffness.values[:]/(1023*9.81)}\n\n'
          f'K_hs: \n {tmp} \n\n'
          )


    




print('\n\nFunction completed.\n')
