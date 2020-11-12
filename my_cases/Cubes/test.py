# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 13:20:30 2020

@author: akeeste
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 13:15:35 2020

@author: akeeste
This script recreates the RM3 model based on sample BEM 
parameters from WEC-Sim (frequency range, directions, etc)

"""

# setup environment
import numpy as np
import os
import sys
import capytaine as cpt
sys.path.insert(1,'C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases')

import call_capytaine as cc # call_capytaine.py has some mods from david's original function
# from david_capytaine_read_function import call_capy


# Load cube mesh files -------------------------------------------------------#
cube_files = (os.getcwd() + os.path.sep + 'r_cube.dat',
             os.getcwd() + os.path.sep + 't_cube.dat')# nemoh
# cube_files = (os.getcwd() + os.path.sep + 'r_cube.gdf',
#              os.getcwd() + os.path.sep + 't_cube.gdf')# wamit


# Set separate cube parameters -----------------------------------------------#
cube_cgs = ((0,0,-2.0),
           (100,0,-2.0))
body_names = ('r_cube_capytaine',
             't_cube_capytaine')


# Set shared cube parameters -------------------------------------------------#
cube_w = np.linspace(0.03, 15, 3) # 500 full, 3 for short test
cube_headings = np.arange(10.0)*10.0
cube_nc = True
cube_ncFile = os.getcwd() + os.path.sep + 'test.nc'


# Check for file overwriting -------------------------------------------------#
if os.path.isfile(cube_ncFile):
    print(f'Output ({cube_ncFile}) file already exists and will be overwritten. '
          'Do you wish to proceed? (y/n)')
    ans = input()
    if ans.lower() != 'y':
        print('\nEnding simulation. file not overwritten')
        sys.exit(0)


# Call capytaine -------------------------------------------------------------#
# cc.call_capy_nbod(meshFName=cube_files,
#                 wCapy=cube_w,
#                 CoG=cube_cgs,
#                 headings=cube_headings,
#                 saveNc=cube_nc,
#                 ncFName=cube_ncFile,
#                 body_name=body_names,
#                 depth = -20.0,
#                 density = 1000.0)

meshFName=cube_files
wCapy=cube_w
CoG=cube_cgs
headings=cube_headings
saveNc=cube_nc
ncFName=cube_ncFile
body_name=body_names
depth = -20.0
density = 1000.0

body_dict = {}
for i in np.arange(0, len(meshFName)):
    body_dict["body{0}".format(i)] = cpt.FloatingBody.from_file(meshFName[i])
    body_dict["body{0}".format(i)].center_of_mass = CoG[i]
    body_dict["body{0}".format(i)].keep_immersed_part()
    if body_name != '':
        body_dict["body{0}".format(i)].name = body_name[i]
    body_dict["body{0}".format(i)].add_all_rigid_body_dofs()
bodies = list(body_dict.values())

# for i in np.arange(0,len(cog)):
#     body[i] = cpt.FloatingBody.from_file(meshFName[i])
#     body[i].add_all_rigid_body_dofs()
#     body[i].center_of_mass = CoG[i]
#     body[i].keep_immersed_part()
#     if body_name != '':
#         body.name[i] = body_name[i]

# define the hydrodynamic problems
problems = [cpt.RadiationProblem(body=body1,
                                 radiating_dof=dof,
                                 omega=w,
                                 sea_bottom=depth,
                                 g=9.81,
                                 rho=density)
                                 for body1 in bodies for dof in body1.dofs for w in wCapy]

problems += [cpt.DiffractionProblem(body=body1,
                                    omega=w,
                                    wave_direction=heading,
                                    sea_bottom=depth,
                                    g=9.81,
                                    rho=density)
                                    for body1 in bodies for w in wCapy for heading in headings]



print('\n\nFunction completed. Data is saved.\n')
