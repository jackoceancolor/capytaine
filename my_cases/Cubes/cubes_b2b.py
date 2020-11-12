# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 14:35:26 2020

@author: akeeste
This script recreates the cube model t_cubed on sample BEM 
parameters from WEC-Sim (frequency range, directions, etc)

"""

# setup environment
import numpy as np
import os
import sys
# from ..call_capytaine import call_capy # call_capytaine.py has some mods from david's original function
import capytaine as cpt
# from david_capytaine_read_function import call_capy

# general settings
ncFName = os.getcwd() + os.path.sep + 'cube_combo_nemoh.nc'
wCapy = np.linspace(0.03, 15, 3) # 500 full, 3 for short test
headings = np.arange(10.0)*10.0
saveNc = True
depth = -20.0
density = 1000.0


# create r_cube from input mesh (WAMIT or NEMOH) 
meshFName = os.getcwd() + os.path.sep + 'r_cube.dat' # nemoh mesh
# meshFName = os.getcwd() + os.path.sep + 'r_cube.gdf' # wamit mesh
r_cube = cpt.FloatingBody.from_file(meshFName)
r_cube.center_of_mass = (0,0,-2.5)
r_cube.keep_immersed_part()
r_cube.name = 'r_cube_capytaine'
r_cube.add_all_rigid_body_dofs()

problems = [cpt.RadiationProblem(body=r_cube,
                                 radiating_dof=dof,
                                 omega=w,
                                 sea_bottom=depth,
                                 g=9.81,
                                 rho=density)
                                 for dof in r_cube.dofs for w in wCapy]

problems += [cpt.DiffractionProblem(body=r_cube,
                                    omega=w,
                                    wave_direction=heading,
                                    sea_bottom=depth,
                                    g=9.81,
                                    rho=density)
                                    for w in wCapy for heading in headings]


# Create t_cube from input mesh (WAMIT or NEMOH)
meshFName = os.getcwd() + os.path.sep + 't_cube.dat' # nemoh mesh
# meshFName = os.getcwd() + os.path.sep + 't_cube.gdf' # wamit mesh
t_cube = cpt.FloatingBody.from_file(meshFName)
t_cube.center_of_mass = (100,0,-2.5)
t_cube.keep_immersed_part()
t_cube.name = 't_cube_capytaine'
t_cube.add_all_rigid_body_dofs()

problems += [cpt.RadiationProblem(body=t_cube,
                                 radiating_dof=dof,
                                 omega=w,
                                 sea_bottom=depth,
                                 g=9.81,
                                 rho=density)
                                 for dof in t_cube.dofs for w in wCapy]

problems += [cpt.DiffractionProblem(body=t_cube,
                                    omega=w,
                                    wave_direction=heading,
                                    sea_bottom=depth,
                                    g=9.81,
                                    rho=density)
                                    for w in wCapy for heading in headings]


if os.path.isfile(ncFName):
    print(f'Output ({ncFName}) file already exists and will be overwritten. '
          'Do you wish to proceed? (y/n)')
    ans = input()
    if ans.lower() != 'y':
        print('\nEnding simulation. file not overwritten')
        sys.exit(0)


print('Solving problems.\n')

solver = cpt.BEMSolver()
results = [solver.solve(problem) for problem in sorted(problems)]
capyData = cpt.assemble_dataset(results)

if saveNc == True:
    cpt.io.xarray.separate_complex_values(capyData).to_netcdf(ncFName)

print('\n\nFunction completed. Data is saved.\n')
# ----------------------------------------------------------------------------#
