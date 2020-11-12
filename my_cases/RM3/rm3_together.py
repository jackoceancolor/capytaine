# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 13:35:00 2020

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
# from ..call_capytaine import call_capy # call_capytaine.py has some mods from david's original function
import capytaine as cpt
# from david_capytaine_read_function import call_capy

# general settings
ncFName = os.getcwd() + os.path.sep + 'rm3_combo_nemoh.nc'
wCapy = np.linspace(0.02, 8.4, 260) # 260 full, 3 or 10 for short test
headings = np.linspace(0,0,1)
saveNc = True
depth = -np.infty
density = 1000.0


# create buoy from input mesh (WAMIT or NEMOH) 
# meshFName = os.getcwd() + os.path.sep + 'float.dat' # nemoh
meshFName = os.getcwd() + os.path.sep + 'float.gdf' # wamit
buoy = cpt.FloatingBody.from_file(meshFName)
buoy.center_of_mass = (0,0,-0.72)
buoy.keep_immersed_part()
buoy.name = 'buoy_capytaine'
buoy.add_all_rigid_body_dofs()

problems = [cpt.RadiationProblem(body=buoy,
                                 radiating_dof=dof,
                                 omega=w,
                                 sea_bottom=depth,
                                 g=9.81,
                                 rho=density)
                                 for dof in buoy.dofs for w in wCapy]

problems += [cpt.DiffractionProblem(body=buoy,
                                    omega=w,
                                    wave_direction=heading,
                                    sea_bottom=depth,
                                    g=9.81,
                                    rho=density)
                                    for w in wCapy for heading in headings]


# Create spar from input mesh (WAMIT or NEMOH)
# meshFName = os.getcwd() + os.path.sep + 'spar.dat' # nemoh
meshFName = os.getcwd() + os.path.sep + 'spar.gdf' # wamit
spar = cpt.FloatingBody.from_file(meshFName)
spar.center_of_mass = (0,0,-21.29)
spar.keep_immersed_part()
spar.name = 'spar_capytaine'
spar.add_all_rigid_body_dofs()

problems += [cpt.RadiationProblem(body=spar,
                                 radiating_dof=dof,
                                 omega=w,
                                 sea_bottom=depth,
                                 g=9.81,
                                 rho=density)
                                 for dof in spar.dofs for w in wCapy]

problems += [cpt.DiffractionProblem(body=spar,
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
