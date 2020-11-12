# -*- coding: utf-8 -*-
"""
Created on Mon Nov  10 13:35:00 2020

@author: akeeste
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 13:15:35 2020

@author: akeeste
This script recreates the OSWEC model based on sample BEM 
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
ncFName = os.getcwd() + os.path.sep + 'oswec_combo_nemoh_test.nc'
wCapy = np.linspace(0.04, 20, 3) # 500 full, 3 for short test
headings = np.arange(10.0)*10.0
saveNc = True
depth = -10.90
density = 1000.0


# create flap from input mesh (WAMIT or NEMOH) 
meshFName = os.getcwd() + os.path.sep + 'flap.dat' # nemoh mesh
# meshFName = os.getcwd() + os.path.sep + 'flap.gdf' # wamit mesh
flap = cpt.FloatingBody.from_file(meshFName)
flap.center_of_mass = (0,0,-3.90)
flap.keep_immersed_part()
flap.name = 'flap_capytaine'
flap.add_all_rigid_body_dofs()

problems = [cpt.RadiationProblem(body=flap,
                                 radiating_dof=dof,
                                 omega=w,
                                 sea_bottom=depth,
                                 g=9.81,
                                 rho=density)
                                 for dof in flap.dofs for w in wCapy]

problems += [cpt.DiffractionProblem(body=flap,
                                    omega=w,
                                    wave_direction=heading,
                                    sea_bottom=depth,
                                    g=9.81,
                                    rho=density)
                                    for w in wCapy for heading in headings]


# Create base from input mesh (WAMIT or NEMOH)
meshFName = os.getcwd() + os.path.sep + 'base.dat' # nemoh mesh
# meshFName = os.getcwd() + os.path.sep + 'base.gdf' # wamit mesh
base = cpt.FloatingBody.from_file(meshFName)
base.center_of_mass = (0,0,-10.90)
base.keep_immersed_part()
base.name = 'base_capytaine'
base.add_all_rigid_body_dofs()

problems += [cpt.RadiationProblem(body=base,
                                 radiating_dof=dof,
                                 omega=w,
                                 sea_bottom=depth,
                                 g=9.81,
                                 rho=density)
                                 for dof in base.dofs for w in wCapy]

problems += [cpt.DiffractionProblem(body=base,
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
