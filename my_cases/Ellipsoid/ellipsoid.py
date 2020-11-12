# -*- coding: utf-8 -*-
"""
Created on Thu Nov  10 13:15:35 2020

@author: akeeste
This script recreates an ellipsoid model based on sample BEM 
parameters from WEC-Sim (frequency range, directions, etc)

"""

# setup environment
import numpy as np
import os
import sys
sys.path.insert(1,'C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases')

from call_capytaine import call_capy # call_capytaine.py has some mods from david's original function
# from david_capytaine_read_function import call_capy


# Load ellipsoid mesh file ------------------------------------------------------#
ellipsoid_file = os.getcwd() + os.path.sep + 'ellipsoid.dat' # nemoh mesh
# ellipsoid_file = os.getcwd() + os.path.sep + 'ellipsoid.gdf' # wamit mesh
ellipsoid_w = np.linspace(0.03, 9.24, 3) # 308 for full, 3 for tests
ellipsoid_cg = (0,0,0)
ellipsoid_headings = np.linspace(0,0,1)
ellipsoid_nc = True
ellipsoid_ncFile = os.getcwd() + os.path.sep + 'ellipsoid_test.nc'
# ----------------------------------------------------------------------------#

if os.path.isfile(ellipsoid_ncFile):
    print(f'Output ({ellipsoid_ncFile}) file already exists and will be overwritten. '
          'Do you wish to proceed? (y/n)')
    ans = input()
    if ans.lower() != 'y':
        print('\nEnding simulation. file not overwritten')
        sys.exit(0)


call_capy(meshFName=ellipsoid_file,
                wCapy=ellipsoid_w,
                CoG=ellipsoid_cg,
                headings=ellipsoid_headings,
                saveNc=ellipsoid_nc,
                ncFName=ellipsoid_ncFile,
                body_name='ellipsoid_capytaine',
                depth = -np.infty,
                density = 1000.0)


print('\n\nFunction completed. Data is saved.\n')
