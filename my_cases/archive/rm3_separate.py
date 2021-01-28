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
from ..call_capytaine import call_capy # call_capytaine.py has some mods from david's original function
# from david_capytaine_read_function import call_capy


# FLOAT ----------------------------------------------------------------------#
# create buoy from input mesh (WAMIT or NEMOH) 
# buoy=float of the RM3, but 'float' is reserved in Python so 'buoy' used)

meshFName = os.getcwd() + os.path.sep + 'float.dat' # nemoh
# meshFName = os.getcwd() + os.path.sep + 'float'.gdf' # wamit

buoy_w = np.linspace(0.02, 8.4, 3) # 420 for full, 3 or 10 for tests
# buoy_w = (np.arange(260)+1)*0.02
buoy_cg = (0,0,-0.72)
buoy_headings = np.linspace(0,0,1)
buoy_nc = True
buoy_ncFile = os.getcwd() + os.path.sep + 'rm3_float_test.nc'

if os.path.isfile(buoy_ncFile):
    print(f'Output ({buoy_ncFile}) file already exists and will be overwritten. '
          'Do you wish to proceed? (y/n)')
    ans = input()
    if ans.lower() != 'y':
        print('\nEnding simulation. file not overwritten')
        sys.exit(0)


call_capy(meshFName=buoy_file,
                wCapy=buoy_w,
                CoG=buoy_cg,
                headings=buoy_headings,
                saveNc=buoy_nc,
                ncFName=buoy_ncFile,
                body_name='buoy_capytaine',
                depth = -np.infty,
                density = 1000.0)


print('\n\nFunction completed. Data is saved.\n')
# ----------------------------------------------------------------------------#



# FLOAT ----------------------------------------------------------------------#
# Create spar from input mesh (WAMIT or NEMOH)

meshFName = os.getcwd() + os.path.sep + 'spar.dat' # nemoh
# meshFName = os.getcwd() + os.path.sep + 'spar'.gdf' # wamit

spar_w = np.linspace(0.02, 8.4, 3) # 420 for full, 3 or 10 for tests
# spar_w = (np.arange(420)+1)*0.02
spar_cg = (0,0,-21.29)
spar_headings = np.linspace(0,0,1)
spar_nc = True
spar_ncFile = os.getcwd() + os.path.sep + 'rm3_spar_test.nc'

if os.path.isfile(spar_ncFile):
    print(f'Output ({spar_ncFile}) file already exists and will be overwritten. '
          'Do you wish to proceed? (y/n)')
    ans = input()
    if ans.lower() != 'y':
        print('\nEnding simulation. file not overwritten')
        sys.exit(0)


call_capy(meshFName=spar_file,
                wCapy=spar_w,
                CoG=spar_cg,
                headings=spar_headings,
                saveNc=spar_nc,
                ncFName=spar_ncFile,
                body_name='spar_capytaine',
                depth = -np.infty,
                density = 1000.0)


print('\n\nFunction completed. Data is saved.\n')
# ----------------------------------------------------------------------------#
