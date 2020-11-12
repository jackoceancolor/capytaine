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
sys.path.insert(1,'C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases')

from call_capytaine import call_capy_nbod # call_capytaine.py has some mods from david's original function
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
cube_ncFile = os.getcwd() + os.path.sep + 'test2.nc'


# Check for file overwriting -------------------------------------------------#
if os.path.isfile(cube_ncFile):
    print(f'Output ({cube_ncFile}) file already exists and will be overwritten. '
          'Do you wish to proceed? (y/n)')
    ans = input()
    if ans.lower() != 'y':
        print('\nEnding simulation. file not overwritten')
        sys.exit(0)


# Call capytaine -------------------------------------------------------------#
cd,probs = call_capy_nbod(meshFName=cube_files,
                wCapy=cube_w,
                CoG=cube_cgs,
                headings=cube_headings,
                saveNc=cube_nc,
                ncFName=cube_ncFile,
                body_name=body_names,
                depth = -20.0,
                density = 1000.0)

# meshFName=cube_files
# wCapy=cube_w
# CoG=cube_cgs
# headings=cube_headings
# saveNc=cube_nc
# ncFName=cube_ncFile
# body_name=body_names
# depth = -20.0
# density = 1000.0


print('\n\nFunction completed. Data is saved.\n')
