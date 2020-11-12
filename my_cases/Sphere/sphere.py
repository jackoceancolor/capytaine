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

import call_capytaine as cc # call_capytaine.py has some mods from david's original function
# from david_capytaine_read_function import call_capy


# Load sphere mesh file ------------------------------------------------------#
sphere_file = ((os.getcwd() + os.path.sep + 'sphere.dat'),) # nemoh
# sphere_file = os.getcwd() + os.path.sep + 'sphere.gdf' # wamit

sphere_cg = ((0,0,-2.0),)
sphere_name = ('sphere_capytaine',)

sphere_w = np.linspace(0.02, 8.4, 10) # 420 for full, 10 for tests
# sphere_w = (np.arange(3)+1)*0.02
sphere_headings = np.linspace(0,0,1)
sphere_nc = True
sphere_ncFile = os.getcwd() + os.path.sep + 'sphere_test.nc'
# ----------------------------------------------------------------------------#

if os.path.isfile(sphere_ncFile):
    print(f'Output ({sphere_ncFile}) file already exists and will be overwritten. '
          'Do you wish to proceed? (y/n)')
    ans = input()
    if ans.lower() != 'y':
        print('\nEnding simulation. file not overwritten')
        sys.exit(0)


cc.call_capy_nbod(meshFName=sphere_file,
                wCapy=sphere_w,
                CoG=sphere_cg,
                headings=sphere_headings,
                saveNc=sphere_nc,
                ncFName=sphere_ncFile,
                body_name=sphere_name,
                depth = -50.0,
                density = 1000.0)


print('\n\nFunction completed. Data is saved.\n')
