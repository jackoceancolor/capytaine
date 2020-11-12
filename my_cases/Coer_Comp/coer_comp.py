# -*- coding: utf-8 -*-
"""
Created on Thu Nov  10 13:15:35 2020

@author: akeeste
This script recreates an coercomp model based on sample BEM 
parameters from WEC-Sim (frequency range, directions, etc)

"""

# setup environment
import numpy as np
import os
import sys
sys.path.insert(1,'C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases')

from call_capytaine import call_capy # call_capytaine.py has some mods from david's original function
# from david_capytaine_read_function import call_capy


# Load coercomp mesh file ------------------------------------------------------#
coercomp_file = os.getcwd() + os.path.sep + 'coer_comp.dat' # nemoh mesh
# coercomp_file = os.getcwd() + os.path.sep + 'coer_comp.gdf' # wamit mesh
coercomp_w = np.linspace(0.1, 80.0, 3) # 800 for full, 3 for tests
coercomp_cg = (0,0,-0.2)
coercomp_headings = np.arange(10.0)*10.0
coercomp_nc = True
coercomp_ncFile = os.getcwd() + os.path.sep + 'coer_comp_test.nc'
# ----------------------------------------------------------------------------#

if os.path.isfile(coercomp_ncFile):
    print(f'Output ({coercomp_ncFile}) file already exists and will be overwritten. '
          'Do you wish to proceed? (y/n)')
    ans = input()
    if ans.lower() != 'y':
        print('\nEnding simulation. file not overwritten')
        sys.exit(0)


call_capy(meshFName=coercomp_file,
                wCapy=coercomp_w,
                CoG=coercomp_cg,
                headings=coercomp_headings,
                saveNc=coercomp_nc,
                ncFName=coercomp_ncFile,
                body_name='coercomp_capytaine',
                depth = -2.20,
                density = 1000.0)


print('\n\nFunction completed. Data is saved.\n')
