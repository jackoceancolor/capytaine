# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 10:21:15 2020

@author: akeeste
"""

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
sys.path.insert(1,'C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases')

import call_capytaine as cc # call_capytaine.py has some mods from david's original function


# Load oswec mesh file ------------------------------------------------------#
oswec_file = (os.getcwd() + os.path.sep + 'flap.dat',
             os.getcwd() + os.path.sep + 'base_shift.stl') # base_cut.stl, base.dat nemoh, .gdf wamit
oswec_cg = ((0,0,-3.90),
            (0,0,-10.90)) # centers of gravity
oswec_name = ('oswec_flap',
              'oswec_base') # body names

oswec_w = np.linspace(0.04, 20.0, 3) # 500 for full, 3 for tests
oswec_headings = np.linspace(0,90,10)
oswec_depth = -10.90

oswec_nc = True
oswec_ncFile = os.getcwd() + os.path.sep + 'oswec_test.nc'
# ----------------------------------------------------------------------------#

# if os.path.isfile(oswec_ncFile):
#     print(f'Output ({oswec_ncFile}) file already exists and will be overwritten. '
#           'Do you wish to proceed? (y/n)')
#     ans = input()
#     if ans.lower() != 'y':
#         print('\nEnding simulation. file not overwritten')
#         sys.exit(0)

# cc.call_capy(meshFName = oswec_file,
#              wCapy     = oswec_w,
#              CoG       = oswec_cg,
#              headings  = oswec_headings,
#              saveNc    = oswec_nc,
#              ncFName   = oswec_ncFile,
#              body_name = oswec_name,
#              depth     = oswec_depth,
#              density   = 1000.0)

# call_capy
import capytaine as cpt
meshFName = oswec_file
wCapy     = oswec_w
CoG       = oswec_cg
headings  = oswec_headings
saveNc    = oswec_nc
ncFName   = oswec_ncFile
body_name = oswec_name
depth     = oswec_depth
density   = 1000.0
        
b0 = cpt.FloatingBody.from_file(meshFName[0])
b0.center_of_mass = CoG[0]
b0.keep_immersed_part(sea_bottom=depth)
b0.name = body_name[0]
b0.add_all_rigid_body_dofs()

b1 = cpt.FloatingBody.from_file(meshFName[1])
b1.mesh.heal_mesh()
b1.center_of_mass = CoG[1]
b1.keep_immersed_part(sea_bottom=depth)
b1.name = body_name[1]
b1.add_all_rigid_body_dofs()

b = b0+b1
# b0 = cc.hydrostatics(b0)
b1 = cc.hydrostatics(b1)

# problem with OSWEC base mesh (BOTH .dat and .gdf)
# NEMOH mesh file still contains object with origin at base, not correct water line
# for nemoh and capy, mesh origin should be located at water line, not cog (this is what wamit does)
# move base mesh down 10.9 + h/2 m



print('\n\nFunction completed. OSWEC data is saved.\n')
