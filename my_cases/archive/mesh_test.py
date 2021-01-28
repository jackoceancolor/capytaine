# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 16:20:29 2020

@author: akeeste
"""


import numpy as np
import os
import sys
sys.path.insert(1,'C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases')

import capytaine as cpt
import call_capytaine as cc # call_capytaine.py has some mods from david's original function


# Load oswec mesh file ------------------------------------------------------#
f1 = os.getcwd() + os.path.sep + 'OSWEC/flap.gdf'
cg1 = (0,0,-3.90)
name1 = 'oswec_flap'
depth1 = -10.90
body1 = cpt.FloatingBody.from_file(f1)
body1.center_of_mass = cg1
# body1.keep_immersed_part(sea_bottom=depth1)

f2 = os.getcwd() + os.path.sep + 'OSWEC/base.gdf'
cg2 = (0,0,-10.90)
name2 = 'oswec_base'
depth2 = -10.90
body2 = cpt.FloatingBody.from_file(f2)
body2.center_of_mass = cg2
# body2.keep_immersed_part(sea_bottom=depth2)

f3 = os.getcwd() + os.path.sep + 'OSWEC/flap.dat'
cg3 = (0,0,-3.90)
name3 = 'oswec_base'
depth3 = -10.90
body3 = cpt.FloatingBody.from_file(f3)
body3.center_of_mass = cg3
# body3.keep_immersed_part(sea_bottom=depth3)
body3.translate_x(10)

f4 = os.getcwd() + os.path.sep + 'OSWEC/base.dat'
cg4 = (0,0,-10.90)
name4 = 'oswec_base'
depth4 = -10.90
body4 = cpt.FloatingBody.from_file(f4)
body4.center_of_mass = cg4
# body4.keep_immersed_part(sea_bottom=depth4)
body4.translate_x(15)

bodies = body1+body2+body3+body4
bodies.show()





