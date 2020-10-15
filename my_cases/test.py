# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 16:07:23 2020

@author: akeeste
"""
import capytaine as cpt

# plate356 = cpt.FloatingBody.from_file('plate356.stl')
# plate772 = cpt.FloatingBody.from_file('plate772.stl')
# plate1649 = cpt.FloatingBody.from_file('plate1649.stl')

# plate356.translate_y(20)
# plate1649.translate_y(40)

# plates = plate356 + plate772 + plate1649
# plates.show()

plate = cpt.FloatingBody.from_file('plate.stl')

plate_ref = cpt.FloatingBody.from_file('plate_ref.stl')
plate_ref.translate_y(30)

plate_test = cpt.FloatingBody.from_file('plate_test.stl')
plate_test.translate_y(60)
# test.show()


buoy = cpt.FloatingBody.from_file('float.stl')
buoy.translate_x(30)

buoy_ref = cpt.FloatingBody.from_file('float_ref.stl')
buoy_ref.translate_y(30)
buoy_ref.translate_x(30)

buoy_test = cpt.FloatingBody.from_file('float_test.stl')
buoy_test.translate_y(60)
buoy_test.translate_x(30)


allbod = plate + plate_ref + plate_test + buoy + buoy_ref + buoy_test
allbod.show()