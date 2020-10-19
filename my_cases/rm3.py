# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 13:41:01 2020

@author: akeeste
This script recreates the RM3 model based on sample BEM 
parameters from WEC-Sim (frequency range, directions, etc)

"""

# setup environment
import numpy as np
import xarray as xr
import capytaine as cpt
import logging
import os

logging.basicConfig(level=logging.WARNING)

# create buoy from input 
# (float of the WEC, but 'float' is reserved in Python)
# buoy_file = os.getcwd() + os.path.sep + 'float_ref.stl'

buoy_file = os.getcwd() + os.path.sep + 'Meshes\\nemoh_float.dat'
form = 'nemoh'

# buoy_file = os.getcwd() + os.path.sep + 'Meshes\\wamit_float.gdp'
# form = 'wamit'

buoy = cpt.FloatingBody.from_file(
    filename=buoy_file,
    file_format='stl',
    name='buoy')
buoy.translate_z(-0.72)
buoy.keep_immersed_part()
buoy.add_translation_dof(name="Surge")
buoy.add_translation_dof(name="Sway")
buoy.add_translation_dof(name="Heave")
buoy.add_rotation_dof(name="Roll")
buoy.add_rotation_dof(name="Pitch")
buoy.add_rotation_dof(name="Yaw")


# create plate from input file
# plate_file = os.getcwd() + os.path.sep + 'plate_ref.stl'
# form = 'stl'

plate_file = os.getcwd() + os.path.sep + 'Meshes\\nemoh_spar.dat'
form = 'nemoh'

# plate_file = os.getcwd() + os.path.sep + 'Meshes\\wamit_spar.gdp'
# form = 'wamit'

plate = cpt.FloatingBody.from_file(
    filename=plate_file,
    file_format=form,
    name='plate')
plate.translate_z(-21.29)
plate.keep_immersed_part()
plate.add_translation_dof(name="Surge")
plate.add_translation_dof(name="Sway")
plate.add_translation_dof(name="Heave")
plate.add_rotation_dof(name="Roll")
plate.add_rotation_dof(name="Pitch")
plate.add_rotation_dof(name="Yaw")

combo = plate+buoy
# combo.show()
###################################################################

all_dofs = ['Surge','Sway','Heave','Roll','Pitch','Yaw']


# Note: radiating_dof can only include up to the body's dofs 
test_matrix = xr.Dataset(coords={
    'omega': np.linspace(0.02, 5.2, 10),
    # 'omega': np.linspace(0.02, 5.2, 5),
    'wave_direction': np.linspace(0, np.pi/2, 2),
    'theta': np.linspace(0, np.pi/2, 2),
    
    'influenced_dof': all_dofs,
    # 'influenced_dof': list(combo.dofs),
    # 'influenced_dof': combo.dofs,
     
    # 'radiating_dof': all_dofs,
    # 'radiating_dof': list(combo.dofs),
    'radiating_dof': combo.dofs,
    
    'water_depth': [np.infty],
})

# solver = cpt.BEMSolver(green_function=cpt.XieDelhommeau(),
#                        engine=cpt.BasicMatrixEngine())
solver = cpt.BEMSolver()

rm3_results = solver.fill_dataset(
    test_matrix, 
    # [buoy, plate],
    [comob],
    wavenumber=True, 
    wavelength=True,
    mesh=True, 
    hydrostatics=True,
    # keep_details=True
    )


# save results in dataset to NetCDF file
rm3_results = cpt.io.xarray.separate_complex_values(rm3_results)
filename = os.getcwd() + os.path.sep + 'rm3.nc'
print(filename)
rm3_results.to_netcdf(filename)
