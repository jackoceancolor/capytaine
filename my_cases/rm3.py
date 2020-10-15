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

logging.basicConfig(level=logging.INFO)

# create buoy (float of the WEC, but 'float' is reserved in Python)
buoy_file = os.getcwd() + os.path.sep + 'float.stl'
buoy = cpt.FloatingBody.from_file(buoy_file)
buoy.add_translation_dof(name="Surge")
buoy.add_translation_dof(name="Sway")
buoy.add_translation_dof(name="Heave")
buoy.add_rotation_dof(name="Roll")
buoy.add_rotation_dof(name="Pitch")
buoy.add_rotation_dof(name="Yaw")
buoy.keep_immersed_part()

# create spar from input stl file
spar_file = os.getcwd() + os.path.sep + 'spar.stl'
spar = cpt.FloatingBody.from_file(spar_file)
spar.add_translation_dof(name="Surge")
spar.add_translation_dof(name="Sway")
spar.add_translation_dof(name="Heave")
spar.add_rotation_dof(name="Roll")
spar.add_rotation_dof(name="Pitch")
spar.add_rotation_dof(name="Yaw")
spar.keep_immersed_part()

all_dofs = ['Surge','Sway','Heave','Roll','Pitch','Yaw']


# Note: radiating_dof can only include up to the body's dofs 
test_matrix = xr.Dataset(coords={
    'omega': np.linspace(0.02, 5.2, 5),
    'wave_direction': np.linspace(0, np.pi/2, 2),
    'theta': np.linspace(0, np.pi/2, 2),
    'influenced_dof': all_dofs,
    # 'radiating_dof': list(body.dofs),  
    'radiating_dof': all_dofs,
    'water_depth': [np.infty],
})

solver = cpt.BEMSolver(green_function=cpt.XieDelhommeau(),
                       engine=cpt.BasicMatrixEngine())
# solver = cpt.BEMSolver()

dataset = solver.fill_dataset(
    test_matrix, 
    [body], 
    wavenumber=True, 
    wavelength=True,
    mesh=False, 
    hydrostatics=True,
    keep_details=True
    )


# save results in dataset to NetCDF file
dataset = cpt.io.xarray.separate_complex_values(dataset)
filename = os.getcwd() + os.path.sep + 'rm3.nc'
print(filename)
dataset.to_netcdf(filename)

