# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 17:31:53 2020

@author: akeeste

capytaine test run

"""

# setup environment
import numpy as np
import xarray as xr
import capytaine as cpt
import logging
import os

logging.basicConfig(level=logging.INFO)

# create sphere
body = cpt.Sphere(radius=2.0, center=(0,0,-1), 
                  name="my sphere", nphi=15, 
                  ntheta=15)
body.add_translation_dof(name="Surge")
body.add_translation_dof(name="Sway")
body.add_translation_dof(name="Heave")
body.add_rotation_dof(name="Roll")
body.add_rotation_dof(name="Pitch")
body.add_rotation_dof(name="Yaw")
body.keep_immersed_part()

all_dofs = ['Surge','Sway','Heave','Roll','Pitch','Yaw']

test_matrix = xr.Dataset(coords={
    'omega': np.linspace(0.1, 4, 20),
    'wave_direction': np.linspace(0, np.pi/2, 20),
    'theta': np.linspace(0, np.pi/2, 20),
    'influenced_dof': all_dofs,
    'radiating_dof': list(body.dofs),  # this can only include up to the body's dofs anyway
    # 'radiating_dof': ['Surge','Heave','Sway'],
    'water_depth': [np.infty],
})

# solver = cpt.BEMSolver(green_function='Delhommeau')
solver = cpt.BEMSolver()

dataset = solver.fill_dataset(
    test_matrix, 
    [body], 
    wavenumber=True, 
    wavelength=True,
    mesh=False, 
    hydrostatics=True
    )


# save results in dataset to NetCDF file
dataset = cpt.io.xarray.separate_complex_values(dataset)
filename = os.getcwd() + os.path.sep + 'capytaine1.nc'
print(filename)
dataset.to_netcdf(filename)

