# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 17:31:53 2020

@author: akeeste

capytaine test run

"""

import numpy as np
import xarray as xr
import capytaine as cpt

body = cpt.Sphere(radius=2.0, center=(0,0,-1), 
                  name="my sphere", nphi=15, 
                  ntheta=15)
body.add_translation_dof(name="Surge")
body.add_translation_dof(name="Heave")
body.keep_immersed_part()

test_matrix = xr.Dataset(coords={
    'omega': np.linspace(0.1, 4, 5),
    'wave_direction': np.linspace(0, np.pi/2, 2),
    'radiating_dof': list(body.dofs),
    'water_depth': [np.infty],
})
dataset = cpt.BEMSolver().fill_dataset(
    test_matrix, [body], wavenumber=True, wavelength=True,
    mesh=False, hydrostatics=True)



