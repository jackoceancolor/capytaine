# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 17:43:45 2020

@author: akeeste
This script recreates the RM3 model based on sample BEM 
parameters from WEC-Sim (frequency range, directions, etc)

"""

# setup environment
import numpy as np
import os
import sys
sys.path.insert(1,'C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases')

import call_capytaine as cc# call_capytaine.py has some mods from david's original function
# from david_capytaine_read_function import call_capy

import capytaine as cpt
import xarray as xr
import pandas as pd
import logging


logging.basicConfig(level=logging.WARNING)


# Load sphere mesh file ------------------------------------------------------#
sphere_file = os.getcwd() + os.path.sep + 'sphere.dat' # nemoh
# sphere_file = os.getcwd() + os.path.sep + 'sphere.gdf' # wamit
sphere_w = np.linspace(0.02, 8.4, 3) # 420 for full, 3 or 10 for tests
sphere_cg = (0,0,-2.0)
sphere_headings = np.linspace(0,0,1)
sphere_nc = True
sphere_ncFile = os.getcwd() + os.path.sep + 'sphere_hs_test.nc'
# ----------------------------------------------------------------------------#


cd = cc.call_capy_hs(meshFName=sphere_file,
                wCapy=sphere_w,
                CoG=sphere_cg,
                headings=sphere_headings,
                saveNc=sphere_nc,
                ncFName=sphere_ncFile,
                body_name='sphere_capytaine',
                depth = -50.0,
                density = 1000.0)
# ----------------------------------------------------------------------------#

# ----------------------------------------------------------------------------#
# # Start copy of call_capy_hs functionality for debugging----------------------#
meshFName=sphere_file
wCapy=sphere_w
CoG=sphere_cg
headings=sphere_headings
saveNc=sphere_nc
ncFName=sphere_ncFile
body_name='sphere_capytaine'
depth = -50.0
density = 1000.0


body = cpt.FloatingBody.from_file(meshFName)
body.center_of_mass = CoG
body.keep_immersed_part()
if body_name != '':
    body.name = body_name
body.add_all_rigid_body_dofs()

# add hydrostatics data from mesh (Vo, cg, cb, K_hs)
# (vo, cg, cb, k_hs) = cc.hydrostatics(body)
# body.hydrostatic_stiffness = body.add_dofs_labels_to_matrix(k_hs)
# body.displaced_volume = xr.DataArray(data=np.asarray(vo))
# body.center_of_mass = xr.DataArray(data=np.asarray(cg), dims=['xyz'],
#                             coords={'xyz': ['x','y','z']},
#                             )
# body.center_of_buoyancy = xr.DataArray(data=np.asarray(cb), dims=['xyz'],
#                             coords={'xyz': ['x','y','z']},
#                             )
body = cc.hydrostatics(body)

# return xr.DataArray(data=np.asarray(vector), dims=['influenced_dof'],
#                             coords={'influenced_dof': list(self.dofs)},
#                             )

# return xr.DataArray(data=np.asarray(matrix), dims=['influenced_dof', 'radiating_dof'],
#                     coords={'influenced_dof': list(self.dofs), 'radiating_dof': list(self.dofs)},
#                     )


# define the hydrodynamic problems
problems = [cpt.RadiationProblem(body=body,
                                  radiating_dof=dof,
                                  omega=w,
                                  sea_bottom=depth,
                                  g=9.81,
                                  rho=density)
                                  for dof in body.dofs for w in wCapy]

problems += [cpt.DiffractionProblem(body=body,
                                    omega=w,
                                    wave_direction=heading,
                                    sea_bottom=depth,
                                    g=9.81,
                                    rho=density)
                                    for w in wCapy for heading in headings]

solver = cpt.BEMSolver()
results = [solver.solve(problem, keep_details=True) for problem in sorted(problems)]
capyData = cpt.assemble_dataset(results,hydrostatics=True)
    

#-----------------------------------------------------------------------------#
# call for HS in cpt.io.xarray
from datetime import datetime
attrs=None
if attrs is None:
    attrs = {}
attrs['creation_of_dataset'] = datetime.now().isoformat()
    
bodies = list({result.body for result in results})
# hs_data = cpt.io.xarray.hydrostatics_dataset(bodies)
    
# replacing hydrostatics_dataset() function call
hs_data = xr.Dataset()
for body_property in ['mass', 'center_of_buoyancy', 'center_of_mass', 'displaced_volume' ,'hydrostatic_stiffness']:
    bodies_properties = {body.name: body.__getattribute__(body_property) for body in bodies if hasattr(body, body_property)}
    if len(bodies_properties) > 0:
        bodies_properties = xr.concat(bodies_properties.values(), pd.Index(bodies_properties.keys(), name='body_name'))
        # bodies_properties = _squeeze_dimensions(bodies_properties, dimensions=['body_name'])
        hs_data = xr.merge([hs_data, {body_property: bodies_properties}])


capyData1 = xr.merge([capyData, hs_data])
capyData1.attrs.update(attrs)

# # End copy of call_capy_hs functionality for debugging -----------------------#
#-----------------------------------------------------------------------------#

print('\n\nFunction completed.\n')
