# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 15:44:14 2020

@author: akeeste
"""

import numpy as np
# import xarray as xr
import capytaine as cpt
# import logging
# import os
# import sys

def __init__(self):
    LOG.info("Capytaine imported.")


def call_capy_gbm(meshFName, wCapy, CoG=[0,0,0], headings=[0.0], saveNc=False,
              ncFName=None, wDes=None, body_name='', depth=-np.infty,
              density=1025.0):
    '''
    call Capytaine for a given mesh, frequency range and wave headings
    This function is taken from David Ogden's work 
    (see david_capytaine_read_function.py for the full function).
    
    May be called with multiple bodies (no B2B interaction). In this case, 
    the meshFName, CoG, body_name should be a tuple of the values for each body

    Parameters
    ----------
    meshFName : str
        string containing path to hydrodynamic mesh.
        mesh must be cropped at waterline (OXY plane) and have no lid
    wCapy: array
        array of frequency points to be computed by Capytaine
    CoG: list
        3x1 vector of body's CoG
    headings: list
        list of wave headings to compute
    saveNc: Bool
        save results to .nc file
    ncFName: str
        name of .nc file
    wDes: array
        array of desired frequency points
        (for interpolation of wCapy-based Capytaine data)
    body_name: str
        String that is the name of the body. Prevents the body name from being 
        a long file path
    depth: float
        Water depth. Should be negative. Use decimal value to prevent 
        Capytaine outputting int32 types. Default is -np.infty
    density: float
        Water density. Use decimal value to prevent Capytaine outputting int32 
        types. Default 1025.0

    Returns
    -------
    hydrodynamic coefficients; as computed or interpolated

    Notes
    -----
    TODO:
    - expand to multibody problems
    '''

    # create capytaine body object
    body = cpt.FloatingBody.from_file(meshFName)
    body.center_of_mass = CoG
    body.keep_immersed_part()
    if body_name != '':
        body.name = body_name
    body.add_all_rigid_body_dofs()
    
    
    
    # for i in np.arange(0,len(cog)):
    #     body[i] = cpt.FloatingBody.from_file(meshFName[i])
    #     body[i].add_all_rigid_body_dofs()
    #     body[i].center_of_mass = CoG[i]
    #     body[i].keep_immersed_part()
    #     if body_name != '':
    #         body.name[i] = body_name[i]

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

    # call Capytaine solver
    print(f'\n-------------------------------\n'
          f'Calling Capytaine BEM solver...\n'
          f'-------------------------------\n'
          f'mesh = {meshFName}\n'
          f'w range = {wCapy[0]:.3f} - {wCapy[-1]:.3f} rad/s\n'
          f'dw = {(wCapy[1]-wCapy[0]):.3f} rad/s\n'
          f'no of headings = {len(headings)}\n'
          f'no of radiation & diffraction problems = {len(problems)}\n'
          f'-------------------------------\n')

    solver = cpt.BEMSolver()
    results = [solver.solve(problem, keep_details=True) for problem in sorted(problems)]
    capyData = cpt.assemble_dataset(results)

    # (optional) save to .nc file
    if saveNc == True:
        cpt.io.xarray.separate_complex_values(capyData).to_netcdf(ncFName)


    