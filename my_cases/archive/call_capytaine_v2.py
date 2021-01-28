# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 10:08:14 2020

@author: akeeste
"""

import numpy as np
# import xarray as xr
import capytaine as cpt
import meshmagick.mesh as mmm
import meshmagick.hydrostatics as mmhs
import xarray as xr
# import logging
# import os
# import sys

def __init__(self):
    LOG.info("Capytaine imported.")

def hydrostatics(myBody):
    '''
    use meshmagick functions to calculate and output the hydrostatic stiffness,
    interia, center of gravity, center of buoyancy and displaced volume of a 
    capytaine body

    Parameters
    ----------
    myBody : TYPE
        DESCRIPTION.

    Returns
    -------
    stiffness: array
        Hydrostatic stiffness of the submerged portion of the body
    
    center_of_mass: array [x,y,z]
        Coordinates of the body's center of gravity. If the body does not have 
        a defined cg, the function assumes constant density and this is 
        equivalent to the center of buoyancy (e.g. the center of volume)
    
    center_of_buoyancy: array [x,y,z]
        Coordinates of the body's center of buoyancy. This is the volumetric 
        center of the submerged part of the body.
    
    disp_vol: float
        Displaced volume of water by the submerged body.
    '''
    body_mesh = mmm.Mesh(myBody.mesh.vertices, myBody.mesh.faces, name= myBody.mesh.name)
    body_hs = mmhs.Hydrostatics(working_mesh=body_mesh,cog=myBody.center_of_mass,
                                           rho_water=1023.0,grav=9.81)
    
    vo = body_hs.displacement_volume
    cg = myBody.center_of_mass
    cb = body_hs.buoyancy_center
    k_hs = [[0, 0, 0,               0,               0,               0],
            [0, 0, 0,               0,               0,               0],
            [0, 0, body_hs.S33/1e4, body_hs.S34/1e4, body_hs.S35/1e4, 0],
            [0, 0, body_hs.S34/1e4, body_hs.S44/1e4, body_hs.S45/1e4, 0],
            [0, 0, body_hs.S35/1e4, body_hs.S45/1e4, body_hs.S55/1e4, 0],
            [0, 0, 0,               0,               0,               0]]
    # k_hs = k_hs/(1e-4)
    # k_hs = body_hs.hs_data['stiffness_matrix']
    
    myBody.hydrostatic_stiffness = myBody.add_dofs_labels_to_matrix(k_hs)
    myBody.displaced_volume = xr.DataArray(data=np.asarray(vo))
    myBody.center_of_mass = xr.DataArray(data=np.asarray(cg), dims=['xyz'],
                                coords={'xyz': ['x','y','z']},
                                )
    myBody.center_of_buoyancy = xr.DataArray(data=np.asarray(cb), dims=['xyz'],
                            coords={'xyz': ['x','y','z']},
                            )
    
    # return vo, cg, cb, k_hs
    return myBody



def call_capy_gbm(meshFName, wCapy, CoG=[0,0,0], headings=[0.0], saveNc=False,
              ncFName=None, wDes=None, body_name='', depth=-np.infty,
              density=1025.0):

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



def call_capy_nbod(meshFName, wCapy, CoG=[0,0,0], headings=[0.0], saveNc=False,
              ncFName=None, wDes=None, body_name='', depth=-np.infty,
              density=1025.0):

    # do not use: results +=; results.merge; results.update(); 
    
    # create capytaine body objects. If one input parameter related to the 
    #     bodies has len > 1, all must have the same length
    # if len(meshFName) > 1:
    body_dict = {}
    for i in np.arange(0, len(meshFName)):
        body_dict["body{0}".format(i)] = cpt.FloatingBody.from_file(meshFName[i])
        body_dict["body{0}".format(i)].center_of_mass = CoG[i]
        body_dict["body{0}".format(i)].keep_immersed_part()
        if body_name != '':
            body_dict["body{0}".format(i)].name = body_name[i]
        body_dict["body{0}".format(i)].add_all_rigid_body_dofs()
        body_dict["body{0}".format(i)] = hydrostatics(body_dict["body{0}".format(i)])
    bodies = list(body_dict.values())
    
    # for i in np.arange(0,len(cog)):
    #     body[i] = cpt.FloatingBody.from_file(meshFName[i])
    #     body[i].add_all_rigid_body_dofs()
    #     body[i].center_of_mass = CoG[i]
    #     body[i].keep_immersed_part()
    #     if body_name != '':
    #         body.name[i] = body_name[i]
    
    
    # define the hydrodynamic problems
    problems = [cpt.RadiationProblem(body=body1,
                                     radiating_dof=dof,
                                     omega=w,
                                     sea_bottom=depth,
                                     g=9.81,
                                     rho=density)
                                     for body1 in bodies for dof in body1.dofs for w in wCapy]

    problems += [cpt.DiffractionProblem(body=body1,
                                        omega=w,
                                        wave_direction=heading,
                                        sea_bottom=depth,
                                        g=9.81,
                                        rho=density)
                                        for body1 in bodies for w in wCapy for heading in headings]

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

    return capyData



def call_capy(meshFName, wCapy, CoG=[0,0,0], headings=[0.0], saveNc=False,
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
    
    # add hydrostatics data (cg, cb, vo, C) to the body
    body = hydrostatics(body)
    
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
    capyData = cpt.assemble_dataset(results,
                                    hydrostatics=True)
    
    
    # add kochin diffraction results
    # kochin = cpt.io.xarray.kochin_data_array(results, headings)
    # capyData.update(kochin)
    
    # use to test read_capytaine_v1() ability to catch and reorder dofs
    # sorted_idofs = ["Heave", "Sway", "Pitch", "Surge", "Yaw", "Roll"]
    # sorted_rdofs = ["Sway", "Heave", "Surge", "Roll", "Pitch", "Yaw"]
    # capyData = capyData.sel(radiating_dof=sorted_idofs, influenced_dof=sorted_rdofs)


    # (optional) save to .nc file
    if saveNc == True:
        cpt.io.xarray.separate_complex_values(capyData).to_netcdf(ncFName)
    
    return capyData
