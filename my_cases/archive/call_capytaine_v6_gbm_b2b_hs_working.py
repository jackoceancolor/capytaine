# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 10:08:14 2020

@author: akeeste

Initial working script by David Ogden from:
https://github.com/mattEhall/FrequencyDomain/blob/b89dd4f4a732fbe4afde56efe2b52c3e32e22d53/FrequencyDomain.py#L842

"""

import numpy as np
# import xarray as xr
import capytaine as cpt
import meshmagick.mesh as mmm
import meshmagick.hydrostatics as mmhs
import xarray as xr
import logging as LOG
import os
# import sys

def __init__(self):
    LOG.info("Capytaine imported.")

"""
def hydrostatics(myBody):
    '''
    use meshmagick functions to calculate and output the hydrostatic stiffness,
    interia, center of gravity, center of buoyancy and displaced volume of a 
    capytaine body

    Parameters
    ----------
    myBody : FloatingBody
        A single capytaine floating body.

    Returns
    -------
    stiffness: array [6x6]
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
    
    cg = myBody.center_of_mass
    
    # meshmagick currently has issue if a body is copmletely submerged (OSWEC base)
    # use try-except statement to catch this error use alternate function for cb/vo
    # if completely submerged, stiffness is 0
    body_mesh = mmm.Mesh(myBody.mesh.vertices, myBody.mesh.faces, name= myBody.mesh.name)
    try:
        body_hs = mmhs.Hydrostatics(working_mesh=body_mesh,cog=myBody.center_of_mass,
                                            rho_water=1023.0,grav=9.81)
        vo = body_hs.displacement_volume
        cb = body_hs.buoyancy_center
        
        # # initialize stiffness matrix. Currently meshmagick on does 5 stiffness components (no gbm)
        # k_hs = np.zeros([6*nbod, 6*nbod])
        # tmp = [[body_hs.S33, body_hs.S34, body_hs.S35],
        #         [body_hs.S34, body_hs.S44, body_hs.S45],
        #         [body_hs.S35, body_hs.S45, body_hs.S55]]
        # k_hs[i*6+2:i*6+5,i*6+2:i*6+5] = tmp
        # myBody.hydrostatic_stiffness = myBody.add_dofs_labels_to_matrix(k_hs)
        
        hs_s = [body_hs.S33, body_hs.S34, body_hs.S35, body_hs.S44, body_hs.S45, body_hs.S55]
    except:
        # TODO: function to calculate the cb from mesh panels
        vo = body_mesh.volume
        cb = cg
        hs_s = [0, 0, 0, 0, 0, 0]
    
    
    hs_s = np.asarray(hs_s)
    myBody.hydrostatic_stiffness = xr.DataArray(data=hs_s, dims=['hydrostatic_S'],
                                coords={'hydrostatic_S': ['S33','S34','S35','S44','S45','S55']},
                                )
    myBody.displaced_volume = xr.DataArray(data=np.asarray(vo))
    myBody.center_of_mass = xr.DataArray(data=np.asarray(cg), dims=['xyz'],
                                coords={'xyz': ['x','y','z']},
                                )
    myBody.center_of_buoyancy = xr.DataArray(data=np.asarray(cb), dims=['xyz'],
                            coords={'xyz': ['x','y','z']},
                            )
        
    # return body with hydrostatics data
    return myBody
"""

def body_hydrostatics(myBodies, savepath=''):
    '''
    use meshmagick functions to calculate and output the hydrostatic stiffness,
    interia, center of gravity, center of buoyancy and displaced volume of a 
    capytaine bodies. Output is saved to Hydrostatics.dat and KH.dat files in 
    the same manner as Nemoh
    
    Example of output format:
    Hydrostatics.dat:
         XF =   0.000 - XG =   0.000
         YF =   0.000 - YG =   0.000
         ZF =  -2.500 - ZG =  -2.500
         Displacement =  0.4999997E+03
         Waterplane area =  0.1000002E+03
         
    KH.dat:
        0.0000000E+00  0.0000000E+00  0.0000000E+00  0.0000000E+00  0.0000000E+00  0.0000000E+00
        0.0000000E+00  0.0000000E+00  0.0000000E+00  0.0000000E+00  0.0000000E+00  0.0000000E+00
        0.0000000E+00  0.0000000E+00  0.9810053E+06 -0.1464844E-01 -0.5859375E-02  0.0000000E+00
        0.0000000E+00  0.0000000E+00 -0.1464844E-01  0.8160803E+07  0.0000000E+00  0.0000000E+00
        0.0000000E+00  0.0000000E+00 -0.5859375E-02  0.0000000E+00  0.8160810E+07  0.0000000E+00
        0.0000000E+00  0.0000000E+00  0.0000000E+00  0.0000000E+00  0.0000000E+00  0.0000000E+00


    Parameters
    ----------
    myBodies : List
        A list of capytaine floating bodies.

    Returns
    -------
    None
    '''
    nbod = len(myBodies)
    
    for i,body in enumerate(myBodies):
        cg = body.center_of_mass
        
        # use meshmagick to compute hydrostatic stiffness matrix
        # NOTE: meshmagick currently has issue if a body is copmletely submerged (OSWEC base)
        # use try-except statement to catch this error use alternate function for cb/vo
        # if completely submerged, stiffness is 0
        body_mesh = mmm.Mesh(body.mesh.vertices, body.mesh.faces, name= body.mesh.name)
        try:
            body_hs = mmhs.Hydrostatics(working_mesh=body_mesh,
                                        cog=body.center_of_mass,
                                        rho_water=1023.0,
                                        grav=9.81)
            vo = body_hs.displacement_volume
            cb = body_hs.buoyancy_center
            khs = body_hs.hydrostatic_stiffness_matrix
        except:
            # Exception if body is fully submerged
            vo = body_mesh.volume
            cb = cg
            hs_s = [0, 0, 0, 0, 0, 0]
            khs = np.zeros((3,3))
        
        # set file index
        fileind = ''
        if nbod != 1:
            fileind = '_' + str(i)
            
        # set file path based on where call_capytaine is called from (using ncFName)
        # filepath = ''
        # if savepath is not '':
        #     filepath,tmp = os.path.split(ncFName)
        
        # Write hydrostatic stiffness to KH.dat file as 
        khs_full = np.zeros((6,6))
        khs_full[2:5, 2:5] += khs
        tmp = savepath + 'KH' + fileind +'.dat'
        np.savetxt(tmp, khs_full)
        
        # Write the other hydrostatics data to Hydrostatics.dat file
        tmp = savepath + 'Hydrostatics' + fileind + '.dat'
        # f = open(rf'{tmp}','w')
        f = open(tmp,'w')
        for j in [0,1,2]:
            line =  f'XF = {cb[j]:7.3f} - XG = {cg[j]:7.3f} \n'
            f.write(line)
        line = f'Displacement = {vo:E}'
        f.write(line)
        f.close()
    
def call_capy_gbm(meshFName, wCapy, CoG=[0,0,0], headings=[0.0],
              ncFName=None, wDes=None, body_name='', depth=-np.infty,
              density=1025.0,additional_dofs_dir=None):
    '''
    Variation of the call_capy function to develop GBM functionality
    
    additional_dofs should be a dict of dicts (one dict of gbm dofs for each body)
    that defines generalized body modes for the body
    '''
    
    bodies = []
    for i in np.arange(0, len(meshFName)):
        bodies.append(cpt.FloatingBody.from_file(meshFName[i]))
        bodies[i].center_of_mass = CoG[i]
        bodies[i].keep_immersed_part()
        if body_name[i] != '':
            bodies[i].name = body_name[i]
        else:
            bodies[i].name = 'body'+str(i)
        # bodies[i].add_all_rigid_body_dofs()
    
    # output hydrostatics data to KH.dat and Hydrostatics.dat files
    path,tmp = os.path.split(ncFName)
    path += os.path.sep
    body_hydrostatics(bodies,path)
            
    
    # add gbm dofs
    # TODO: might need to call extra file to do this right. creating a custom gbm dof is easiest if the mesh is present
    # 1. pass flag or gbm_dofs.py file name to call_capy
    # 2. pass mesh to a local gbm_dofs.py script
    # 3. in the local gbm_dofs.py script, create dofs based on body mesh that is passed in
    # 4. pass dof dict back to call_capy and continue adding to body
    if additional_dofs_dir is not None:
        old_dir = os.getcwd()
        os.chdir(additional_dofs_dir)
        import gbm_dofs
        additional_dofs = gbm_dofs.new_dofs(body_name, bodies)
        
        for i in np.arange(0, len(meshFName)):
            if body_name[i] in additional_dofs:
                for k,v in additional_dofs[body_name[i]].items():
                    bodies[i].dofs[k] = v
        
        os.chdir(old_dir)
    
    for i in np.arange(0, len(meshFName)):
        bodies[i].add_all_rigid_body_dofs()
    
    
    
    # combine all bodies to account for B2B interactions
    combo = bodies[0]
    for i in np.arange(1,len(bodies),1):
        combo += bodies[i]
    
    # call Capytaine solver
    print(f'\n-------------------------------\n'
          f'Calling Capytaine BEM solver...\n'
          f'-------------------------------\n'
          f'mesh = {meshFName}\n'
          f'w range = {wCapy[0]:.3f} - {wCapy[-1]:.3f} rad/s\n'
          f'dw = {(wCapy[1]-wCapy[0]):.3f} rad/s\n'
          f'no of headings = {len(headings)}\n'
          # f'no of radiation & diffraction problems = {len(problems)}\n'
          f'-------------------------------\n')

    #    Create a dataset of parameters. 
    #    'fill_dataset()' automatically creates problems and solves them.
    problems = xr.Dataset(coords={
        'omega': wCapy,
        'wave_direction': headings,
        'radiating_dof': list(combo.dofs),
        'water_depth': [depth],
        'rho': [density],
        })
    
    solver = cpt.BEMSolver()
    capyData = solver.fill_dataset(problems, [combo], hydrostatics=False)
    
    # add kochin diffraction results
    # kochin = cpt.io.xarray.kochin_data_array(results, headings)
    # capyData.update(kochin)
    
    # use to test read_capytaine_v1() ability to catch and reorder dofs
    # sorted_idofs = ["Heave", "Sway", "Pitch", "Surge", "Yaw", "Roll"]
    # sorted_rdofs = ["Sway", "Heave", "Surge", "Roll", "Pitch", "Yaw"]
    # capyData = capyData.sel(radiating_dof=sorted_idofs, influenced_dof=sorted_rdofs)

    # save to .nc file
    cpt.io.xarray.separate_complex_values(capyData).to_netcdf(ncFName)
    
    print('\n\nCapytaine call complete. Data saved to \n' + ncFName)
    
    return capyData, problems


def call_capy(meshFName, wCapy, CoG=([0,0,0],), headings=[0.0],ncFName=None, 
              wDes=None, body_name=('',), depth=np.infty, density=1025.0):
    '''
    call Capytaine for a given mesh, frequency range and wave headings
    This function is modified from David Ogden's work 
    (see https://github.com/mattEhall/FrequencyDomain/blob/b89dd4f4a732fbe4afde56efe2b52c3e32e22d53/FrequencyDomain.py#L842 for the original function).
    
    May be called with multiple bodies (automatically implements B2B). In this case, 
    the meshFName, CoG, body_name should be a tuple of the values for each body
    
    Parameters
    ----------
    meshFName : tuple of strings
        Tuple containing a string for the path to each body's hydrodynamic mesh.
        mesh must be cropped at waterline (OXY plane) and have no lid
    wCapy: array
        array of frequency points to be computed by Capytaine
    CoG: tuple of lists
        tuple contains a 3x1 list of each body's CoG
    headings: list
        list of wave headings to compute
    saveNc: Bool
        save results to .nc file
    ncFName: str
        name of .nc file
    wDes: array
        array of desired frequency points
        (for interpolation of wCapy-based Capytaine data)
    body_name: tuple of strings
        Tuple containing strings. Strings are the names of each body. 
        Prevent the body name from being a long file path
    depth: float
        Water depth. Should be positive downwards. Use decimal value to prevent 
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
    - expand to generalized body modes using an additional_dof parameter
    '''
    
    # Create Capytaine floating bodies form each mesh file and calculate 
    # additional body properties (cg, dofs, hydrostatics).
    # body_dict = {}
    # for i in np.arange(0, len(meshFName)):
    #     body_dict["body{0}".format(i)] = cpt.FloatingBody.from_file(meshFName[i])
    #     body_dict["body{0}".format(i)].center_of_mass = CoG[i]
    #     body_dict["body{0}".format(i)].keep_immersed_part()
    #     if body_name != '':
    #         body_dict["body{0}".format(i)].name = body_name[i]
    #     body_dict["body{0}".format(i)].add_all_rigid_body_dofs()
    # bodies = list(body_dict.values())
    
    bodies = []
    for i in np.arange(0, len(meshFName)):
        bodies.append(cpt.FloatingBody.from_file(meshFName[i]))
        bodies[i].center_of_mass = CoG[i]
        bodies[i].keep_immersed_part()
        if body_name[i] != '':
            bodies[i].name = body_name[i]
        bodies[i].add_all_rigid_body_dofs()
    
    
    # output hydrostatics data to KH.dat and Hydrostatics.dat files
    path,tmp = os.path.split(ncFName)
    path += os.path.sep
    body_hydrostatics(bodies,path)
            
    
    # add gbm dofs here or before
    
    
    # combine all bodies to account for B2B interactions
    combo = bodies[0]
    for i in np.arange(1,len(bodies),1):
        combo += bodies[i]
    
    # call Capytaine solver
    print(f'\n-------------------------------\n'
          f'Calling Capytaine BEM solver...\n'
          f'-------------------------------\n'
          f'mesh = {meshFName}\n'
          f'w range = {wCapy[0]:.3f} - {wCapy[-1]:.3f} rad/s\n'
          f'dw = {(wCapy[1]-wCapy[0]):.3f} rad/s\n'
          f'no of headings = {len(headings)}\n'
          # f'no of radiation & diffraction problems = {len(problems)}\n'
          f'-------------------------------\n')

    #    Create a dataset of parameters. 
    #    'fill_dataset()' automatically creates problems and solves them.
    problems = xr.Dataset(coords={
        'omega': wCapy,
        'wave_direction': headings,
        'radiating_dof': list(combo.dofs),
        'water_depth': [depth],
        'rho': [density],
        })
    
    solver = cpt.BEMSolver()
    capyData = solver.fill_dataset(problems, [combo], hydrostatics=False)
    
    
    ###########################################################################
    # problems = [cpt.RadiationProblem(body=combo,
    #                                   radiating_dof=dof,
    #                                   omega=w,
    #                                   sea_bottom=depth,
    #                                   g=9.81,
    #                                   rho=density)
    #                                   for dof in combo.dofs for w in wCapy]
    # problems += [cpt.DiffractionProblem(body=combo,
    #                                     omega=w,
    #                                     wave_direction=heading,
    #                                     sea_bottom=depth,
    #                                     g=9.81,
    #                                     rho=density)
    #                                     for w in wCapy for heading in headings]
    
    # results = [solver.solve(problem, keep_details=True) for problem in sorted(problems)]
    # capyData = cpt.assemble_dataset(results,
    #                                 hydrostatics=True)
    ###########################################################################
    
    # add kochin diffraction results
    # kochin = cpt.io.xarray.kochin_data_array(results, headings)
    # capyData.update(kochin)
    
    # use to test read_capytaine_v1() ability to catch and reorder dofs
    # sorted_idofs = ["Heave", "Sway", "Pitch", "Surge", "Yaw", "Roll"]
    # sorted_rdofs = ["Sway", "Heave", "Surge", "Roll", "Pitch", "Yaw"]
    # capyData = capyData.sel(radiating_dof=sorted_idofs, influenced_dof=sorted_rdofs)

    # save to .nc file
    cpt.io.xarray.separate_complex_values(capyData).to_netcdf(ncFName)
    
    
    print('\n\nCapytaine call complete. Data saved to \n' + ncFName)
    
    return capyData, problems
