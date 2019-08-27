#!/usr/bin/env python
# coding: utf-8
"""Solver for the BEM problem.

Example
-------

::

    problem = RadiationProblem(...)
    result = BEMSolver(green_functions=..., engine=...).solve(problem)

"""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging

import numpy as np

from datetime import datetime

from capytaine.green_functions import Delhommeau
from capytaine.bem.engines import BasicEngine, HierarchicalToeplitzMatrices
from capytaine.io.xarray import problems_from_dataset, assemble_dataset, kochin_data_array

LOG = logging.getLogger(__name__)

class BEMSolver:

    def __init__(self, green_function=Delhommeau(), engine=BasicEngine()):
        self.green_function = green_function
        self.engine = engine

        self.exportable_settings = {
            **self.green_function.exportable_settings,
            **self.engine.exportable_settings
        }

    def solve(self, problem, keep_details=True):
        """Solve the linear potential flow problem.

        Parameters
        ----------
        problem: LinearPotentialFlowProblem
            the problem to be solved
        keep_details: bool, optional
            if True, store the sources and the potential on the floating body in the output object
            (default: True)

        Returns
        -------
        LinearPotentialFlowResult
            an object storing the problem data and its results
        """
        LOG.info("Solve %s.", problem)

        if problem.wavelength < 8*problem.body.mesh.faces_radiuses.max():
            LOG.warning(f"Resolution of the mesh (8×max_radius={8*problem.body.mesh.faces_radiuses.max():.2e}) "
                        f"might be insufficient for this wavelength (wavelength={problem.wavelength:.2e})!")

        S, K = self.engine.build_matrices(
            problem.body.mesh, problem.body.mesh,
            problem.free_surface, problem.sea_bottom, problem.wavenumber,
            self.green_function
        )
        sources = self.engine.linear_solver(K, problem.boundary_condition)
        potential = S @ sources

        result = problem.make_results_container()
        if keep_details:
            result.sources = sources
            result.potential = potential

        for influenced_dof_name, influenced_dof_vectors in problem.influenced_dofs.items():
            # Scalar product on each face:
            influenced_dof_normal = np.sum(influenced_dof_vectors * problem.body.mesh.faces_normals, axis=1)
            # Sum over all faces:
            integrated_potential = - problem.rho * np.sum(potential * influenced_dof_normal * problem.body.mesh.faces_areas)
            # Store result:
            result.store_force(influenced_dof_name, integrated_potential)
            # Depending of the type of problem, the force will be kept as a complex-valued Froude-Krylov force
            # or stored as a couple of added mass and radiation damping coefficients.

        LOG.debug("Done!")

        return result

    def solve_all(self, problems, **kwargs):
        """Solve several problems.
        Optional keyword arguments are passed to `Nemoh.solve`.

        Parameters
        ----------
        problems: list of LinearPotentialFlowProblem
            several problems to be solved

        Returns
        -------
        list of LinearPotentialFlowResult
            the solved problems
        """
        return [self.solve(problem, **kwargs) for problem in sorted(problems)]

    def fill_dataset(self, dataset, bodies, **kwargs):
        """Solve a set of problems defined by the coordinates of an xarray dataset.

        Parameters
        ----------
        dataset : xarray Dataset
            dataset containing the problems parameters: frequency, radiating_dof, water_depth, ...
        bodies : list of FloatingBody
            the bodies involved in the problems

        Returns
        -------
        xarray Dataset
        """
        attrs = {'start_of_computation': datetime.now().isoformat(),
                 **self.exportable_settings}
        problems = problems_from_dataset(dataset, bodies)
        if 'theta' in dataset.coords:
            results = self.solve_all(problems, keep_details=True)
            kochin = kochin_data_array(results, dataset.coords['theta'])
            dataset = assemble_dataset(results, attrs=attrs, **kwargs)
            dataset['kochin'] = kochin
        else:
            results = self.solve_all(problems, keep_details=False)
            dataset = assemble_dataset(results, attrs=attrs, **kwargs)
        return dataset

    def get_potential_on_mesh(self, result, mesh, chunk_size=50):
        """Compute the potential on a mesh for the potential field of a previously solved problem.
        Since the interaction matrix does not need to be computed in full to compute the matrix-vector product,
        only a few lines are evaluated at a time to reduce the memory cost of the operation.

        Parameters
        ----------
        result : LinearPotentialFlowResult
            the return of Nemoh's solver
        mesh : Mesh or CollectionOfMeshes
            a mesh
        chunk_size: int, optional
            Number of lines to compute in the matrix.
            (legacy, should be passed as an engine setting instead).

        Returns
        -------
        array of shape (mesh.nb_faces,)
            potential on the faces of the mesh

        Raises
        ------
        Exception: if the :code:`Result` object given as input does not contain the source distribution.
        """
        LOG.info(f"Compute potential on {mesh.name} for {result}.")

        if result.sources is None:
            raise Exception(f"""The values of the sources of {result} cannot been found.
            They probably have not been stored by the solver because the option keep_details=True have not been set.
            Please re-run the resolution with this option.""")

        if chunk_size > mesh.nb_faces:
            S = self.engine.build_S_matrix(
                mesh,
                result.body.mesh,
                result.free_surface, result.sea_bottom, result.wavenumber,
                self.green_function
            )
            phi = S @ result.sources

        else:
            phi = np.empty((mesh.nb_faces,), dtype=np.complex128)
            for i in range(0, mesh.nb_faces, chunk_size):
                S = self.engine.build_S_matrix_for_reconstruction(
                    mesh.extract_faces(list(range(i, i+chunk_size))),
                    result.body.mesh,
                    result.free_surface, result.sea_bottom, result.wavenumber,
                    self.green_function
                )
                phi[i:i+chunk_size] = S @ result.sources

        LOG.debug(f"Done computing potential on {mesh.name} for {result}.")

        return phi

    def get_free_surface_elevation(self, result, free_surface, keep_details=False):
        """Compute the elevation of the free surface on a mesh for a previously solved problem.

        Parameters
        ----------
        result : LinearPotentialFlowResult
            the return of the solver
        free_surface : FreeSurface
            a meshed free surface
        keep_details : bool, optional
            if True, keep the free surface elevation in the LinearPotentialFlowResult (default:False)

        Returns
        -------
        array of shape (free_surface.nb_faces,)
            the free surface elevation on each faces of the meshed free surface

        Raises
        ------
        Exception: if the :code:`Result` object given as input does not contain the source distribution.
        """
        fs_elevation = 1j*result.omega/result.g * self.get_potential_on_mesh(result, free_surface.mesh)
        if keep_details:
            result.fs_elevation[free_surface] = fs_elevation
        return fs_elevation


# LEGACY INTERFACE

def _arguments(f):
    """Returns the name of the arguments of the function f"""
    return f.__code__.co_varnames[:f.__code__.co_argcount]

class Nemoh(BEMSolver):
    """Solver for the BEM problem based on Nemoh's Green function. Legacy API.
    Parameters are dispatched to the Delhommeau class and to the engine
    (BasicEngine or HierarchicalToeplitzMatrices).
    """

    def __init__(self, **params):
        green_function = Delhommeau(
           **{key: params[key] for key in params if key in _arguments(Delhommeau.__init__)}
	)
        if 'hierarchical_matrices' in params and params['hierarchical_matrices']:
            engine = HierarchicalToeplitzMatrices(
               **{key: params[key] for key in params if key in _arguments(HierarchicalToeplitzMatrices.__init__)}
            )
        else:
            engine = BasicEngine(
               **{key: params[key] for key in params if key in _arguments(BasicEngine.__init__)}
            )

        super().__init__(green_function, engine)

    def build_matrices(self, *args, **kwargs):
        """Legacy API."""
        return self.engine.build_matrices(*args, **kwargs)

