"""
Components for torsiondrives of molecules for FF benchmarking.

These components leverage the QCArchive software stack for execution.

"""

import io
import os
import shutil
import contextlib

from ..geometry_optimizations.seasons import SEASONS
from ..exceptions import QCEngineError


class TorsiondriveExecutor:
    def __init__(self):
        self.stop = False

    def execute_torsiondrive_single(self,
            offmol, dihedrals, grid_spacing, dihedral_ranges,
            program, method, basis, season, 
            ncores=1, memory=2,
            scf_maxiter=200, geometric_maxiter=300, geometric_coordsys='dlc',
            geometric_qccnv=False):
        """Torsiondrive executor for a single molecule.

        Parameters
        ----------
        offmol : openff.toolkit.topology.Molecule
            Molecule with 

        dihedrals : List of tuples
            A list of the dihedrals to scan over.
        grid_spacing : List of int
            The grid seperation for each dihedral angle.
        dihedral_ranges : List of [low, high] pairs
            Pairs of ranges, with each range corresponding to each
            dihedral in `dihedrals`, e.g. `[[-120, 120], [-90, 150]]`.
            The dihedrals in given `offmol` conformation must be within these
            ranges.
        program : str
            Program to use for calculation.
            Either `season` or (`program`, `method`, `basis`) must be given, not both.
        method : str
            Method to use within program.
            Either `season` or (`program`, `method`, `basis`) must be given, not both.
        basis : str
            Basis to use for method.
            Either `season` or (`program`, `method`, `basis`) must be given, not both.
        season : str
            Season identifier specifying compute selections applied to molecules.
            See `..geometry_optimizations.seasons.SEASONS` for available keys.
            Either `season` or (`program`, `method`, `basis`) must be given, not both.
        ncores : int
            Number of threads to use for each gradient calculation step by the underlying `program`.
        memory : int
            Maximum memory in GiB to allocate for each gradient calculation.
            Set at least 5% *below* desired maximum.
        scf_maxiter : int
            Maximum iterations to allow for SCF convergence, if applicable.
        geometric_maxiter : int
            Maximum iterations to allow for geometry optimization convergence.
        geometric_coordsys : str
            Internal coordinate scheme to use for geometry optimization.
        geometric_qccnv : bool
            If set, use QChem-style convergence criteria for geometry optimization.

        Returns
        -------
        opts : Dict[str, Dict[str, OptimizationResult]]
            Optimization results, grouped by compute spec name, with dihedral angle as keys.

        """
        import numpy as np
        from qcelemental.models import Molecule
        from torsiondrive import td_api
        from geometric.nifty import bohr2ang, ang2bohr
        from openff.qcsubmit.factories import TorsiondriveDatasetFactory

        local_options={"ncores": ncores,
                       "memory": memory}

        factory = TorsiondriveDatasetFactory()

        # This block will fail for OFF Toolkit 0.8.4, but succeed for 0.8.4rc1
        try:
            attributes = factory.create_cmiles_metadata(offmol)
            qcmol = offmol.to_qcschema(extras=attributes)
        # In OFFTK 0.8.4, the CMILES field is automatically populated by this method
        except:
            qcmol = offmol.to_qcschema()

        # this is the object we care about returning
        opts = dict()

        if (season is not None) and any((program, method, basis)):
            raise ValueError("Cannot specify both `season` and any of (`program`, `method`, `basis`)")
        elif season is not None:
            specs = SEASONS[season].items()
        elif all((program, method)):
            specs = {'single': {
                        "program": program,
                        "method": method,
                        "basis": basis,
                        "spec_name": 'single',
                        "spec_description": 'user-specified compute specification'
                        },
                    }
        else:
            raise ValueError("Must specify at least (`program`, `method`) or `season`")

        for spec_name, compute_spec in specs.items():
            print("Processing spec: '{}'".format(spec_name))
            opts[spec_name] = dict()

            # generate initial torsiondrive state
            td_stdout = io.StringIO()
            with contextlib.redirect_stdout(td_stdout):
                
                # this object gets updated through the whole execution,
                # but it is an implementation detail; we do not return it
                state = td_api.create_initial_state(
                        dihedrals=dihedrals,
                        grid_spacing=grid_spacing,
                        elements=qcmol.symbols,
                        init_coords=[qcmol.geometry.flatten().tolist()],
                        dihedral_ranges=dihedral_ranges,
                        )

            # generate new grid points and optimize them until exhausted
            while True:
                next_jobs = td_api.next_jobs_from_state(state, verbose=True)

                # if no next jobs, we are finished
                if len(next_jobs) == 0:
                    break

                task_results = dict()
                for gridpoint in next_jobs:
                    opts[spec_name][gridpoint] = list()
                    task_results[gridpoint] = list()

                    # each gridpoint may have more than one job
                    for job in next_jobs[gridpoint]:
                        print(f"... '{gridpoint}'")

                        # set geometry of molecule to that of next job
                        qcmol_s = qcmol.copy(deep=True).dict()
                        qcmol_s['geometry'] = np.array(job).reshape(len(qcmol_s['symbols']), 3)
                        qcmol_s = Molecule.from_data(qcmol_s)

                        # generate constraints input
                        angles = gridpoint.split()
                        constraints = {"set":
                                [{"type": "dihedral", "indices": dihedral, "value": int(angle)}
                                    for dihedral, angle in zip(dihedrals, angles)]}


                        # perform optimization
                        input_data = self._generate_optimization_input(qcmol_s, compute_spec, constraints)
                        result = self._execute_qcengine(input_data,
                                                        local_options=local_options,
                                                        scf_maxiter=scf_maxiter,
                                                        geometric_maxiter=geometric_maxiter,
                                                        geometric_coordsys=geometric_coordsys,
                                                        geometric_qccnv=geometric_qccnv)
                        
                        if result.success:

                            # TODO: consider if we need to do multiple optimizations per grid point to
                            # get robust results?
                            task_results[gridpoint].append((result.initial_molecule.geometry.flatten().tolist(),
                                                            result.final_molecule.geometry.flatten().tolist(),
                                                            result.energies[-1]))

                            # at the moment, we only hold on to a single result per gridpoint
                            opts[spec_name][gridpoint] = result
                        else:
                            raise QCEngineError(f"QCEngine failure: {result.error.error_message}")

                td_api.update_state(state, task_results)

        # when complete, return output as a JSON-serializable blob
        return opts

    def _execute_qcengine(
            self, input_data, procedure="geometric", local_options=None,
            scf_maxiter=None, geometric_maxiter=None, geometric_coordsys=None,
            geometric_qccnv=None):
        import qcengine

        # inject reset=True fix, set coordsys to 'dlc'
        input_data['keywords']['reset'] = True
        input_data['keywords']['coordsys'] = 'dlc'

        # inject exposed convergence parameters, if specified
        if scf_maxiter is not None:
            input_data['input_specification']['keywords']['maxiter'] = scf_maxiter
        if geometric_maxiter is not None:
            input_data['keywords']['maxiter'] = geometric_maxiter
        if geometric_coordsys is not None:
            input_data['keywords']['coordsys'] = geometric_coordsys
        if geometric_qccnv is not None:
            input_data['keywords']['qccnv'] = geometric_qccnv

        return qcengine.compute_procedure(
                input_data, procedure=procedure, local_options=local_options)

    @staticmethod
    def _generate_optimization_input(qcmol, compute_spec, constraints):
        import ast
        from qcelemental.models import OptimizationInput

        method = compute_spec['method']
        basis = compute_spec['basis']
        program = compute_spec['program']

        # generate optimization inputs
        input_data = OptimizationInput(
                keywords={
                      "coordsys": "dlc",
                      "enforce": 0,
                      "epsilon": 1e-05,
                      "reset": True,
                      "qccnv": True,
                      "molcnv": False,
                      "check": 0,
                      "trust": 0.1,
                      "tmax": 0.3,
                      "maxiter": 300,
                      "convergence_set": "gau",
                      "program": program,
                      "constraints": constraints,
                      },
                extras={},
                protocols={},
                input_specification={
                    "driver": "gradient",
                    "model": {
                      "method": method,
                      "basis": basis,
                    },
                    "keywords": {
                      "maxiter": 200,
                      "scf_properties": [
                        "dipole",
                        "quadrupole",
                        "wiberg_lowdin_indices",
                        "mayer_indices"
                      ]
                    },
                },
                initial_molecule=qcmol)

        return input_data.dict()
