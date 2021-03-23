"""
Components for torsiondrives of molecules for FF benchmarking.

These components leverage the QCArchive software stack for execution.

"""

import io
import os
import shutil
import contextlib

from .seasons import SEASONS

class TorsiondriveExecutor:
    def __init__(self):
        self.stop = False

    def execute_torsiondrive_oneshot(self,
            offmol, dihedrals, grid_spacing, dihedral_ranges,
            season, ncores=1, memory=2,
            scf_maxiter=200, geometric_maxiter=300, geometric_coordsys='dlc',
            geometric_qccnv=False):
        """Torsiondrive executor for a single molecule.

        Parameters
        ----------
        dihedrals : List of tuples
            A list of the dihedrals to scan over.
        grid_spacing : List of int
            The grid seperation for each dihedral angle
        dihedral_ranges: (Optional) List of [low, high] pairs
            consistent with launch.py, e.g. [[-120, 120], [-90, 150]]
        energy_decrease_thresh: (Optional) Float
            Threshold of an energy decrease to triggle activate new grid point. Default is 1e-5
        energy_upper_limit: (Optional) Float
            Upper limit of energy relative to current global minimum to spawn new optimization tasks.

        """
        from openff.qcsubmit.factories import OptimizationDatasetFactory
        from torsiondrive import td_api

        local_options={"ncores": ncores,
                       "memory": memory}

        # get qcelemental molecule from openff molecule
        qcmol = offmol.to_qcschema()

        # generate initial torsiondrive state
        td_stdout = io.StringIO()
        with contextlib.redirect_stdout(td_stdout):
            
            # this object gets updated through the whole execution,
            # but it is an implementation detail; we do not return it
            state = td_api.create_initial_state(
                    dihedrals=dihedrals,
                    grid_spacing=grid_spacing,
                    elements=qcmol.symbols,
                    init_coords=qcmol.geometry,
                    dihedral_ranges=dihedral_ranges,
                    )

            # this is the object we care about returning
            opts = dict()

        for spec_name, compute_spec in SEASONS[season].items():
            opts[spec_name] = dict()
            while True:
                next_jobs = td_api.next_jobs_from_state(state)

                # if no next jobs, we are finished
                if len(next_jobs) == 0:
                    break

                task_results = {}
                for gridpoint in next_jobs:
                    task_results[gridpoint] = []

                    qcmol_s = qcmol.copy(deep=True)

                    # perform initial optimization
                    input_data = self._generate_optimization_input(qcmol, compute_spec)
                    result = self._execute_qcengine(input_data,
                                                    local_options=local_options,
                                                    scf_maxiter=scf_maxiter,
                                                    geometric_maxiter=geometric_maxiter,
                                                    geometric_coordsys=geometric_coordsys,
                                                    geometric_qccnv=geometric_qccnv)

                    # TODO: consider if we need to do multiple optimizations per grid point to
                    # get robust results?
                    task_results[gridpoint].append((result['initial_molecule'].geometry,
                                                    result['final_molecule'].geometry,
                                                    result['energies'][-1]))

                    opts[spec_name][gridpoint] = result

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
    def _generate_optimization_input(qcmol, compute_spec):
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
                      "qccnv": False,
                      "molcnv": False,
                      "check": 0,
                      "trust": 0.1,
                      "tmax": 0.3,
                      "maxiter": 300,
                      "convergence_set": "gau",
                      "program": program
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

    def execute_torsiondrive_old(self):
        """
        Author : Josh Horton (@jthorton)
        Posted with permission
        Small additions by David Dotson (@dotsdl)
        """

        from typing import Tuple, List, Optional
        
        import click
        import time
        from pprint import pprint
        
        from qcsubmit.datasets import TorsiondriveDataset
        from qcsubmit.factories import TorsiondriveDatasetFactory
        from qcsubmit.common_structures import QCSpec
        from qcsubmit.results import TorsionDriveCollectionResult
        
        from qcfractal.interface import FractalClient
        from qcfractal import FractalSnowflakeHandler
        from qcengine import list_available_programs
        
        from openforcefield.topology import Molecule
        
        
        def run_torsiondrive(molecule: Molecule, dihedral: Tuple[int, int, int, int], qcspec: QCSpec, workers: int = 2):
            """ This is the main run function of the script witch takes the molecule and
            the dihedral and runs the torsiondrive at the given spec using qcfractal.
        
            Parameters
            ----------
            molecule: The openforcefield molecule instance
            dihedral: The 0 indexed dihedral which should be driven
            qcspec: The QCspecification which should be used to compute the torsiondirve, note the environment must have the program installed or it will hang.
        
            """
            print("building the dataset for submission ...")
            dataset = TorsiondriveDataset()
            factory = TorsiondriveDatasetFactory()
        
            attributes = factory.create_cmiles_metadata(molecule=molecule)
        
            dataset.clear_qcspecs()
            dataset.add_qc_spec(**qcspec.dict())
        
            index = molecule.to_smiles()
        
            # add conformers if non supplied
            if molecule.n_conformers <= 1:
                # in total this should be 3 conformers
                molecule.generate_conformers(n_conformers=3, clear_existing=False)
                print(len(molecule.conformers))
        
            # add the molecule to the dataset
            dataset.add_molecule(index=index, molecule=molecule, attributes=attributes, dihedrals=dihedral)
        
            # set the dataset name
            dataset.dataset_name = "Bespoke Torsiondrive"
            dataset.metadata.long_description_url = "https://test.org"
            print("Dataset validated and ready")
            print(f"Starting qcarchive server with {workers} workers.")
            server = FractalSnowflakeHandler(ncores=workers)
            client = server.client()
        
        
            # this is ready to submit
            print("submitting calculation")
            response = dataset.submit(client=client)
        
            print(f"response from calculations {response}")
        
            collection = client.get_collection("TorsionDriveDataset", dataset.dataset_name)
        
            # now start to monitor the task.
            while True:
                print("waiting for torsiondrive tasks to finish...")
                #server.await_services()
                #server.await_results()
        
                # check if the task is done
                tasks = get_tasks(client=client, dataset=dataset)
        
                # check that the torsiondrives are done else restart them
                task = tasks[0]
                if task.status.value == "ERROR":
                    print(f"Torsiondrive found with error restarting job {task}")
                    
                    # get the optimization ids
                    td_opts = [optimization for optimization in task.optimization_history.values()]
                    opt_restarts = [opt.id for opt in td_opts if opt.status == "ERROR"]
                    error_cycle_optimizations(client=client, optimization_ids=opt_restarts)
                    error_cycle_torsions(client=client, torsiondrive_ids=[task.id, ])
                
                elif task.status in ("INCOMPLETE", "RUNNING"):
                    running = [task for task in tasks if task.status == 'RUNNING']
                    incomplete = [task for task in tasks if task.status == 'INCOMPLETE']
                    print(f"incomplete:\t{len(incomplete)}")
                    print(f"running:\t{len(running)}")
                    time.sleep(10)
                    continue
                
                elif task.status == "COMPLETE":
                    print("torsiondrive complete extracting results.")
        
                    # the task is complete gather for collection
                    result = TorsionDriveCollectionResult.from_server(
                            client=client,
                            spec_name="default",
                            dataset_name=dataset.dataset_name,
                            include_trajectory=False,
                            final_molecule_only=True)
                    break
            
            # now the task is done lets shut down the server
            server.stop()
            print("exporting results to result.json")
            result.export_results("result.json")
        
            return result
        
        
        def build_qcspec(method: str, program: str, basis: Optional[str] = None) -> QCSpec:
            """ Build the qcspec for the user provided inputs.
            """
            spec = QCSpec(method=method, basis=basis, program=program, spec_name="default", spec_description="Bespoke torsiondrive spec.")
            
            return spec
        
        
        def get_tasks(client: FractalClient, dataset: TorsiondriveDataset):
            """ Get a task from the torsion drive dataset.
            """
            tasks = []
            ds = client.get_collection(dataset.dataset_type, dataset.dataset_name)
            ids = dataset.dataset.keys()
        
            for index in ids:
                record = ds.get_record(index, "default")
                tasks.append(record)
            
            return tasks
        
        
        def error_cycle_torsions(client: FractalClient, torsiondrive_ids: List[int]) -> bool:
            """ Error cycle the given torsion ids.
            """
            for td in torsiondrive_ids:
                client.modify_services(operation="restart", procedure_id=td)
            
            return True 
        
        
        def error_cycle_optimizations(client: FractalClient, optimization_ids: List[int]) -> bool:
            """ Error cycle the given optimization ids.
            """
            for opt in optimization_ids:
                client.modify_tasks(operation="restart", base_result=opt)
            
            return True
        
        
        @click.command()
        @click.option("--molecule", type=str, required=True, help="The name of the input file containing the molecule.")
        @click.option("--dihedral", required=True, type=(int, int, int, int), help="The zero indexed torsion tuple")
        @click.option("--method", type=str, default="ani2x")
        @click.option("--basis", type=str, default=None)
        @click.option("--program", type=click.Choice(["torchani", "xtb", "openmm", "psi4", "rdkit"], case_sensitive=False), default="torchani")
        @click.option("--workers", type=int, default=2)
        def main(molecule, dihedral, method, basis, program, workers):
            """ Run a torsiondrive on the selected dihedral in the given molecule, using the specified compute spec.
            """
            # check that we have the specified program
            # error if QCEngine can't find it
            if program.lower() not in list_available_programs():
                raise ValueError(f"Program '{program}' not found in this environment")
        
            # load the molecule
            molecule = Molecule.from_file(molecule)
        
            # if we have multipule conformers they will be seperate molecules so combine
            if isinstance(molecule, list):
                mol = molecule.pop()
                for m in molecule:
                    mol.add_conformer(m.conformer[0])
                molecule = mol
            
            qcspec = build_qcspec(method=method, basis=basis, program=program)
            
            # now run the main function 
            print(dihedral)
            result = run_torsiondrive(molecule=molecule, dihedral=[dihedral, ], qcspec=qcspec, workers=workers)
            
            return result
