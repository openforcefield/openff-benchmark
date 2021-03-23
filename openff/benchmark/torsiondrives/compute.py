"""
Components for torsiondrives of molecules for FF benchmarking.

These components leverage the QCArchive software stack for execution.

"""

import os
import shutil

class TorsiondriveExecutor:

    def __init__(self):
        self.stop = False

    def execute_torsiondrive_oneshot(self,
            offmol, season, ncores=1, memory=2):
        """Torsiondrive executor for a single molecule.

        """
        from torsiondrive import td_api

        # get qcelemental molecule from openff molecule


        # generate initial torsiondrive state
        initial_state = td_api.create_initial_state(
                dihedrals=dihedrals,
                grid_spacing=grid_spacing,
                elements=elements,
                init_coords=init_coords,
                dihedral_ranges=dihedral_ranges,
                )

        # perform initial optimization


        # generate next optimization(s), then execute in a loop


        # when complete, return output as a JSON-serializable blob



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
