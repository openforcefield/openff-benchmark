"""
Components for energy minimization / geometry optimization of molecules for FF benchmarking.

These components leverage the QCArchive software stack for execution.

"""

import os
import shutil
import logging
from collections import defaultdict

# DO NOT USE; MESSES UP LOGGING 
# consider overriding the loglevel of any particularly noisy loggers individual
#logging.disable(logging.WARNING) 

from qcportal import FractalClient

from .seasons import SEASONS
from ..utils.io import mols_from_paths


class OptimizationExecutor:

    def __init__(self):
        self.stop = False

    @staticmethod
    def _mol_to_id(mol):
        id = "{org}-{molecule:05}-{conformer:02}".format(
                org=mol.properties['group_name'],
                molecule=mol.properties['molecule_index'],
                conformer=mol.properties['conformer_index'])
        return id
    
    def submit_molecules(self, fractal_uri, input_paths, season,
            dataset_name, recursive=False):
        """Submit SDF molecules from given directory to the target QCFractal server.
    
        Parameters
        ----------
        fractal_uri : str
            Target QCFractal server URI.
        input_paths : iterable of Path-like
            Paths to SDF files or directories; for directories, all SDF files are loaded.
        season : str
            Benchmark season identifier. Indicates the mix of compute specs to utilize.
        dataset_name : str
            Dataset name to use for submission on the QCFractal server.
        recursive : bool
            If True, recursively load SDFs from any directories given in `input_paths`.
    
        """
        from openff.qcsubmit.factories import OptimizationDataset, OptimizationDatasetFactory
    
        # extract molecules from SDF inputs
        mols = mols_from_paths(input_paths, recursive=recursive)
    
        # create OptimizationDataset with QCSubmit
        ds = OptimizationDataset(dataset_name=dataset_name)
        factory = OptimizationDatasetFactory()
    
        for mol in mols:
            id = self._mol_to_id(mol)

            attributes = factory.create_cmiles_metadata(mol)
            ds.add_molecule(index=id, molecule=mol, attributes=attributes)
    
        ds.qc_specifications = SEASONS[season]
    
        ds.metadata.long_description_url = "https://localhost.local/null"

        # add in known modifications to `OptimizationDataset` defaults
        ds.optimization_procedure.coordsys = 'dlc'
        ds.optimization_procedure.reset = True
    
        print("Submitting...")
        client = FractalClient(fractal_uri, verify=False)
        ds.submit(verbose=True, client=client)
        print("Submitted!")
    
    @staticmethod
    def _mol_from_qcserver(record):
        from openforcefield.topology import Molecule
        from openforcefield.utils import toolkits
        
        # make sure we deregister OpenEye, if it is present
        try:
            toolkits.GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(toolkits.OpenEyeToolkitWrapper)
        except toolkits.ToolkitUnavailableException:
            pass
    
        return Molecule.from_qcschema(record)
    
    def export_molecule_data(self, fractal_uri, output_directory, dataset_name,
                             delete_existing=False, keep_existing=True):
        """Export all molecule data from target QCFractal server to the given directory.
    
        Parameters
        ----------
        fractal_uri : str
            Target QCFractal server URI.
        output_directory : str
            Directory path to deposit exported data.
        dataset_name : str
            Dataset name to extract from the QCFractal server.
        delete_existing : bool (False)
            If True, delete existing directory if present.
        keep_existing : bool (True)
            If True, keep existing files in export directory.
            Files corresponding to server data will not be re-exported.
            Relies *only* on filepaths of existing files for determining match.
    
        """
        # get dataset
        client = FractalClient(fractal_uri, verify=False)
        optds = client.get_collection("OptimizationDataset", dataset_name)
        optds.status()
    
        try:
            os.makedirs(output_directory)
        except OSError:
            if delete_existing:
                shutil.rmtree(output_directory)
            elif keep_existing:
                pass
            else:
                raise Exception(f'Output directory {output_directory} already exists. '
                                 'Specify `delete_existing=True` to remove, or `keep_existing=True` to tolerate')
    
        # for each compute spec, create a folder in the output directory
        # deposit SDF giving final molecule, energy
        specs = optds.list_specifications().index.tolist()
        for spec in specs:
            os.makedirs(os.path.join(output_directory, spec), exist_ok=True)
            optentspec = optds.get_specification(spec)
    
            records = optds.data.dict()['records']
    
            for id, opt in optds.df[spec].iteritems():
    
                # skip incomplete cases
                if opt.status != 'COMPLETE':
                    continue
    
                # fix to ensure output fidelity of ids; losing 02 padding on conformer
                org, molecule, conformer = id.split('-')
                output_id = "{org}-{molecule:05}-{conformer:02}".format(org=org,
                                                                        molecule=int(molecule),
                                                                        conformer=int(conformer))
    
                # subfolders for each compute spec, files named according to molecule ids
                outfile = "{}.sdf".format(
                        os.path.join(output_directory, spec, output_id))
    
                # if we did not delete everything at the start and the path already exists,
                # skip this one; reduces processing and writes to filesystem
                if (not delete_existing) and os.path.exists(outfile):
                    continue
    
                mol = self._mol_from_qcserver(records[id.lower()])
    
                # set conformer as final, optimized geometry
                final_qcmol = opt.get_final_molecule()
                mol = self._process_final_mol(output_id,
                                              mol,
                                              final_qcmol,
                                              optentspec.qc_spec.method,
                                              optentspec.qc_spec.basis,
                                              optentspec.qc_spec.program,
                                              opt.energies)

                mol.to_file(outfile, file_format='sdf')

    def get_optimization_status(self, fractal_uri, dataset_name, client=None,
            compute_specs=None, molids=None):
        """Get status of optimization for each molecule ID.
    
        """
        if client is None:
            client = FractalClient(fractal_uri, verify=False)

        optds = client.get_collection("OptimizationDataset", dataset_name)
        optds.status()
        
        df = optds.df.sort_index(ascending=True)

        if (molids is not None) and (len(molids) != 0):
            df = df.loc[list(molids)]

        if compute_specs is not None:
            df = df[compute_specs]

        return df

    def set_optimization_priority(self, fractal_uri, priority, dataset_name):
        from qcportal.models.task_models import PriorityEnum

        client = FractalClient(fractal_uri, verify=False)

        optds = client.get_collection("OptimizationDataset", dataset_name)
        optds.status()
        opts = optds.df.values.flatten()

        priority_map = {"high": PriorityEnum.HIGH,
                        "normal": PriorityEnum.NORMAL,
                        "low": PriorityEnum.LOW}

        client.modify_tasks(operation='modify',
                            base_result=[opt.id for opt in opts if opt.status != 'COMPLETE'],
                            new_priority=priority_map[priority])

    def errorcycle_optimizations(self, fractal_uri, dataset_name, client=None,
            compute_specs=None, molids=None):
        """Restart optimizations that have failed.

        Parameters
        ----------
        compute_specs : iterable 
            Iterable of compute spec names to error cycle only.
        molids : iterable 
            Iterable of molecule ids to error cycle only.
    
        """
        if client is None:
            client = FractalClient(fractal_uri, verify=False)

        optds = client.get_collection("OptimizationDataset", dataset_name)
        optds.status()

        df = optds.df

        if (molids is not None) and (len(molids) != 0):
            df = df.loc[list(molids)]

        if compute_specs is not None:
            df = df[compute_specs]

        for opt in df.values.flatten():
            if opt.status == 'ERROR':
                print(f"Restarted ERRORed optimization `{opt.id}`")
                client.modify_tasks(operation='restart', base_result=opt.id)
    
    def get_optimization_from_server(self, fractal_uri, dataset_name, client=None,
            compute_specs=None, molids=None):
        """Get full optimization data from the given molecules.
    
        """

        if client is None:
            client = FractalClient(fractal_uri, verify=False)

        optds = client.get_collection("OptimizationDataset", dataset_name)
        optds.status()

        df = optds.df

        if (molids is not None) and (len(molids) != 0):
            df = df.loc[list(molids)]

        if compute_specs is not None:
            df = df[compute_specs]

        results = []
        for opt in df.values.flatten():
            results.append(opt)

        return results

    def get_optimization_tracebacks(self, fractal_uri, dataset_name, client=None,
            compute_specs=None, molids=None):

        if client is None:
            client = FractalClient(fractal_uri, verify=False)

        optds = client.get_collection("OptimizationDataset", dataset_name)
        optds.status()
        
        df = optds.df.sort_index(ascending=True)

        if (molids is not None) and (len(molids) != 0):
            df = df.loc[list(molids)]

        if compute_specs is not None:
            df = df[compute_specs]

        errors = df.applymap(lambda x: x.get_error().error_message if x.status == 'ERROR' else None)

        # filter down to only those rows with errors
        errors = errors.dropna(how='all')
            
        return errors

    def list_optimization_datasets(self, fractal_uri, client=None):
        if client is None:
            client = FractalClient(fractal_uri, verify=False)

        datasets = client.list_collections('OptimizationDataset')

        return datasets.reset_index()['name'].to_list()

    def delete_optimization_datasets(self, fractal_uri, dataset_names, client=None):
        if client is None:
            client = FractalClient(fractal_uri, verify=False)

        for dataset_name in dataset_names:
            client.delete_collection('OptimizationDataset', dataset_name)
    
    def _source_specs_output_paths(self, input_paths, specs, output_directory, recursive):
        mols = mols_from_paths(input_paths, sourcefile_keys=True, recursive=recursive)
    
        in_out_path_map = defaultdict(list)
        for sourcefile, mol in mols.items():
            id = self._mol_to_id(mol)
            for spec in specs:
                in_out_path_map[sourcefile].append("{}.sdf".format(os.path.join(output_directory, spec, id)))
    
        return in_out_path_map
    
    def execute_optimization_from_server(
            self, fractal_uri, dataset_name, output_directory=None, ncores=1,
            memory=2, client=None, compute_specs=None, molids=None,
            scf_maxiter=200, geometric_maxiter=300, geometric_coordsys='dlc',
            geometric_qccnv=False):
        """Execute optimization from the given molecule locally on this host.

        Will not send results back to the server; this is purely for debugging.

        TODO: make this send results back to server using same API as manager does.
              then merge with `execute_...` above.
    
        """
        if client is None:
            client = FractalClient(fractal_uri, verify=False)

        optds = client.get_collection("OptimizationDataset", dataset_name)
        optds.status()

        df = optds.df

        if (molids is not None) and (len(molids) != 0):
            df = df.loc[list(molids)]

        if compute_specs is not None:
            df = df[compute_specs]

        local_options={"ncores": ncores,
                       "memory": memory}

        results = []
        for spec_name in df:

            if output_directory is not None:
                os.makedirs(os.path.join(output_directory, spec_name), exist_ok=True)

            print("Processing spec: '{}'".format(spec_name))
            for id, opt in df[spec_name].iteritems():

                # fix to ensure output fidelity of ids; losing 02 padding on conformer
                org, molecule, conformer = id.split('-')
                output_id = "{org}-{molecule:05}-{conformer:02}".format(org=org,
                                                                        molecule=int(molecule),
                                                                        conformer=int(conformer))

                # subfolders for each compute spec, files named according to molecule ids
                if output_directory is not None:
                    outfile = "{}".format(
                            os.path.join(output_directory, spec_name, output_id))
    
                print("... '{}'".format(id))
                task = client.query_tasks(base_result=opt.id)[0]

                # execute optimization
                result = self._execute_qcengine(task.spec.args[0],
                                                local_options=local_options,
                                                scf_maxiter=scf_maxiter,
                                                geometric_maxiter=geometric_maxiter,
                                                geometric_coordsys=geometric_coordsys,
                                                geometric_qccnv=geometric_qccnv)

                if result.success:

                    # process final molecule
                    final_molecule = self._process_optimization_result(output_id, result)
                    final_molecule.to_file("{}.sdf".format(outfile), file_format='sdf')
                else:
                    print("Optimization failed for '{}'; check JSON results output".format(output_id))

                with open("{}.json".format(outfile), 'w') as f:
                    f.write(result.json())

                results.append(result)

        return results

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

    def _process_optimization_result(self, output_id, result):
        from openforcefield.topology import Molecule

        cmiles = result.initial_molecule.extras['canonical_isomeric_explicit_hydrogen_mapped_smiles']
        
        final_molecule = Molecule.from_mapped_smiles(cmiles)
        final_qcmol = result.final_molecule
        
        method = result.input_specification.model.method
        basis = result.input_specification.model.basis
        program = result.keywords['program']

        final_molecule = self._process_final_mol(output_id,
                                       final_molecule,
                                       final_qcmol,
                                       method,
                                       basis,
                                       program,
                                       result.energies)

        return final_molecule

    @staticmethod
    def _process_final_mol(output_id, offmol, qcmol, method, basis, program, energies):
        from openforcefield.topology.molecule import unit
        import numpy as np
        import pint
    
        punit = pint.UnitRegistry()
    
        # set conformer qcmol geometry
        geometry = unit.Quantity(
                np.array(qcmol.geometry, np.float), unit.bohr
            )
        offmol._add_conformer(geometry.in_units_of(unit.angstrom))
    
        # set molecule metadata
        offmol.name = output_id
    
        (offmol.properties['group_name'],
         offmol.properties['molecule_index'],
         offmol.properties['conformer_index']) = output_id.split('-')
    
        # SDF key-value pairs should be used for method, basis, program, provenance, `openff-benchmark` version
        offmol.properties['method'] = method
        offmol.properties['basis'] = basis
        offmol.properties['program'] = program
    
        # SDF key-value pairs also used for final energies
        offmol.properties['initial_energy'] = (energies[0] * punit.hartree * punit.avogadro_constant).to('kcal/mol')
        offmol.properties['final_energy'] = (energies[-1] * punit.hartree * punit.avogadro_constant).to('kcal/mol')

        return offmol
    
    def execute_optimization_from_molecules(
            self, input_paths, output_directory, season, ncores=1, memory=2, 
            delete_existing=False, keep_existing=True, recursive=False,
            scf_maxiter=200, geometric_maxiter=300, geometric_coordsys='dlc',
            geometric_qccnv=False):
        """Execute optimizations from the given SDF molecules locally on this host.

        Optimizations are performed in series for the molecules given,
        with `ncores` and `memory` setting the resource constraints each optimization.

        Parameters
        ----------
        input_paths : iterable of Path-like
            Paths to SDF files or directories; if directories, all files SDF files in are loaded, recursively.
        output_directory : str
            Directory path to deposit exported data.
        season : str
            Benchmark season identifier. Indicates the mix of compute specs to utilize.
        ncores : int
            Number of concurrent cores to use for each optimization.
        memory : float
            Amount of memory in GiB to allow for each optimization.
        delete_existing : bool (False)
            If True, delete existing directory if present.
        keep_existing : bool (True)
            If True, keep existing files in export directory.
            Files corresponding to server data will not be re-exported.
            Relies *only* on filepaths of existing files for determining match.
        recursive : bool
            If True, recursively load SDFs from any directories given in `input_paths`.
    
        """
        from datetime import datetime
        import json

        import numpy as np
        import pint
    
        from openff.qcsubmit.factories import OptimizationDatasetFactory
        from openforcefield.topology import Molecule
        from openforcefield.topology.molecule import unit

        punit = pint.UnitRegistry()
    
        # fail early if output_directory already exists and we aren't deleting it
        if os.path.isdir(output_directory):
            if delete_existing:
                shutil.rmtree(output_directory)
            elif keep_existing:
                pass
            else:
                raise Exception(f'Output directory {output_directory} already exists. '
                                 'Specify `delete_existing=True` to remove, or `keep_existing=True` to tolerate')

        # get paths to submit, using output directory contents to inform choice
        # for the given specs, if *any* expected output files are not present, we submit corresponding input file
        if keep_existing:
            in_out_path_map = self._source_specs_output_paths(
                    input_paths, SEASONS[season], output_directory, recursive=recursive)
    
            input_paths = []
            for input_file, output_files in in_out_path_map.items():
                if not all(map(os.path.exists, output_files)):
                    input_paths.append(input_file)


        # extract molecules from SDF inputs
        mols = mols_from_paths(input_paths, recursive=recursive)
        factory = OptimizationDatasetFactory()

        local_options={"ncores": ncores,
                       "memory": memory}

        results = []
        for spec_name, compute_spec in SEASONS[season].items():
            print("Processing spec: '{}'".format(spec_name))

            os.makedirs(os.path.join(output_directory, spec_name), exist_ok=True)

            for mol in mols:
                id = self._mol_to_id(mol)

                # fix to ensure output fidelity of ids; losing 02 padding on conformer
                org, molecule, conformer = id.split('-')
                output_id = "{org}-{molecule:05}-{conformer:02}".format(org=org,
                                                                        molecule=int(molecule),
                                                                        conformer=int(conformer))

                # subfolders for each compute spec, files named according to molecule ids
                outfile = "{}".format(
                        os.path.join(output_directory, spec_name, output_id))
    
                print("... '{}'".format(id))

                input_data = self._generate_optimization_input(mol, compute_spec, factory)

                # execute optimization
                start_dt = datetime.utcnow()
                result = self._execute_qcengine(input_data,
                                                local_options=local_options,
                                                scf_maxiter=scf_maxiter,
                                                geometric_maxiter=geometric_maxiter,
                                                geometric_coordsys=geometric_coordsys,
                                                geometric_qccnv=geometric_qccnv)

                end_dt = datetime.utcnow()
                if result.success:
                    # process final molecule
                    final_molecule = self._process_optimization_result(output_id, result)
                    final_molecule.to_file("{}.sdf".format(outfile), file_format='sdf')
                else:
                    print("Optimization failed for '{}'; check JSON results output".format(output_id))

                with open("{}.json".format(outfile), 'w') as f:
                    f.write(result.json())

                with open("{}.perf.json".format(outfile), 'w') as f:
                    json.dump({'start': start_dt.isoformat(), 'end': end_dt.isoformat()}, f)

                results.append(result)

        return results

    @staticmethod
    def _generate_optimization_input(mol, compute_spec, factory):
        import ast
        from qcelemental.models import OptimizationInput

        # TODO: bug report in openff where `atom_map` is a string
        if isinstance(mol.properties.get('atom_map'), str):
            mol.properties['atom_map'] = ast.literal_eval(mol.properties['atom_map'])
    
        attributes = factory.create_cmiles_metadata(mol)
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
                initial_molecule=mol.to_qcschema(extras=attributes)
            )

        return input_data.dict()

    def execute_with_snowflake(
            self, input_paths, output_directory, season, ncores=None, memory=None, 
            dataset_name='Benchmark Scratch', delete_existing=False, keep_existing=True,
            recursive=False):
        """Execute optimizations from the given SDF molecules locally on this host.

        Optimizations are performed in series for the molecules given,
        with `ncores` and `memory` setting the resource constraints each optimization.

        Parameters
        ----------
        input_paths : iterable of Path-like
            Paths to SDF files or directories; if directories, all files SDF files in are loaded, recursively.
        output_directory : str
            Directory path to deposit exported data.
        season : str
            Benchmark season identifier. Indicates the mix of compute specs to utilize.
        ncores : int
            Number of concurrent cores to use for each optimization.
        dataset_name : str
            Dataset name to extract from the QCFractal server.
        delete_existing : bool (False)
            If True, delete existing directory if present.
        keep_existing : bool (True)
            If True, keep existing files in export directory.
            Files corresponding to server data will not be re-exported.
            Relies *only* on filepaths of existing files for determining match.
        recursive : bool
            If True, recursively load SDFs from any directories given in `input_paths`.
    
        """
        from openff.qcsubmit.factories import OptimizationDatasetFactory
    
        # fail early if output_directory already exists and we aren't deleting it
        if os.path.isdir(output_directory):
            if delete_existing:
                shutil.rmtree(output_directory)
            elif keep_existing:
                pass
            else:
                raise Exception(f'Output directory {output_directory} already exists. '
                                 'Specify `delete_existing=True` to remove, or `keep_existing=True` to tolerate')

        # get paths to submit, using output directory contents to inform choice
        # for the given specs, if *any* expected output files are not present, we submit corresponding input file
        if keep_existing:
            in_out_path_map = self._source_specs_output_paths(
                    input_paths, SEASONS[season], output_directory, recursive=recursive)
    
            input_paths = []
            for input_file, output_files in in_out_path_map.items():
                if not all(map(os.path.exists, output_files)):
                    input_paths.append(input_file)


        from time import sleep
        import psutil

        from tqdm import trange
        from qcfractal import FractalSnowflake

        # start up Snowflake
        server = FractalSnowflake(max_workers=ncores)

        client = server.client()
        fractal_uri = server.get_address()
    
        # get paths to submit, using output directory contents to inform choice
        # for the given specs, if *any* expected output files are not present, we submit corresponding input file
        if keep_existing:
            in_out_path_map = self._source_specs_output_paths(
                    input_paths, SEASONS[season], output_directory, recursive=recursive)
    
            input_paths = []
            for input_file, output_files in in_out_path_map.items():
                if not all(map(os.path.exists, output_files)):
                    input_paths.append(input_file)
    
        # submit molecules
        self.submit_molecules(fractal_uri, input_paths, season, dataset_name=dataset_name)

        df = self.get_optimization_status(fractal_uri, dataset_name, client=client)
        progbar = trange(df.size)
    
        complete = 0
        while not self.stop:
            df = self.get_optimization_status(fractal_uri, dataset_name, client=client)
    
            # write out what we can
            self.export_molecule_data(fractal_uri, output_directory, dataset_name=dataset_name,
                    delete_existing=False, keep_existing=True)
    
            # break if complete
            complete_i = df.applymap(lambda x: x.status.value == 'COMPLETE').sum().sum()
            progbar.update(complete_i-complete)
            complete = complete_i
            if complete == df.size:
                break
            sleep(10)
    
        # one final export, just in case some completed since last write
        self.export_molecule_data(fractal_uri, output_directory, dataset_name=dataset_name, 
                delete_existing=False, keep_existing=True)

        # stop the server and all its processes
        #parent = psutil.Process(server._qcfractal_proc.pid)
        #for child in parent.children(recursive=True):
        #    child.kill()
        #parent.kill()

        server.stop()
