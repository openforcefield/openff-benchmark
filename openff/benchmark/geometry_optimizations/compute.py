"""
Components for energy minimization / geometry optimization of molecules for FF benchmarking.

These components leverage the QCArchive software stack for execution.

"""

import os
import shutil
import logging
from collections import defaultdict
logging.disable(logging.WARNING) 

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
    
    def submit_molecules(self, fractal_uri, input_paths, season, dataset_name="Benchmark Optimizations", 
            replace=True):
        """Submit SDF molecules from given directory to the target QCFractal server.
    
        Parameters
        ----------
        fractal_uri : str
            Target QCFractal server URI.
        input_paths : iterable of Path-like
            Paths to SDF files or directories; if directories, all files SDF files in are loaded, recursively.
        season : str
            Benchmark season identifier. Indicates the mix of compute specs to utilize.
        dataset_name : str
            Dataset name to use for submission on the QCFractal server.
    
        """
        from qcsubmit.factories import OptimizationDataset, OptimizationDatasetFactory
    
        # extract molecules from SDF inputs
        mols = mols_from_paths(input_paths)
    
        # TODO: check existence of dataset, consult with JH on best way to add to existing
    
        # create OptimizationDataset with QCSubmit
        ds = OptimizationDataset(dataset_name=dataset_name)
        factory = OptimizationDatasetFactory()
    
        for mol in mols:
            id = self._mol_to_id(mol)
    
            attributes = factory.create_cmiles_metadata(mol)
            ds.add_molecule(index=id, molecule=mol, attributes=attributes)
    
        ds.qc_specifications = SEASONS[season]
    
        ds.metadata.long_description_url = "https://localhost.local/null"
    
        logging.info("Submitting...")
        client = FractalClient(fractal_uri, verify=False)
        ds.submit(client=client)
        logging.info("Submitted!")
    
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
    
    def export_molecule_data(self, fractal_uri, output_directory, dataset_name="Benchmark Optimizations",
                             delete_existing=False, keep_existing=False):
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
    
        """
        from openforcefield.topology.molecule import unit
        import numpy as np
        import pint
    
        punit = pint.UnitRegistry()
    
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
                geometry = unit.Quantity(
                        np.array(final_qcmol.geometry, np.float), unit.bohr
                    )
                mol._add_conformer(geometry.in_units_of(unit.angstrom))
    
                # set molecule metadata
                mol.name = output_id
    
                (mol.properties['group_name'],
                 mol.properties['molecule_index'],
                 mol.properties['conformer_index']) = output_id.split('-')
    
                # SDF key-value pairs should be used for method, basis, program, provenance, `openff-benchmark` version
                mol.properties['method'] = optentspec.qc_spec.method
                mol.properties['basis'] = optentspec.qc_spec.basis
                mol.properties['program'] = optentspec.qc_spec.program
    
                # SDF key-value pairs also used for final energies
                mol.properties['initial_energy'] = (opt.energies[0] * punit.hartree * punit.avogadro_constant).to('kcal/mol')
                mol.properties['final_energy'] = (opt.energies[-1] * punit.hartree * punit.avogadro_constant).to('kcal/mol')
    
                mol.to_file(outfile, file_format='sdf')
    
    def get_optimization_status(self, fractal_uri, dataset_name="Benchmark Optimizations", client=None):
        """Get status of optimization for each molecule ID
    
        """
        if client is None:
            client = FractalClient(fractal_uri, verify=False)
    
        optds = client.get_collection("OptimizationDataset", dataset_name)
        optds.status()
        return optds.df.sort_index(ascending=True)

    def set_optimization_priority(self, fractal_uri, priority, dataset_name="Benchmark Optimizations"):
        from qcportal.models.task_models import PriorityEnum

        client = FractalClient(fractal_uri, verify=False)

        optds = client.get_collection("OptimizationDataset", dataset_name)
        optds.status()
        opts = optds.df.values.flatten()

        priority_map = {"high": PriorityEnum.HIGH,
                        "normal": PriorityEnum.NORMAL,
                        "low": PriorityEnum.LOW}

        for opt in opts:
            if opt.status != 'COMPLETE':
                client.modify_tasks(operation='modify', base_result=opt.id, new_priority=priority_map[priority])

    def errorcycle_optimizations(self, fractal_uri, dataset_name="Benchmark Optimizations", client=None,
            compute_specs=None, molids=None):
        """Restart optimizations that have failed.

        Parameters
        ----------
        compute_specs : list
            List of compute specs to error cycle only.
        molids : list
            List of molecule ids to error cycle only.
    
        """
        if client is None:
            client = FractalClient(fractal_uri, verify=False)

        optds = client.get_collection("OptimizationDataset", dataset_name)
        optds.status()

        df = optds.df

        if (molids is not None) and (len(molids) != 0):
            df = df.loc[molids]

        if compute_specs is not None:
            df = df[compute_specs]

        for opt in df.values.flatten():
            if opt.status == 'ERROR':
                print(f"Restarted ERRORed optimization `{opt.id}`")
                client.modify_tasks(operation='restart', base_result=opt.id)

    
    def execute_optimization_from_server(self, molecule_id, push_complete=False):
        """Execute optimization from the given molecule locally on this host.
    
        """
        pass
    
    def _source_specs_output_paths(self, input_paths, specs, output_directory):
        mols = mols_from_paths(input_paths, sourcefile_keys=True)
    
        in_out_path_map = defaultdict(list)
        for sourcefile, mol in mols.items():
            id = self._mol_to_id(mol)
            for spec in specs:
                in_out_path_map[sourcefile].append("{}.sdf".format(os.path.join(output_directory, spec, id)))
    
        return in_out_path_map
    
    def execute_optimization_from_molecules(
            self, input_paths, output_directory, season, ncores=1,
            dataset_name='Benchmark Scratch', delete_existing=False, keep_existing=False):
        """Execute optimization from the given SDF molecules locally on this host.
    
        """
        from time import sleep
        import psutil

        from tqdm import trange
        from qcfractal import FractalSnowflake
    
        # fail early if output_directory already exists and we aren't deleting it
        if os.path.isdir(output_directory):
            if delete_existing:
                shutil.rmtree(output_directory)
            elif keep_existing:
                pass
            else:
                raise Exception(f'Output directory {output_directory} already exists. '
                                 'Specify `delete_existing=True` to remove, or `keep_existing=True` to tolerate')
    
        # start up Snowflake
        server = FractalSnowflake(max_workers=ncores)

        client = server.client()
        fractal_uri = server.get_address()
    
        # get paths to submit, using output directory contents to inform choice
        # for the given specs, if *any* expected output files are not present, we submit corresponding input file
        if keep_existing:
            in_out_path_map = self._source_specs_output_paths(input_paths, SEASONS[season], output_directory)
    
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
            #total = df.size
            #print(f"{complete}/{total} COMPLETE", end='\r')
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
