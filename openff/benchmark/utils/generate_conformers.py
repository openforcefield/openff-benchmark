from openforcefield.topology import Molecule
from openforcefield.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, OpenEyeToolkitWrapper
import glob
import os
import re
from simtk import unit

oetk_loaded = False
for tkw in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:
    if isinstance(tkw, OpenEyeToolkitWrapper):
        oetk_loaded = True
if oetk_loaded:
    GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(OpenEyeToolkitWrapper)


def generate_conformers(input_directory, output_directory):
    if os.path.exists(output_directory):
        raise Exception(f'Output directory {output_directory} already exists. '
                         'The user must delete this manually. '
                         'This script will not overwrite it.')
    


    input_3d_files = glob.glob(os.path.join(input_directory,'*.sdf'))
    input_graph_files = glob.glob(os.path.join(input_directory,'*.smi'))

    group2idx2mols = {}
    for input_graph_file in input_graph_files:
        mapped_smiles = open(input_graph_file).read()
        mol = Molecule.from_mapped_smiles(mapped_smiles)
        match = re.findall('([a-zA-Z0-9]{3})-([0-9]{5}).smi',
                           input_graph_file)
        assert len(match) == 1
        group_id = match[0][0]
        mol_idx = match[0][1]
        group2idx2mols[group_id] = group2idx2mols.get(group_id, {})
        #try:
        assert mol_idx not in group2idx2mols[group_id]
        #except:
        #    raise Exception(input_graph_file)
        group2idx2mols[group_id][mol_idx] = mol
        

    for input_3d_file in input_3d_files:
        #try:
        # We have to allow undefined stereo here because of some silliness with
        # RDKitToolkitWrapper percieveing stereo around terminal ethene groups
        mol = Molecule.from_file(input_3d_file, allow_undefined_stereo=True)
        #except:
        #    raise Exception(input_3d_file)
        match = re.findall('([a-zA-Z0-9]{3})-([0-9]{5})-([0-9]{2}).sdf',
                           input_3d_file)
        assert len(match) == 1
        group_id = match[0][0]
        mol_idx = match[0][1]
        conf_id = match[0][2]
        group2idx2mols[group_id][mol_idx].add_conformer(mol.conformers[0])
        #raise Exception((group_id, mol_id))

    for group_id in group2idx2mols:
        for mol_idx in group2idx2mols[group_id]:
            mol = group2idx2mols[group_id][mol_idx]
            mol.generate_conformers(n_conformers=10,
                                    rms_cutoff=2.*unit.angstrom,
                                    clear_existing=False)
            
            # Trim down to a max of 10 conformers
            # TODO: Look at pairwise RMSD and pick most diverse if there is a mix
            # of pre-defined and generated conformers
            group2idx2mols[group_id][mol_idx]._conformers = mol._conformers[:10]
            #raise Exception((group2idx2mols[group_id][mol_idx].n_conformers,
            #                 group2idx2mols))

    os.makedirs(output_directory)
    
    for group_id in group2idx2mols:
        for mol_idx in group2idx2mols[group_id]:
            mol_name = f'{group_id}-{int(mol_idx):05d}'
            group2idx2mols[group_id][mol_idx].properties['group_name'] = group_id
            group2idx2mols[group_id][mol_idx].properties['molecule_index'] = mol_idx
            group2idx2mols[group_id][mol_idx].name = mol_name
            
            mol_copy = Molecule(group2idx2mols[group_id][mol_idx])
            
            # Write the graph representation as a fully mapped SMILES
            with open(os.path.join(output_directory, f'{mol_name}.smi'), 'w') as of:
                cmiles = mol_copy.to_smiles(mapped=True)
                of.write(cmiles)
            
            for conf_index, conformer in enumerate(group2idx2mols[group_id][mol_idx].conformers):
                mol_copy._conformers = None
                #raise Exception(conformer)
                mol_copy.add_conformer(conformer)
                mol_copy.properties['conformer_index'] = conf_index
                out_file_name = f'{mol_copy.name}-{conf_index:02d}.sdf'
                mol_copy.to_file(os.path.join(output_directory, out_file_name), file_format='sdf')
            

        
