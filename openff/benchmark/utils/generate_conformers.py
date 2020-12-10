from openforcefield.topology import Molecule
from openforcefield.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, OpenEyeToolkitWrapper
import glob
import os
import re
from simtk import unit
import shutil

oetk_loaded = False
for tkw in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:
    if isinstance(tkw, OpenEyeToolkitWrapper):
        oetk_loaded = True
if oetk_loaded:
    GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(OpenEyeToolkitWrapper)


def generate_conformers(input_directory, output_directory, delete_existing=False):

    try:
        os.makedirs(output_directory)
    except OSError:
        if delete_existing:
            shutil.rmtree(output_directory)
            os.makedirs(output_directory)
        else:
            raise Exception(f'Output directory {output_directory} already exists. '
                             'Specify `delete_existing=True` to remove.')

    input_3d_files = glob.glob(os.path.join(input_directory,'*.sdf'))
    input_graph_files = glob.glob(os.path.join(input_directory,'*.smi'))

    group2idx2mols2confs = {}
    # for input_graph_file in input_graph_files:
    #     mapped_smiles = open(input_graph_file).read()
    #     mol = Molecule.from_mapped_smiles(mapped_smiles)
    #     match = re.findall('([a-zA-Z0-9]{3})-([0-9]{5}).smi',
    #                        input_graph_file)
    #     assert len(match) == 1
    #     group_id = match[0][0]
    #     mol_idx = match[0][1]
    #     group2idx2mols2confs[group_id] = group2idx2mols2confs.get(group_id, {})
    #     #try:
    #     assert mol_idx not in group2idx2mols2confs[group_id]
    #     #except:
    #     #    raise Exception(input_graph_file)
    #     group2idx2mols2confs[group_id][mol_idx] = mol
        

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
        if group_id not in group2idx2mols2confs:
            group2idx2mols2confs[group_id] = {}
        if mol_idx not in group2idx2mols2confs[group_id]:
            group2idx2mols2confs[group_id][mol_idx] = {}
        assert conf_id not in group2idx2mols2confs[group_id][mol_idx]
            
        group2idx2mols2confs[group_id][mol_idx][conf_id] = mol
        #raise Exception((group_id, mol_id))

    for group_id in group2idx2mols2confs:
        for mol_idx in group2idx2mols2confs[group_id]:
            n_user_confs = len(group2idx2mols2confs[group_id][mol_idx].keys())
            mol = group2idx2mols2confs[group_id][mol_idx]['00']
            mol.generate_conformers(n_conformers=10,
                                    rms_cutoff=2.*unit.angstrom,
                                    clear_existing=False)
            # Iterate through newly generated confs and add each to the dict as a new mol, to a max of 10
            
            # Get rid of the first conformer, since that was there before
            #mol._conformers = mol._conformers[1:]
            
            for i in range(n_user_confs, mol.n_conformers):
                mol_copy = Molecule(mol)
                mol_copy._conformers = None
                mol_copy.add_conformer(mol.conformers[i])
                group2idx2mols2confs[group_id][mol_idx][f'{i:02d}'] = mol_copy
                
                #filename = f'{group_id}_{mol_idx}_{i:02d}.sdf'
                #full_filename = os.path.join(output_dir
            
            
            
            # Trim down to a max of 10 conformers
            # TODO: Look at pairwise RMSD and pick most diverse if there is a mix
            # of pre-defined and generated conformers
            #group2idx2mols2confs[group_id][mol_idx]._conformers = mol._conformers[:10]
            #raise Exception((group2idx2mols2confs[group_id][mol_idx].n_conformers,
            #                 group2idx2mols2confs))

    # Write outputs
    for group_id in group2idx2mols2confs:
        for mol_idx in group2idx2mols2confs[group_id]:
            for conf_idx in group2idx2mols2confs[group_id][mol_idx]:
                mol_name = f'{group_id}-{int(mol_idx):05d}'
                this_conf = group2idx2mols2confs[group_id][mol_idx][conf_idx]
                this_conf.properties['group_name'] = group_id
                this_conf.properties['molecule_index'] = mol_idx
                this_conf.name = mol_name

                this_conf.properties['conformer_index'] = conf_idx
                out_file_name = f'{this_conf.name}-{int(conf_idx):02d}.sdf'
                this_conf.to_file(os.path.join(output_directory, out_file_name), file_format='sdf')

                #mol_copy = Molecule(group2idx2mols2confs[group_id][mol_idx])

                # Write the graph representation as a fully mapped SMILES
                #with open(os.path.join(output_directory, f'{mol_name}.smi'), 'w') as of:
                #    cmiles = mol_copy.to_smiles(mapped=True)
                #    of.write(cmiles)

                #for conf_index, conformer in enumerate(group2idx2mols2confs[group_id][mol_idx].conformers):
                #    mol_copy._conformers = None
                #raise Exception(conformer)
                #mol_copy.add_conformer(conformer)


        
