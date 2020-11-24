from openforcefield.topology import Molecule
from openforcefield.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, OpenEyeToolkitWrapper
import glob
import os

oetk_loaded = False
for tkw in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:
    if isinstance(tkw, OpenEyeToolkitWrapper):
        oetk_loaded = True
if oetk_loaded:
    GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(OpenEyeToolkitWrapper)
    
def validate_and_assign(input_graph_files,
                        input_3d_files,
                        output_directory,
                        group_name):
    """
    Load a molecule dataset from SDF, validate it for common 
    issues, and assign it unique identifiers.
    """
    
    # Load all the graph molecules
    #    delete conformers, if present
    #    raise an error on any duplicates
    # Load all the 3d molecules, see if they're identical to any 2d molecules
    #    if not, make a new 2d molecule
    #    if so, add the coordinates as a new molecule (performing an isomorphism)

    if os.path.exists(output_directory):
        raise Exception(f'Output directory {output_directory} already exists. '
                         'The user must delete this manually. '
                         'This script will not overwrite it.')
    
    # Write all molecules to .smi files
    # Write all conformers to 3D sdfs with numerical conformer

    smiles2mol = {}

    # Handle graph molecules
    #molecule_graph_files = glob.glob(input_graph_files)
    
    for molecule_graph_file in input_graph_files:
        # TODO: Have an option for permissive stereochemistry?
        loaded_mols = Molecule.from_file(molecule_graph_file,
                                         file_format='sdf')
        if not isinstance(loaded_mols, list):
            loaded_mols = [loaded_mols]
            
        # Process each graph molecule and check for duplicates
        for mol_index, mol in enumerate(loaded_mols):
            # Sanitize any unwanted original information
            mol.name = None
            keys = list(mol.properties.keys())
            for key in keys:
                mol.properties.pop(key)
            mol.partial_charges = None
            # Since this is a graph molecule input, clear any conformers
            mol._conformers = []
            
            # Keep a record of the context from which this was loaded
            mol.properties['original_file'] = molecule_graph_file
            mol.properties['original_file_index'] = mol_index + 1

            # Check whether this is a duplicate
            smiles = mol.to_smiles()
            if smiles in smiles2mol:
                other_mol = smiles2mol[smiles]
                other_file = other_mol.properties['original_file']
                other_index = other_mol.properties['original_file_index']
                msg = 'Duplicate graph molecule detected:\n'
                msg += f'Molecule {other_index}: {other_mol}\n'
                msg += f'Molecule {mol_index}: {mol}'
                raise Exception(msg)
                
            smiles2mol[smiles] = mol
                 
            
    # Handle 3d molecules
    #molecule_3d_files = glob.glob(input_3d_files)
    for molecule_3d_file in input_3d_files:
        loaded_mols = Molecule.from_file(molecule_3d_file,
                                         file_format='sdf')
        if not isinstance(loaded_mols, list):
            loaded_mols = [loaded_mols]

        for mol_index, mol in enumerate(loaded_mols):
            # Sanitize any information that might already be present
            mol.name = None
            keys = list(mol.properties.keys())
            for key in keys:
                mol.properties.pop(key)
            mol.partial_charges = None
            
            # Keep a record of the context from which this was loaded
            mol.properties['original_file'] = molecule_3d_file
            mol.properties['original_file_index'] = mol_index + 1

            # If this graph molecule IS already known, add this 3d information as a conformer
            smiles = mol.to_smiles()
            if smiles in smiles2mol:
                orig_mol = smiles2mol[smiles]
                _, atom_map = Molecule.are_isomorphic(mol,
                                                      orig_mol,
                                                      return_atom_map=True,
                                                      formal_charge_matching=False,
                                                      aromatic_matching=False,
                                                      #atom_stereochemistry_matching=False,
                                                      #bond_stereochemistry_matching=False,
                                                      )
                reordered_mol = mol.remap(atom_map)
                smiles2mol[smiles].add_conformer(reordered_mol.conformers[0])
                
            # If this graph molecule ISN'T already known, then add
            # this representation as a new molecule
            else:
                smiles2mol[smiles] = mol                
                
    os.makedirs(output_directory)
    
    # Assign names and write out files
    for unique_mol_index, smiles in enumerate(smiles2mol.keys()):
        mol_name = f'{group_name}-{unique_mol_index:05d}'
        smiles2mol[smiles].properties['group_name'] = group_name
        smiles2mol[smiles].properties['group_index'] = unique_mol_index
        smiles2mol[smiles].name = mol_name
        mol_copy = Molecule(smiles2mol[smiles])
        # Pop off now-nonessential metadata
        mol_copy.properties.pop('original_file')
        mol_copy.properties.pop('original_file_index')
        # Write the graph representation as a fully mapped SMILES
        with open(os.path.join(output_directory, f'{mol_name}.smi'), 'w') as of:
            cmiles = smiles2mol[smiles].to_smiles(mapped=True)
            of.write(cmiles)

        # Write conformers
        # Don't try to write conformers if there aren't any
        if mol_copy.conformers is None:
            continue
        
        for conf_index, conformer in enumerate(smiles2mol[smiles].conformers):
            mol_copy._conformers = None
            #raise Exception(conformer)
            mol_copy.add_conformer(conformer)
            mol_copy.properties['conformer_index'] = conf_index
            out_file_name = f'{mol_copy.name}-{conf_index:02d}.sdf'
            mol_copy.to_file(os.path.join(output_directory, out_file_name), file_format='sdf')
            
    
if __name__ == '__main__':
    validate_and_assign()
