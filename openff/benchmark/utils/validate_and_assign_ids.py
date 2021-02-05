from tqdm import tqdm
import os
import logging
import tempfile
import copy
import shutil
from openff.benchmark.utils.generate_conformers import align_offmol_conformers, greedy_conf_deduplication

logger = logging.getLogger('openforcefield.utils.toolkits')
prev_log_level = logger.getEffectiveLevel()
logger.setLevel(logging.ERROR)

from openforcefield.topology import Molecule
from openforcefield.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, OpenEyeToolkitWrapper

logger.setLevel(prev_log_level)

oetk_loaded = False
for tkw in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:
    if isinstance(tkw, OpenEyeToolkitWrapper):
        oetk_loaded = True
if oetk_loaded:
    GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(OpenEyeToolkitWrapper)


def validate_and_assign(loaded_mols,
                        group_name,
                        add,
                        existing_output_mols,
                        name_assignments=None):
    """
    Parameters
    ----------

    """
    if name_assignments is None:
        name_assignments = []

    logging.basicConfig(filename='log.txt',
                        #level=logging.DEBUG
                        )
    #logging.setLevel(logging.INFO)
    #this_logger = logging.getLogger()
    #this_logger.setLevel(logging.INFO)
        
    smiles_to_success_mol = {}
    error_mols = []

    existing_smiles_to_mol = {}
    for mol in existing_output_mols:
        existing_smiles_to_mol[mol.to_smiles()] = mol

    # Handle 3d molecules
    print("Validating input molecules and grouping by connection table")
    logging.info("Validating input molecules and grouping by connection table")
    for mol_index, mol in enumerate(tqdm(loaded_mols)):
        # Simulate a SDF file roundtrip to check for errors such as undefined stereochemistry
        try:
            with tempfile.NamedTemporaryFile(suffix='.sdf') as of:
                mol.to_file(of.name, file_format='sdf')
                of.seek(0)
                test_loaded_mol = Molecule.from_file(of.name, file_format='sdf')
                test_loaded_mol.to_rdkit()
        except Exception as e:
            error_mols.append((f'{mol.properties["original_file"]}:{mol.properties["original_file_index"]}', mol, e))
            continue

        # See whether this graph is already in the existing outputs
        smiles = mol.to_smiles()
        if smiles in existing_smiles_to_mol:
            msg = f'Input molecule graph is already present in output.\n'
            msg += f'{mol.name} from {mol.properties["original_file"]}:{mol.properties["original_file"]} '
            msg += f'has an equivalent connection table to existing output'
            msg += f'{existing_smiles_to_mol[smiles]}'
            logging.debug(msg)
            error_mols.append((f'{mol.properties["original_file"]}:{mol.properties["original_file_index"]}', mol, msg))
            continue

        # If we've reached here, then the molecule is validated

        # Pop off now-nonessential metadata
        allowed_properties = ['original_file',
                              'original_file_index',
                              'original_name']
        keys = list(mol.properties.keys())
        for key in keys:
            if key not in allowed_properties:
                mol.properties.pop(key)
        mol.partial_charges = None

        # If this graph molecule IS already known, add this 3d information as a conformer
        if smiles in smiles_to_success_mol:
            try:
                orig_mol = smiles_to_success_mol[smiles]
                _, atom_map = Molecule.are_isomorphic(mol,
                                                      orig_mol,
                                                      return_atom_map=True,
                                                      formal_charge_matching=False,
                                                      aromatic_matching=False,
                                                      #atom_stereochemistry_matching=False,
                                                      #bond_stereochemistry_matching=False,
                                                      )
                reordered_mol = mol.remap(atom_map)
                # Make a temporary copy of the parent mol for conformer alignment and deduplication
                temp_mol = copy.deepcopy(orig_mol)
                temp_mol.add_conformer(reordered_mol.conformers[0])
                temp_mol, _rmslist = align_offmol_conformers(temp_mol)
                # Don't trust rmslist above for deduplication -- It doesn't take into
                # account multiple atom mappings
                rms_cutoff = 0.2 # angstrom
                confs_to_delete = greedy_conf_deduplication(temp_mol,
                                                            rms_cutoff)
                if len(confs_to_delete) > 0:
                    msg = f'Duplicate molecule conformer input detected.\n'
                    msg += f'{mol.properties["original_file"]}:{mol.properties["original_file_index"]} has an RMSD within {rms_cutoff} A '
                    msg += f'to a conformer of the molecule originally loaded from the following file(s) and indices:'
                    msg += f'{orig_mol.properties["original_file"]}:{orig_mol.properties["original_file_index"]}'
                    logging.debug(msg)
                    temp_mol._conformers = [temp_mol.conformers[-1]]
                    error_mols.append((f'{mol.properties["original_file"]}:{mol.properties["original_file_index"]}',
                                       mol, msg))
                    continue
                temp_mol.properties['original_file'].append(mol.properties['original_file'])
                temp_mol.properties['original_file_index'].append(mol.properties['original_file_index'])
                temp_mol.properties['original_name'].append(mol.properties['original_name'])
                smiles_to_success_mol[smiles] = temp_mol
            # This is a catch-all for any unexpected processing errors that are encountered above
            except Exception as e:
                error_mols.append(
                    (f'{mol.properties["original_file"]}:{mol.properties["original_file_index"]}', mol, e))

        # If this graph molecule ISN'T already known, then add
        # this representation as a new molecule
        else:
            # Change the metadata into lists so that we can record it for each conformer
            mol.properties['original_file'] = [mol.properties['original_file']]
            mol.properties['original_file_index'] = [mol.properties['original_file_index']]
            mol.properties['original_name'] = [mol.properties['original_name']]
            smiles_to_success_mol[smiles] = mol

    # Assign names and write out files
    # Preserve a mapping of input filename/mol index to output name
    success_mols = []
    print("Assigning IDs and preparing molecules for output")
    logging.info("Assigning IDs and preparing molecules for output")
    # Determine the highest previously-existing mol index if the output directory is already populated
    if len(existing_output_mols) == 0:
        output_mol_index_start_value = 0
    else:
        existing_mol_indices = [int(mol.properties['molecule_index']) for mol in existing_output_mols]
        output_mol_index_start_value= max(existing_mol_indices) + 1

    for success_mol_index, smiles in enumerate(tqdm(smiles_to_success_mol.keys())):
        unique_mol_index = success_mol_index + output_mol_index_start_value
        mol_name = f'{group_name}-{unique_mol_index:05d}'
        smiles_to_success_mol[smiles].properties['group_name'] = group_name
        smiles_to_success_mol[smiles].properties['molecule_index'] = unique_mol_index
        smiles_to_success_mol[smiles].name = mol_name
        mol_copy = copy.deepcopy(smiles_to_success_mol[smiles])

        # Write conformers
        for conf_index, conformer in enumerate(smiles_to_success_mol[smiles].conformers):
            mol_copy2 = copy.deepcopy(mol_copy)
            mol_copy2.name = f'{mol_copy.name}-{conf_index:02d}'

            orig_file = smiles_to_success_mol[smiles].properties['original_file'][conf_index]
            orig_file_index = smiles_to_success_mol[smiles].properties['original_file_index'][conf_index]
            orig_name = smiles_to_success_mol[smiles].properties['original_name'][conf_index]
            msg = f'Molecule with name {orig_name} from '
            msg += f'file:position {orig_file}:{orig_file_index}'
            msg += f' has passed validation '
            msg += f'and is being renamed to {mol_copy2.name}.'
            logging.debug(msg)

            name_assignments.append((orig_name, orig_file, orig_file_index, mol_copy2.name))
            mol_copy2._conformers = None
            mol_copy2.add_conformer(conformer)
            mol_copy2.properties['conformer_index'] = conf_index
            # Sanitize last remaining metadata
            mol_copy2.properties.pop('original_file')
            mol_copy2.properties.pop('original_file_index')
            mol_copy2.properties.pop('original_name')
            success_mols.append(mol_copy2)
    return success_mols, error_mols, name_assignments


                  
        
    
if __name__ == '__main__':
    validate_and_assign()
