import glob
import os
import copy
import re
from simtk import unit
import shutil
import logging
from rdkit import Chem
from simtk import unit
import numpy as np
from tqdm import tqdm
import io

off_logger = logging.getLogger('openforcefield.utils.toolkits')
prev_log_level = off_logger.getEffectiveLevel()
off_logger.setLevel(logging.ERROR)

from openforcefield.topology import Molecule
from openforcefield.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, OpenEyeToolkitWrapper

off_logger.setLevel(prev_log_level)

#logger = logging.getLogger(__name__)


oetk_loaded = False
for tkw in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:
    if isinstance(tkw, OpenEyeToolkitWrapper):
        oetk_loaded = True
if oetk_loaded:
    GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(OpenEyeToolkitWrapper)

RMS_MATRIX_CACHE = {}

def greedy_conf_deduplication(offmol, rms_cutoff, user_confs=None):
    """
    Parameters
    ----------
    offmol : An openforcefield.topology.Molecule object
        A molecule with one or more conformers
    rms_cutoff : float
        The RMSD cutoff (in angstroms) to apply during deduplication
    user_confs : iterable of integer
        The indices of existing molecule conformers which should never be removed.
        
    Returns
    -------
    confs_to_delete : set of int
        The indices of conformers that will be deleted by this cutoff
    
    """
    # Try to pull the rms matrix out of the cache. Calculate it if it doesn't yet exist
    mol_hash = hash(str(offmol.to_dict()))
    #print('a', rms_cutoff)
    if mol_hash not in RMS_MATRIX_CACHE:
        rdmol = offmol.to_rdkit()
        rdmol = Chem.RemoveHs(rdmol)
        rms_matrix = np.zeros((rdmol.GetNumConformers(), rdmol.GetNumConformers()))
        for i in range(rdmol.GetNumConformers()):
            for j in range(i+1, rdmol.GetNumConformers()):
                #print(i,j)
                rmsd = Chem.rdMolAlign.GetBestRMS(rdmol,
                                               rdmol,
                                               prbId=j,
                                               refId=i)
                rms_matrix[i,j] = rmsd
                rms_matrix[j,i] = rmsd
        RMS_MATRIX_CACHE[mol_hash] = rms_matrix
    else:
        rms_matrix = RMS_MATRIX_CACHE[mol_hash]
        
    if user_confs == None:
        user_confs = []
    #print('b')
    # Do greedy deduplication
    i = 0
    confs_to_delete = set()
    while i < offmol.n_conformers:
        if i not in confs_to_delete:
            #rmslist = []
            #Chem.rdMolAlign.AlignMolConformers(rdmol, RMSlist=rmslist)
            #print(rmslist)
            #for idx, rmsd in enumerate(rmslist):
                # +1 because the reference conf is not included in the rmslist
                # so when i is 0 and idx is 0, the first entry in rmslist is actually conf 1
                #compare_idx = i + idx + 1 
                #if (rmsd < rms_cutoff) and (compare_idx not in user_confs):
                #    confs_to_delete.add(compare_idx)
            
            # Use CalcRMS because alignMol doesn't try multiple atom mappings
            for j in range(i+1, offmol.n_conformers):
                rmsd = rms_matrix[i,j]
                if (rmsd < rms_cutoff) and (j not in user_confs):
                    confs_to_delete.add(j)
        #rdmol.RemoveConformer(i)        
        i += 1
    return confs_to_delete


def align_offmol_conformers(offmol):
    from rdkit.Chem import rdMolAlign
    rdmol = offmol.to_rdkit()
    rmslist = []
    rdMolAlign.AlignMolConformers(rdmol, RMSlist=rmslist)
    offmol2 = Molecule.from_rdkit(rdmol)
    # The RDKit roundtrip above may have messed with the properties dict,
    # so transfer all the aligned confs to a copy of the original mol.
    return_mol = copy.deepcopy(offmol)
    return_mol._conformers = []
    for aligned_conf in offmol2.conformers:
        return_mol.add_conformer(aligned_conf)
    return return_mol, rmslist
    

def gen_confs_preserving_orig_confs(conformer_mols, 
                                    target_n_confs=10, 
                                    min_rmsd=2, 
                                    rmsd_step=0.1):
    """
    Parameters
    ----------
    conformer_mols : dict of openforcefield.topology.Molecule
        A dict of single-conformer openFF Molecule objects. This must contain
        at least one item, and all items in this list must be of the same molecule.
        Keys are strings of the form '00', '01', ...
        The atom indexing must be the same in all of these molecules.
        These conformers will be used during RMSD deduplication to prune out nearby
        generated conformers. These conformers are guaranteed to be in the output.
    target_n_confs : int
        The number of conformers desired for each molecule. The RMS cutoff will be 
        selected as the LOWEST RMSD cutoff which produces target_n_confs 
        conformers or less.
    min_rmsd : float
        The minimum RMSD cutoff (in angstroms) to permit between output conformers. 
    rmsd_step : float
        The amount to increate the RMSD cutoff each iteration
        
    Returns
    -------
    final_confs : dict of {str: openforcefield.topology.Molecule}
        A dict of single-conformer openFF Molecule objects. 
        Keys are strings of the form '00', '01', ...
        The keys and values from the conformer_mols input are guaranteed to appear here.
        Note: If the number of conformers in conformer_mols is already 
        equal to or larger than target_n_confs, then conformer_mols
        will be returned unmodified
        
    Raises
    ------
    
    """
    
    logging.info(f"Generating conformers for {conformer_mols['00'].name}")
    
    # If there are already enough conformers, just return immediately
    if len(conformer_mols) >= target_n_confs:
        return conformer_mols
    
    n_user_confs = len(conformer_mols)
    
    # Make a single offmol with all conformers 
    mol = copy.deepcopy(conformer_mols['00'])
    sorted_conf_idxs = sorted([i for i in conformer_mols.keys()])
    for conf_idx in sorted_conf_idxs[1:]:
        mol.add_conformer(conformer_mols[conf_idx].conformers[0])
    
    # Generate a large number of candidate conformers
    mol.generate_conformers(n_conformers=5 * target_n_confs,
                            # Multiply by 1.5 because this RMSD include hydrogens
                            rms_cutoff=min_rmsd * 1.5 * unit.angstrom, 
                            clear_existing=False)
    
    mol, rmslist = align_offmol_conformers(mol)
    
    # Run the "RMSD ratchet", increasing the RMS cutoff until we have a number
    # of conformers equal to or less than the target number
    logging.info(f'initially {mol.n_conformers} conformers, {n_user_confs} are from user')
    rmsd_cutoff = min_rmsd
    confs_to_delete = None
    while True:
        confs_to_delete = greedy_conf_deduplication(mol, 
                                                    rmsd_cutoff, 
                                                    user_confs=range(n_user_confs))
        n_remaining_confs = mol.n_conformers - len(confs_to_delete)
        logging.info(f'rms_cutoff {rmsd_cutoff:0.2f}: {n_remaining_confs} confs remain')
        if n_remaining_confs <= target_n_confs:
            logging.info(f'Target number of conformers reached. Stopping iteration.')
            break
        rmsd_cutoff += rmsd_step
    
    # Unpack the conformers back to OFFMols, skipping those which should 
    # be deleted from the output
    output_mols = {}
    confs_to_write = [i for i in range(mol.n_conformers) if i not in confs_to_delete]
    logging.info(f'conformers to write are {confs_to_write}')
    for output_conf_index, generated_conf_index in enumerate(confs_to_write):
        new_mol = copy.deepcopy(mol)
        new_mol._conformers = []
        new_mol.add_conformer(mol.conformers[generated_conf_index])
        new_mol.name = mol.name
        new_mol.properties['conformer_index'] = output_conf_index
        output_mols[f'{output_conf_index:02d}'] = new_mol
    
    return output_mols    
        
    

def generate_conformers(group2idx2mols2confs):

    #logging.basicConfig(filename=os.path.join(output_directory,'log.txt'),
    #                    level=logging.DEBUG)

    # Generate new conformers
    print('Generating new conformers')

    # TODO: Make tests for error logic

    success_mols = []
    error_mols = []

    sorted_group_ids = sorted(list(group2idx2mols2confs.keys()))
    for group_id in sorted_group_ids:
        print(f'Processing group {group_id}')
        sorted_mol_ids = sorted(list(group2idx2mols2confs[group_id].keys()))
        for mol_idx in tqdm(sorted_mol_ids):
            try:
                conf_dict = group2idx2mols2confs[group_id][mol_idx]
                output_confs = gen_confs_preserving_orig_confs(conf_dict,
                                                               target_n_confs=10,
                                                               min_rmsd=1.5)
                group2idx2mols2confs[group_id][mol_idx] = output_confs
            except Exception as e:
                mol_id = f'{group_id}-{mol_idx}'
                logging.info(f'Error [{e}] for molecule {mol_id}. Writing to error_mols')
                try:
                    for conf_idx in group2idx2mols2confs[group_id][mol_idx]:
                        conf_file = f'{mol_id}-{int(conf_idx):02d}.sdf'
                        conf_mol = group2idx2mols2confs[group_id][mol_idx][conf_idx]
                        error_mols.append((conf_mol, e))
                        #conf_mol.to_file(os.path.join(error_dir, conf_file), file_format='sdf')
                except Exception as e2:
                    logging.info(f'Unable to write all error structures to file. Encountered error [{e}]')
                    e = str(e) + '\n Then, when trying to write out the conformers:\n' + str(e2)
                    error_mols.append((Molecule(), e))
                #with open(os.path.join(error_dir, f'{mol_id}.txt'), 'w') as of:
                #    of.write(str(e))
                del group2idx2mols2confs[group_id][mol_idx]
                continue

    # Write outputs
    #for group_id in group2idx2mols2confs:
    #    for mol_idx in group2idx2mols2confs[group_id]:
            for conf_idx in group2idx2mols2confs[group_id][mol_idx]:
                mol_name = f'{group_id}-{int(mol_idx):05d}-{int(conf_idx):02d}'
                this_conf = group2idx2mols2confs[group_id][mol_idx][conf_idx]
                this_conf.properties['group_name'] = group_id
                this_conf.properties['molecule_index'] = mol_idx
                this_conf.name = mol_name
                
                this_conf.properties['conformer_index'] = conf_idx
                out_file_name = f'{this_conf.name}.sdf'

                # Ensure this conformer is loadable
                import tempfile
                try:
                    # Do a roundtrip save/load to ensure this won't crash in subsequent steps
                    with tempfile.NamedTemporaryFile(suffix='.sdf') as of:
                        this_conf.to_file(of.name, file_format='sdf')
                        of.seek(0)
                        test_loaded_mol = Molecule.from_file(of.name, file_format='sdf')
                        test_loaded_mol.to_rdkit()
                    success_mols.append(this_conf)
                except Exception as e:
                    error_mols.append((this_conf, e))
                #this_conf.to_file(os.path.join(output_directory, out_file_name), file_format='sdf')

    return success_mols, error_mols

        
