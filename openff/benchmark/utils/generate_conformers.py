import glob
import os
import re
from simtk import unit
import shutil
import logging
from rdkit import Chem
from simtk import unit
import numpy as np

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
    rdmol = offmol.to_rdkit()
    if user_confs == None:
        user_confs = []
    i = 0
    confs_to_delete = set()
    while i < rdmol.GetNumConformers():
        if i not in confs_to_delete:
            # TODO: Always keep user inputs
            rmslist = []
            Chem.rdMolAlign.AlignMolConformers(rdmol, RMSlist=rmslist)
            #print(rmslist)
            for idx, rmsd in enumerate(rmslist):
                # +1 because the reference conf is not included in the rmslist
                # so when i is 0 and idx is 0, the first entry in rmslist is actually conf 1
                compare_idx = i + idx + 1 
                if (rmsd < rms_cutoff) and (compare_idx not in user_confs):
                    confs_to_delete.add(compare_idx)
        rdmol.RemoveConformer(i)        
        i += 1
    return confs_to_delete


def align_offmol_conformers(offmol):
    from rdkit.Chem import rdMolAlign
    rdmol = offmol.to_rdkit()
    rmslist = []
    rdMolAlign.AlignMolConformers(rdmol, RMSlist=rmslist)
    offmol = Molecule.from_rdkit(rdmol)
    return offmol, rmslist
    

def gen_confs_preserving_orig_confs(conformer_mols, 
                                    target_n_confs=10, 
                                    min_rmsd=1.5, 
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
    mol = Molecule(conformer_mols['00'])
    sorted_conf_idxs = sorted([i for i in conformer_mols.keys()])
    for conf_idx in sorted_conf_idxs[1:]:
        mol.add_conformer(conformer_mols[conf_idx].conformers[0])
    
    # Generate a large number of candidate conformers
    mol.generate_conformers(n_conformers=10 * target_n_confs,
                            rms_cutoff=min_rmsd * unit.angstrom,
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
        new_mol = Molecule(mol)
        new_mol._conformers = []
        new_mol.add_conformer(mol.conformers[generated_conf_index])
        new_mol.name = mol.name
        new_mol.properties['conformer_index'] = output_conf_index
        output_mols[f'{output_conf_index:02d}'] = new_mol
    
    return output_mols    
        
    

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

    logging.basicConfig(filename=os.path.join(output_directory,'log.txt'),
                        #encoding='utf-8',
                        level=logging.DEBUG)

    input_3d_files = glob.glob(os.path.join(input_directory,'*.sdf'))
    input_graph_files = glob.glob(os.path.join(input_directory,'*.smi'))

    group2idx2mols2confs = {}
        

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

    # Generate new conformers
    error_mols = []
    sorted_group_ids = sorted(list(group2idx2mols2confs.keys()))
    for group_id in sorted_group_ids:
        sorted_mol_ids = sorted(list(group2idx2mols2confs[group_id].keys()))
        for mol_idx in sorted_mol_ids:
            try:
                conf_dict = group2idx2mols2confs[group_id][mol_idx]
                output_confs = gen_confs_preserving_orig_confs(conf_dict,
                                                               target_n_confs=10,
                                                               min_rmsd=2)
                group2idx2mols2confs[group_id][mol_idx] = output_confs
            except Exception as e:
                error_mols.append((group_id, mol_idx), e)
                del group2idx2mols2confs[group_id][mol_idx]
            # TODO: Add tqdm progress bar

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

    # TODO: Make tests for error logic
    error_dir = os.path.join(output_directory, 'error_mols')
    os.makedirs(error_dir)
    for (group_id, mol_idx), e in error_mols:
        mol_id = f'{group_id}-{mol_idx}'
        with open(os.path.join(error_dir, f'{mol_id}.txt'), 'w') as of:
            of.write(str(e))

        
