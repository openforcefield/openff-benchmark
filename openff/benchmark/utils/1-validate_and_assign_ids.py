import click
from openforcefield.topology import Molecule
from openforcefield.utils.toolkits import RDKitToolkitWrapper, AmberToolsToolkitWrapper, ToolkitRegistry


@click.command()
@click.option('-m',
              '--molecules',
              help='SDF file to read containing input molecules.')
@click.option('--groupname',
              help='Group name for assigning IDs to the molecules.')
def validate_and_assign(molecules, groupname):
    """
    Load a molecule dataset from SDF, validate it for common 
    issues, and assign it unique identifiers.
    """

    tkr = ToolkitRegistry([RDKitToolkitWrapper, AmberToolsToolkitWrapper])
    mols = Molecule.from_file(molecules)
    
    # If there's only one  molecule in the file, turn it into a list
    if not isinstance(mols, list):
        mols = [mols]
    
    # Deduplicate multiple instances of the same molecule
    unique_mols = []
    for idx, mol in enumerate(mols):
        this_smiles = mol.to_smiles(toolkit_registry=tkr)
        if this_smiles in unique_mols:
            other_idx = unique_mols.index(this_smiles)
            other_mol = mols[other_idx]
            msg = 'Duplicate molecule detected:\n'
            msg += f'Molecule {other_idx}: {other_mol}\n'
            msg += f'Molecule {idx}: {mol}'
            raise Exception(msg)
        unique_mols.append(this_smiles)
    
        
    # Assign unique identifiers
    for idx, mol in enumerate(mols):
        mol.name = f'{groupname}-{idx:05d}'
        mol.to_file(mol.name + '-001.sdf', file_format='sdf')
        
        

        
if __name__ == '__main__':
    validate_and_assign()
