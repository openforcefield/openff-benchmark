import numpy as np
import openeye as oechem
from openforcefield.topology.molecule import Molecule
from cmiles.utils import _symbols, BOHR_2_ANGSTROM, ANGSROM_2_BOHR


def run_through_openeye(record):
    coords = record.get_final_molecule().dict()['geometry']
    geometry = np.asarray(coords, dtype=float)*BOHR_2_ANGSTROM
    symbols = record.get_final_molecule().dict()['symbols']

    oemol = oechem.OEMol()

    for s in symbols:
        oemol.NewAtom(_symbols[s])

    oemol.SetCoords(oechem.OEFloatArray(geometry.reshape(-1)))
    oemol.SetDimension(3)

    # Let OpenEye try to determine connectivity instead of trusting the connectivity in the record
    oechem.OEDetermineConnectivity(oemol)
    oechem.OEFindRingAtomsAndBonds(oemol)
    oechem.OEPerceiveBondOrders(oemol)
    oechem.OEAssignImplicitHydrogens(oemol)
    oechem.OEAssignFormalCharges(oemol)
    Molecule.from_openeye(oemol).to_rdkit()


def sanity_check_bond_order(val, ds=None):
    # TODO: get the smiles/cmiles in from the 'extras' in the record
    entry = ds.get_entry(val)
    smiles_in = entry.dict()['attributes']['canonical_isomeric_explicit_hydrogen_smiles']
    mol_in = Molecule.from_smiles(smiles_in)
    mol_in.assign_fractional_bond_orders(bond_order_model='am1-wiberg')
    expected_bond_orders = [b.bond_order for b in mol_in.bonds]

    record = ds.get_record(val, specification='default')

    coords = record.get_final_molecule().dict()['geometry']
    geometry = np.asarray(coords, dtype=float)*BOHR_2_ANGSTROM
    symbols = record.get_final_molecule().dict()['symbols']

    oemol = oechem.OEMol()

    for s in symbols:
        oemol.NewAtom(_symbols[s])

    oemol.SetCoords(oechem.OEFloatArray(geometry.reshape(-1)))
    oemol.SetDimension(3)
    oechem.OEDetermineConnectivity(oemol)
    oechem.OEDetermineConnectivity(oemol)
    oechem.OEFindRingAtomsAndBonds(oemol)
    oechem.OEAssignImplicitHydrogens(oemol)
    oechem.OEPerceiveBondOrders(oemol)

    found_bond_orders = [b.GetOrder() for b in oemol.GetBonds()]

    for x, y in zip(expected_bond_orders, found_bond_orders):
        assert x == y
