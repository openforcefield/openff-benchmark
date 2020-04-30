import numpy as np
import pandas as pd
import mdtraj as md
from openforcefield.topology.molecule import Molecule
import openeye as oechem
import cmiles
from cmiles.utils import _symbols, BOHR_2_ANGSTROM


#  TODO: Cleaner handling of units

class QCValidator:
    """
    Run throughValidate that a QCArchive entry against its CMILES representation in entry
    """
    def __init__(
        self,
        dataset=None,
        check_connectivity_changes_heuristic=True,
        check_connectivity_changes_wbo=True,
        check_intra_h_bonds=True,
        check_proton_transfer=True,
        check_stereochemical_changes=True,
        check_bond_order_changes=True,
        num_iter=1000,
    ):

        # Should these be set-able by the user?
        self._dataset = dataset
        self._check_connectivity_changes_heuristic = check_connectivity_changes_heuristic,
        self._check_connectivity_changes_wbo = check_connectivity_changes_wbo,
        self._check_intra_h_bonds = check_intra_h_bonds
        self._check_proton_transfer = check_proton_transfer
        self._check_stereochemical_changes = check_stereochemical_changes
        self._check_bond_order_changes = check_bond_order_changes
        self._num_iter = num_iter

    def run_validation(self):

        validation = dict()
        if self._check_connectivity_changes_heuristic:
            validation['connectivity_changes'] = dict()
        if self._check_connectivity_changes_wbo:
            validation['wbo_changes'] = dict()

        for key, record in self._dataset.collection.items():
            if self._check_connectivity_changes_wbo:
                wbo_changes = record.detect_connectivity_changes_wbo(0.65)
                # TODO: QCSubmit collapses multiple molecules into one entry and this
                # does not account for diffconthis looks at sever
                check = any(wbo_changes.values())
                validation['wbo_changes'][key] = check
            if self._check_connectivity_changes_heuristic:
                heuristic_changes = record.detect_connectivity_changes_heuristic()
                check = any(heuristic_changes)
                validation['connectivity_changes'][key] = check

        self._validation = validation

    def to_pandas(self):
        return pd.DataFrame(self._validation_data)

    def check_intra_h_bonds(self, key=None):
        record = self._ds.get_record(key, specification='default')
        qc_mol = record.get_final_molecule()
        qcjson_mol = qc_mol.dict(encoding='json')
        oemol = cmiles.utils.load_molecule(qcjson_mol)
        mol = Molecule.from_openeye(oemol)
        mol.to_file(f'{key}.mol2', file_format='mol2')
        trj = md.load(f'{key}.mol2')
        h_bonds = md.baker_hubbard(trj)
        num_h_bonds = h_bonds.shape[0]
        if num_h_bonds > 0:
            return True
        else:
            return False

    def check_stereochemical_changes(self, key=None):
        return False

    def check_proton_transfer(self, key=None):
        return False

    def check_bond_order_changes(self, key=None):
        # TODO: Replace this with a proper WBO-based inference
        record = self._ds.get_record(key, specification='default')
        qcmol = record.get_final_molecule().dict()
        geometry = np.asarray(qcmol['geometry'], dtype=float) * BOHR_2_ANGSTROM
        symbols = qcmol['symbols']

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

        inferred_bond_orders = [b.GetOrder() for b in oemol.GetBonds()]

        smiles_in_entry = self._ds.get_entry(key).dict()['attributes']['canonical_isomeric_smiles']
        mol_in = Molecule.from_smiles(smiles_in_entry)
        expected_bond_orders = [b.bond_order for b in mol_in.bonds]

        for initial, inferred in zip(expected_bond_orders, inferred_bond_orders):
            assert initial == inferred
            return False


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
    oechem.OEFindRingAtomsAndBonds(oemol)
    oechem.OEAssignImplicitHydrogens(oemol)
    oechem.OEPerceiveBondOrders(oemol)

    found_bond_orders = [b.GetOrder() for b in oemol.GetBonds()]

    for x, y in zip(expected_bond_orders, found_bond_orders):
        assert x == y
