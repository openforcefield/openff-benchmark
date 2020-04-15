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
