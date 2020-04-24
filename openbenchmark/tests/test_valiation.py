import pytest

import qcportal as ptl
from openforcefield.utils.toolkits import UndefinedStereochemistryError

from openbenchmark.validation import QCValidator, run_through_openeye, sanity_check_bond_order


client = ptl.FractalClient()

ds = client.get_collection('OptimizationDataset', 'OpenFF Full Optimization Benchmark 1')

good_smiles_in = 'CC(C)O-0'
bad_smiles_in = 'CO/N=C/1\C[N@](C[C@H]1C[NH3+])c2c(cc3c(=O)c(cn(c3n2)C4CC4)C(=O)[O-])F-3'


def test_run_through_openeye():
    run_through_openeye(ds.get_record(good_smiles_in, specification='default'))

    with pytest.raises(UndefinedStereochemistryError):
        run_through_openeye(ds.get_record(bad_smiles_in, specification='default'))


def test_bond_santity_checker():
    sanity_check_bond_order(good_smiles_in, ds=ds)

    with pytest.raises(AssertionError):
        sanity_check_bond_order(bad_smiles_in, ds=ds)


def test_validator_basic():
    qcv = QCValidator(
        ds,
        check_intra_h_bonds=False,
        check_proton_transfer=True,
        check_stereochemical_changes=True,
        check_bond_order_changes=True
    )
    
    qcv.run_validation()
    qcv.to_pandas().head()
