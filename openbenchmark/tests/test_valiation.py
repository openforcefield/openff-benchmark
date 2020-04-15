import pytest

import qcportal as ptl
from openforcefield.utils.toolkits import UndefinedStereochemistryError

from openbenchmark.validation import run_through_openeye


client = ptl.FractalClient()

ds = client.get_collection('OptimizationDataset', 'OpenFF Full Optimization Benchmark 1')


def test_run_through_openeye_good():
    good_smiles_in = 'CC(C)O-0'
    run_through_openeye(ds.get_record(good_smiles_in, specification='default'))

def test_run_through_openeye_bad():
    # TODO: Split this out into various tests capturing specific
    bad_smiles_in = 'CO/N=C/1\C[N@](C[C@H]1C[NH3+])c2c(cc3c(=O)c(cn(c3n2)C4CC4)C(=O)[O-])F-3'
    with pytest.raises(UndefinedStereochemistryError):
        run_through_openeye(ds.get_record(bad_smiles_in, specification='default'))
