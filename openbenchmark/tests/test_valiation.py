import pytest
import qcportal as ptl

from openbenchmark.validation import QCValidator


class TestValidator:

    @pytest.fixture
    def dataset(self):
        client = ptl.FractalClient()
        dataset = client.get_collection('OptimizationDataset', 'OpenFF Full Optimization Benchmark 1')
        return dataset

    def test_validator(self, dataset):
        QCValidator(dataset=dataset)
