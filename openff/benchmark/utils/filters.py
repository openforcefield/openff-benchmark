from openff.qcsubmit.workflow_components import SmartsFilter
from typing import TYPE_CHECKING, List, Optional

if TYPE_CHECKING:
    from openforcefield.topology import Molecule
    from openff.qcsubmit.datasets import ComponentResult


def smirks_filter(input_molecules: List["Molecule"], filtered_smirks: List[str], processors: Optional[int] = None) -> "ComponentResult":
    """
    Filter a list of openforcefield.topology.Molecules based on a list of unwanted smirks patterns.

    Parameters
    ----------
    input_molecules: A list of molecules to be processed by the filter.
    filtered_smirks: A list of smirks to be queried against the molecules. Molecules matching any of these smirks will be filtered out.
    processors: The number of processors that should be used when filtering.

    Returns
    -------
        A component result object which contains a list of passed and failed molecules .
    """

    smarts_fil = SmartsFilter(filtered_substructures=filtered_smirks)
    result = smarts_fil.apply(molecules=input_molecules, processors=processors, verbose=True)
    return result
