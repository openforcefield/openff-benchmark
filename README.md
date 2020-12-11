OpenFF Benchmark
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.com/openforcefield/openff-benchmark.svg?branch=master)](https://travis-ci.com/openforcefield/openff-benchmark)
[![codecov](https://codecov.io/gh/openforcefield/openff-benchmark/branch/master/graph/badge.svg)](https://codecov.io/gh/openforcefield/openff-benchmark/branch/master)

Comparison benchmarks between public force fields and Open Force Field Initiative force fields

### Input/Output specs

* `validate`
    * Inputs: 
        * group_name: A string, should be 3 characters
        * input-3d-mols: A filename or list of filenames of 3D SDFs, which may contain one or more molecules. MAY have duplicate unique molecules.
        * ~input-graph-mols: A filename or list of filenames of 3D SDFs, which may contain one or more molecules. MUST NOT have duplicate unique molecules (any duplicates are treated as error mols).~
        * output_directory: A directory for the output of this step. Must not exist before this is run.
    * Outputs:
        * ~"graph files": Files of the form `OFF-00000.smi`, containing the graph for each unique input molecule~
        * "3d files": Files of the form `OFF-00000-00.sdf`, which may or may not exist depending on whether input conformers were provided for a molecule. All conformers of the same molecule are guaranteed to have the same indexing, even if they didn't initially.
        * "error_mols/*sdf": SDFs of each molecule that can't be handled by OFF's infrastructure
        * "error_mols/*txt": some output that may help explain why a molecule can't be handled by OFF's infrastructure

* `generate-conformers`
    * Inputs:
        * input_dir: directory path of the output directory from the step above
            * ~This MUST contain `smi` files for each unique molecule, and MAY contain `sdf` files for some molecules~
    * Outputs:
        * output conformers: SDFs of the form `OFF-00000-00.sdf`. 
        * No SMILES are output from this step


### Copyright

Copyright (c) 2020, Open Forcefield Consortium


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
