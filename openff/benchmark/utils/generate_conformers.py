from openforcefield.topology import Molecule
import glob
import os
import re

def generate_conformers(input_directory, output_directory):
    input_3d_files = glob.glob(os.path.join(input_directory,'*.sdf'))
    input_graph_files = glob.glob(os.path.join(input_directory,'*.smi'))

    group2mols = {}
    for input_graph_file in input_graph_files:
        mols = Molecule.from_file(input_graph_file)
        match = re.findall('([a-zA-Z0-9]{3})-([0-9]{5}).smi',
                           input_graph_file)
        match = re.findall('([a-zA-Z0-9]{3})-([0-9]{5})-([0-9]{2}).sdf',
                           input_graph_file)
        print(match)
        
