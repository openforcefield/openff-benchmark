import os

def mols_from_paths(input_paths, sourcefile_keys=False, recursive=False):
    from openforcefield.topology import Molecule
    from openforcefield.utils import toolkits
    
    # make sure we deregister OpenEye, if it is present
    try:
        toolkits.GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(toolkits.OpenEyeToolkitWrapper)
    except toolkits.ToolkitUnavailableException:
        pass

    def read_file(path, mols, sourcefile_keys):
        if (os.path.isfile(path) and path.split('.')[-1].lower() == 'sdf'):
            if sourcefile_keys:
                mols[path] = Molecule.from_file(path, file_format='sdf', allow_undefined_stereo=True)
            else:
                mols.append(Molecule.from_file(path, file_format='sdf', allow_undefined_stereo=True))
        return mols

    if sourcefile_keys:
        mols = {}
    else:
        mols = []

    for path in input_paths:
        if os.path.isdir(path):
            if recursive:
                for root, dirs, files in os.walk(path):
                    for file in files:
                        filepath = os.path.join(root, file)
                        read_file(filepath, mols, sourcefile_keys)
            else:
                for file in os.listdir(path):
                    filepath = os.path.join(path, file)
                    read_file(filepath, mols, sourcefile_keys)
        else:
            read_file(path, mols, sourcefile_keys)

    return mols
