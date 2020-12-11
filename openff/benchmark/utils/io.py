import os

def mols_from_paths(input_paths):
    from openforcefield.topology import Molecule
    from openforcefield.utils import toolkits
    
    # make sure we deregister OpenEye, if it is present
    try:
        toolkits.GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(toolkits.OpenEyeToolkitWrapper)
    except toolkits.ToolkitUnavailableException:
        pass

    mols = []
    for path in input_paths:
        if (os.path.isfile(path) and path.split('.')[-1].lower() == 'sdf'):
            mols.append(Molecule.from_file(path, file_format='sdf', allow_undefined_stereo=True))
        elif os.path.isdir(path):
            for root, dirs, files in os.walk(path):
                for file in files:
                    path = os.path.join(root, file)
                    if (os.path.isfile(path) and path.split('.')[-1].lower() == 'sdf'):
                        mols.append(Molecule.from_file(path, file_format='sdf', allow_undefined_stereo=True))

    return mols
