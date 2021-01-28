import contextlib


def get_data_file_path(relative_path):
    """
    Get the full path to one of the reference files.
    In the source distribution, these files are in ``openff/benchmark/data/``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.

    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the repex folder).
    """

    import os

    from pkg_resources import resource_filename

    fn = resource_filename("openff.benchmark", os.path.join("data", relative_path))

    if not os.path.exists(fn):
        raise ValueError(
            f"Sorry! {fn} does not exist. If you just added it, you'll have to re-install"
        )

    return fn


@contextlib.contextmanager
def temporary_cd(dir_path):
    """Context to temporary change the working directory.
    Parameters
    ----------
    dir_path : str
        The directory path to enter within the context
    Examples
    --------
    >>> dir_path = '/tmp'
    >>> with temporary_cd(dir_path):
    ...     pass  # do something in dir_path
    """
    import os

    prev_dir = os.getcwd()
    os.chdir(os.path.abspath(dir_path))
    try:
        yield
    finally:
        os.chdir(prev_dir)


def prepare_folders(output_directory: str, delete_existing: bool, add: bool) -> str:
    """
    A general folder prep function to help with the preprocess steps, this will handle general file prep and delete_existing and add conflicts.
    Parameters
    ----------
    output_directory: The name of the output directory that should be made.
    delete_existing: If the existing output directory should be overwritten
    add: If the user is trying to add molecules to an existing output folder.

    Returns
    -------
    error_dir: The name of the error directory that has been made/found.
    """
    import os
    import shutil

    error_dir = os.path.join(output_directory, 'error_mols')
    if delete_existing and add:
        raise Exception("Can not specify BOTH --delete-existing AND --add flags")
    # Delete pre-existing output mols if requested
    elif delete_existing and not (add):
        if os.path.exists(output_directory):
            shutil.rmtree(output_directory)
        os.makedirs(output_directory)
        # Create error directory
        error_dir = os.path.join(output_directory, 'error_mols')
        os.makedirs(error_dir)
    elif not (delete_existing) and add:
        if not os.path.exists(output_directory):
            raise Exception(f'--add flag was specified but directory {output_directory} not found')
    # If ADD is FALSE, make new output dir
    elif not (delete_existing) and not (add):
        if os.path.exists(output_directory):
            raise Exception(f'Output directory {output_directory} already exists. '
                            f'Specify `--delete-existing` to remove.')
        os.makedirs(output_directory)
        # Create error directory
        os.makedirs(error_dir)

    return error_dir
