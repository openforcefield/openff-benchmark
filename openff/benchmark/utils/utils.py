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
