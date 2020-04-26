# Included sample data

The files hereby provided are CSV files exported from SDF files provided by Victoria Lim.

Each dataset is compressed in a `tar.gz` file, which gets uncompressed on the fly when needed.
The TAR contains several CSV files:

* Each `refdata_{dataset}_full_{forcefield}.csv` contains energies and conformations reported by
  each forcefield against the `trim3` QCArchive dataset.
* The `ddE.csv` provides convenient energy minima comparisons of all the forcefields against the QM energies.

Check `sdf_to_csv.py` for details on how the SDF to CSV conversion took place.
