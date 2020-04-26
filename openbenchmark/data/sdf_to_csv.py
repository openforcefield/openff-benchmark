"""
Process SDF data into single CSV

Check each function for the expected formats
"""

import sys
import pandas as pd


def parse_sdf(filename):
    """
    Retrieves all tags from the SDF files and dumps them in a CSV file

    Expected tags are:

    * > <Energy FORCEFIELD-NAME>
    * > <SMILES QCArchive>
    * > <Energy QCArchive>
    * > <RMSD to ../some/path>
    * > <TFD to ../some/path>

    This function will also take the record name, assuming it's preceded by
    a `full_*` card title.
    """
    TAGS = {
        "ff_energy": None,
        "qm_energy": None,
        "smiles": None,
        "RMSD": None,
        "TFD": None,
        "name": None,
    }
    records = [["ff_energy", "qm_energy", "smiles", "RMSD", "TFD", "name"]]
    tags = TAGS.copy()
    with open(filename) as f:
        for line in f:
            if line.startswith("> <"):
                if "SMILES" in line:
                    tags["smiles"] = next(f).strip()
                elif "Energy QCArchive" in line:
                    tags["qm_energy"] = float(next(f).strip())
                elif "Energy" in line:
                    tags["ff_energy"] = float(next(f).strip())
                elif "RMSD" in line:
                    tags["RMSD"] = float(next(f).strip())
                elif "TFD" in line:
                    tags["TFD"] = float(next(f).strip())
            elif line.startswith("full_"):
                tags["name"] = next(f).strip()
            elif "$$$$" in line:  # end of record
                records.append(
                    [
                        tags["ff_energy"],
                        tags["qm_energy"],
                        tags["smiles"],
                        tags["RMSD"],
                        tags["TFD"],
                        tags["name"],
                    ]
                )
                tags = TAGS.copy()
    return pd.DataFrame.from_records(records[1:], columns=records[0])


def ddE_to_csv(path):
    """
    This functions takes a DAT file with the following structure

    * One section per forcefield, which starts with # Relative energies for FF $NAME
    * Rows that contain SMILES and energies, separated by whitespace
    """
    from collections import defaultdict

    smiles_to_energies = defaultdict(dict)
    with open(path) as f:
        for line in f:
            if line.startswith("# Relative energies for FF"):
                ff_name = line.strip().split()[-1]
                smiles, energy = None, None
            elif not line.startswith("#"):
                smiles, energy = line.strip().split()
                energy = float(energy)
                smiles_to_energies[smiles][ff_name] = energy

    return pd.DataFrame.from_dict(smiles_to_energies, orient="index")


if __name__ == "__main__":
    from pathlib import Path

    path = Path(sys.argv[1])
    if sys.argv[1].endswith("ddE.dat"):
        df = ddE_to_csv(path)
    else:
        df = parse_sdf(path)

    df.to_csv(path.stem + ".csv")
