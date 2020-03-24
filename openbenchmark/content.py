import tarfile
import pandas as pd

from .utils import datapath

DATASETS = {"trim3": "TRIM3 (sample)"}
FORCEFIELDS = {
    "smirnoff": "Smirnoff",
    "gaff": "GAFF",
    "gaff2": "GAFF2",
    "mmff94": "MMFF94",
    "mmff94s": "MMFF94s",
    "opls3e": "OPLS3e",
    "parsley": "Parsley 1.0",
    "parsley110": "Parsley 1.1",
}


def available_datasets():
    return [{"label": label, "value": key} for (key, label) in DATASETS.items()]


def available_forcefields():
    return [{"label": label, "value": key} for (key, label) in FORCEFIELDS.items()]


def get_dataframe(dataset, forcefields):
    all_datasets = [ds["value"] for ds in available_datasets()]
    if dataset not in all_datasets:
        raise ValueError(f"Dataset {dataset} must be one of {all_datasets}")
    all_forcefields = [ff["value"] for ff in available_forcefields()]
    if any(ff not in all_forcefields for ff in forcefields):
        raise ValueError(
            f"All forcefields selected ({forcefields}) must be one of {all_forcefields}"
        )

    tar_filename = datapath(f"{dataset}.tar.gz")
    with tarfile.open(tar_filename, mode="r:gz") as tar:
        for ff in forcefields:
            csv_filename = f"refdata_{dataset}_full_{ff}.csv"
            csv = tar.extractfile(csv_filename)
            yield pd.read_csv(csv)
