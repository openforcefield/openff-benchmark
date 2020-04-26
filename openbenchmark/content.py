import tarfile
import json
import logging
from functools import lru_cache

import pandas as pd

from .utils import datapath

logger = logging.getLogger(__name__)

DATASETS = {
    "trim3": "TRIM3 (sample)",
    # "OpenFF-FullOpt1-sample100": "OpenFF Full Optimization Benchmark 1 (sample)",
}
FORCEFIELDS = {
    "smirnoff": "Smirnoff",
    "gaff": "GAFF",
    "gaff2": "GAFF2",
    "mmff94": "MMFF94",
    "mmff94s": "MMFF94s",
    # "opls3e": "OPLS3e",
    "parsley": "Parsley 1.0",
    "parsley110": "Parsley 1.1",
}


def available_datasets():
    return [{"label": label, "value": key} for (key, label) in DATASETS.items()]


def available_forcefields():
    return [{"label": label, "value": key} for (key, label) in FORCEFIELDS.items()]


@lru_cache(maxsize=10)
def get_dataframe_from_trim3(dataset, forcefields):
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
            df = pd.read_csv(csv)
            yield df


def get_json_from_off(dataset, forcefields=("dataset",)):
    if dataset not in DATASETS:
        raise ValueError(f"Dataset {dataset} must be one of {list(DATASETS.keys())}")

    tar_filename = datapath(f"{dataset}.tar.gz")
    dicts = []
    forcefields = list(forcefields) + ["dataset"]
    with tarfile.open(tar_filename, mode="r:gz") as tar:
        for forcefield in forcefields:
            try:
                f = tar.extractfile(f"{forcefield}.json")
                d = json.load(f)
                d["forcefield"] = forcefield
                dicts.append(d)
            except KeyError:
                logger.warning(
                    "Dataset %s does not contain info for forcefield %s", dataset, forcefield
                )
    return dicts
