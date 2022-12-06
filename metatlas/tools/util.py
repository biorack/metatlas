""" stand alone utility functions """

from pathlib import Path
from typing import Optional, TypeVar

import pandas as pd

Generic = TypeVar("Generic")


def or_default(none_or_value: Optional[Generic], default: Generic) -> Generic:
    """
    inputs:
        none_or_value: variable to test
        default: value to return if none_or_value is None
    """
    return none_or_value if none_or_value is not None else default


def repo_path() -> Path:
    """returns Path object pointing to root directory of the git repo"""
    return Path(__file__).resolve().parent.parent.parent


def is_atlas_df_subset(large: pd.DataFrame, small: pd.DataFrame) -> bool:
    """
    Returns True if RT-compound-adduct space in small is a subset of the
    RT-compound-adduct space in large. Otherwise return False

    A subset contains only rows where a matching row occurs in the
    superset where the following columns have equal values in the
    subset row and the superset row:
    ['inchi_key', 'compound_name', 'rt_max', 'detected_polarity', 'mz',
     'mz_tolerance', 'mz_tolerance_units', 'mono_isotopic_molecular_weight',
     'pubchem_compound_id', 'inchi', 'adduct']
    and the RT range [rt_min, rt_max] in the subset row is a subset of the
    RT range in the superset row
    """
    if len(large) == 0 and len(small) > 0:
        return False
    if len(large) >= 0 and len(small) == 0:
        return True
    comp_cols = [
        "inchi_key",
        "compound_name",
        "detected_polarity",
        "mz",
        "mz_tolerance",
        "mz_tolerance_units",
        "mono_isotopic_molecular_weight",
        "pubchem_compound_id",
        "inchi",
        "adduct",
    ]
    range_cols = ["rt_min", "rt_max"]
    small_comp = small[comp_cols + range_cols].drop_duplicates()
    large_comp = large[comp_cols + range_cols].drop_duplicates()
    merged = large_comp.merge(small_comp, on=comp_cols, suffixes=("_L", "_S"))
    if len(small_comp) != len(merged):
        return False
    return pd.Series.all(merged.loc[:, "rt_min_L"] <= merged.loc[:, "rt_min_S"]) and pd.Series.all(
        merged.loc[:, "rt_max_L"] >= merged.loc[:, "rt_max_S"]
    )
