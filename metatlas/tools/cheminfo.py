"""cheminformatics related functions"""

import matchms
import numpy as np


def get_parent_mass(precursor_mz: float, adduct: str) -> float:
    """Returns the mass of the input molecule that would result in the supplied precursor_mz and adduct"""
    dummy = matchms.Spectrum(mz=np.array([]),
                             intensities=np.array([]),
                             metadata={"precursor_mz": precursor_mz, "adduct": adduct})
    updated = matchms.filtering.add_parent_mass(dummy)
    return updated.metadata['parent_mass']


def get_precursor_mz(parent_mass: float, adduct: str) -> float:
    """For an input molecule with parent_mass that generates adduct, return the resutling precursor_mz"""
    adducts = matchms.importing.load_adducts_dict()
    if adduct not in adducts:
        raise KeyError("Adduct '%s' is not supported")
    multiplier = adducts[adduct]["mass_multiplier"]
    correction_mass = adducts[adduct]["correction_mass"]
    return (parent_mass + correction_mass) / multiplier
