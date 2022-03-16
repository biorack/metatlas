"""Atlas minipulation functions"""
from copy import deepcopy


def _get_rt(cid):
    """Extract peak RT time from a compound identification or 0 if not set"""
    try:
        return cid.rt_references[0].rt_peak
    except (AttributeError, IndexError):
        return 0


def sort_atlas(atlas):
    """sort atlas based on rt_peak"""
    new_atlas = deepcopy(atlas)
    new_atlas.compound_identifications = sorted(new_atlas.compound_identifications, key=_get_rt)
    return new_atlas
