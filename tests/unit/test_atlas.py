# pylint: disable=missing-function-docstring, missing-module-docstring, line-too-long

import metatlas.datastructures.atlas as mda


def test_get_rt01(compound_identification):
    assert f"{mda._get_rt(compound_identification):0.3f}" == "2.196"  # pylint: disable=protected-access


def test_sort_atlas(atlas_with_2_cids):
    atlas_with_2_cids.compound_identifications[0].rt_references[0].rt_peak = 999.0
    sorted_atlas = mda.sort_atlas(atlas_with_2_cids)
    assert atlas_with_2_cids.compound_identifications[0].name == sorted_atlas.compound_identifications[1].name
    assert atlas_with_2_cids.compound_identifications[1].name == sorted_atlas.compound_identifications[0].name
