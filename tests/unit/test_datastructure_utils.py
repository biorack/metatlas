""" tests for metatlas.datastructures.utils"""
# pylint: disable=missing-function-docstring,protected-access,unused-argument,too-many-arguments

from metatlas.datastructures.utils import get_atlas


def test_get_atlas01(sqlite_with_atlas, username):
    query = "HILICz150_ANT20190824_PRD_EMA_Unlab_POS"
    atlas = get_atlas(query, username)
    assert atlas.name == query
    assert len(atlas.compound_identifications) == 1
    assert atlas.compound_identifications[0].compound[0].name == "2'-deoxyadenosine"
