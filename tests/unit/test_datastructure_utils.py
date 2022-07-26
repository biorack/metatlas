""" tests for metatlas.datastructures.utils"""
# pylint: disable=missing-function-docstring,protected-access,unused-argument,too-many-arguments

from metatlas.datastructures.utils import get_atlas


def test_get_atlas01(sqlite_with_atlas):
    query = "HILICz150_ANT20190824_PRD_EMA_Unlab_POS"
    result = get_atlas(query)
    assert result.name == query
    assert len(result.compound_identifications) == 1
    assert result.compound_identifications[0].compound[0].name == "2'-deoxyadenosine"


def test_get_atlas02(sqlite_with_atlas, username):
    query = "HILICz150_ANT20190824_PRD_EMA_Unlab_POS"
    result = get_atlas(name=query, username=username)
    assert result.name == query


def test_get_atlas03(sqlite_with_atlas, atlas):
    query_name = "HILICz150_ANT20190824_PRD_EMA_Unlab_POS"
    result = get_atlas(name=query_name, unique_id=atlas.unique_id)
    assert result.unique_id == atlas.unique_id


def test_get_atlas04(sqlite_with_atlas, atlas):
    result = get_atlas(unique_id=atlas.unique_id)
    assert result.unique_id == atlas.unique_id


def test_get_atlas05(sqlite_with_atlas, username, atlas):
    query_name = "HILICz150_ANT20190824_PRD_EMA_Unlab_POS"
    result = get_atlas(name=query_name, unique_id=atlas.unique_id, username=username)
    assert result.unique_id == atlas.unique_id
