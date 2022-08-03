""" tests for metatlas.datastructures.utils"""
# pylint: disable=missing-function-docstring,protected-access,unused-argument,too-many-arguments

import pytest

from metatlas.datastructures.utils import get_atlas


def test_get_atlas01(sqlite_with_atlas, atlas):
    result = get_atlas(atlas.unique_id)
    assert result.name == atlas.name
    assert len(result.compound_identifications) == 1
    assert result.compound_identifications[0].compound[0].name == "2'-deoxyadenosine"


def test_get_atlas02(sqlite_with_atlas, atlas):
    result = get_atlas(atlas.unique_id, atlas.name)
    assert result.name == atlas.name


def test_get_atlas03(sqlite_with_atlas, atlas):
    with pytest.raises(ValueError):
        get_atlas(atlas.unique_id, "FOOBAR_NOT_ATLAS_NAME")
