"""Tests of utility functions"""
# pylint: disable=missing-function-docstring,unused-argument,use-implicit-booleaness-not-comparison

import pandas as pd

from metatlas.tools.util import or_default, is_atlas_df_subset


def test_or_default01():
    assert or_default("foo", "bar") == "foo"
    assert or_default(None, "bar") == "bar"
    assert or_default([], [1]) == []
    assert or_default(None, [1]) == [1]


def test_is_atlas_df_subset01(atlas_df):
    assert not is_atlas_df_subset(pd.DataFrame(), atlas_df)


def test_is_atlas_df_subset02(atlas_df):
    mod = atlas_df.copy()
    mod.loc[0, "rt_max"] = 999
    assert is_atlas_df_subset(mod, atlas_df)
    assert not is_atlas_df_subset(atlas_df, mod)


def test_is_atlas_df_subset03(atlas_df):
    mod = atlas_df.copy()
    mod.loc[1, :] = mod.loc[0, :]
    mod.loc[1, "inchi_key"] = "not_real_inchi_key"
    assert is_atlas_df_subset(mod, atlas_df)
    assert not is_atlas_df_subset(atlas_df, mod)
