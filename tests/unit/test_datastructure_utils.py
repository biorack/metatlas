""" tests for metatlas.datastructures.utils"""
# pylint: disable=missing-function-docstring,protected-access,unused-argument,too-many-arguments

import os
import pytest

from metatlas.datastructures import metatlas_objects as metob
from metatlas.datastructures.analysis_identifiers import AnalysisIdentifiers
from metatlas.datastructures.metatlas_dataset import MetatlasDataset
from metatlas.datastructures.utils import get_atlas, set_atlas_mz_tolerance
from metatlas.tools import config


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


def test_set_atlas_mz_tolerance01(atlas):
    original_mz_tolerance = atlas.compound_identifications[0].mz_references[0].mz_tolerance
    new_tol = 100
    set_atlas_mz_tolerance(atlas, new_tol, override=True)
    assert original_mz_tolerance != new_tol
    assert atlas.compound_identifications[0].mz_references[0].mz_tolerance == new_tol


def test_set_atlas_mz_tolerance02(atlas):
    original_mz_tolerance = atlas.compound_identifications[0].mz_references[0].mz_tolerance
    new_tol = 100
    set_atlas_mz_tolerance(atlas, new_tol, override=False)
    assert original_mz_tolerance != new_tol
    assert atlas.compound_identifications[0].mz_references[0].mz_tolerance == original_mz_tolerance


def test_set_atlas_mz_tolerance03(atlas):
    atlas.compound_identifications[0].mz_references[0].mz_tolerance = None
    new_tol = 100
    set_atlas_mz_tolerance(atlas, new_tol, override=False)
    assert atlas.compound_identifications[0].mz_references[0].mz_tolerance == new_tol


def test_set_atlas_mz_tolerance04(sqlite_with_atlas_with_2_cids, analysis_ids, analysis_parameters):
    experiment = "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583"
    new_tol = 100
    analysis_parameters["mz_tolerance_override"] = new_tol
    configuration, _, _ = config.get_config(analysis_parameters)
    ids = AnalysisIdentifiers(os.getcwd(), experiment, configuration, "Test-HILIC", "EMA-POS")
    metatlas_dataset = MetatlasDataset(ids=ids, save_metadata=False)
    assert metatlas_dataset.atlas.compound_identifications[0].mz_references[0].mz_tolerance == new_tol


def test_set_atlas_mz_tolerance05(sqlite, atlas_with_2_cids, analysis_ids, analysis_parameters):
    atlas_with_2_cids.compound_identifications[0].mz_references[0].mz_tolerance = None
    metob.store(atlas_with_2_cids)
    experiment = "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583"
    new_tol = 100
    analysis_parameters["mz_tolerance_default"] = new_tol
    analysis_parameters["source_atlas_unique_id"] = atlas_with_2_cids.unique_id
    analysis_parameters["copy_atlas"] = False
    configuration, _, _ = config.get_config(analysis_parameters)
    ids = AnalysisIdentifiers(
        os.getcwd(),
        experiment,
        configuration,
        "Test-HILIC",
        "EMA-POS",
        source_atlas_unique_id=atlas_with_2_cids.unique_id,
    )
    metatlas_dataset = MetatlasDataset(ids=ids, save_metadata=False)
    assert metatlas_dataset.atlas.unique_id == atlas_with_2_cids.unique_id
    assert metatlas_dataset.atlas.compound_identifications[0].mz_references[0].mz_tolerance == new_tol
    assert metatlas_dataset.atlas.compound_identifications[1].mz_references[0].mz_tolerance != new_tol
