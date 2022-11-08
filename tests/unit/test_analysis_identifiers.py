""" tests for AnalysisIds """
# pylint: disable=missing-function-docstring,unused-argument

import os
import pytest
import traitlets

from metatlas.datastructures.analysis_identifiers import AnalysisIdentifiers


def test_analysis_identifiers01(sqlite_with_test_config_atlases, configuration):
    with pytest.raises(traitlets.traitlets.TraitError, match=r"Database does not contain an atlas.*"):
        AnalysisIdentifiers(
            source_atlas_unique_id="Not_A_Real_Atlas_Unique_ID",
            experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
            analysis_number=1,
            google_folder="0B-ZDcHbPi-aqZzE5V3hOZFc0dms",
            project_directory="/foo/bar",
            configuration=configuration,
            workflow="Test-HILIC",
            analysis="EMA-POS",
        )


def test_analysis_identifiers04(sqlite_with_test_config_atlases, configuration):
    with pytest.raises(
        traitlets.traitlets.TraitError,
        match="The 'analysis_number' trait of an AnalysisIdentifiers instance expected an int, not",
    ):
        AnalysisIdentifiers(
            source_atlas_unique_id=configuration.workflows[0].rt_alignment.atlas.unique_id,
            experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
            analysis_number="this is a string",
            google_folder="0B-ZDcHbPi-aqZzE5V3hOZFc0dms",
            project_directory="/foo/bar",
            configuration=configuration,
            workflow="Test-HILIC",
            analysis="EMA-POS",
        )


def test_analysis_identifiers05(sqlite_with_test_config_atlases, configuration):
    with pytest.raises(
        traitlets.traitlets.TraitError,
        match="The 'analysis_number' trait of an AnalysisIdentifiers instance expected an int, not",
    ):
        AnalysisIdentifiers(
            source_atlas_unique_id=configuration.workflows[0].rt_alignment.atlas.unique_id,
            experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
            analysis_number="1",
            google_folder="0B-ZDcHbPi-aqZzE5V3hOZFc0dms",
            project_directory="/foo/bar",
            configuration=configuration,
            workflow="Test-HILIC",
            analysis="EMA-POS",
        )


def test_analysis_identifiers06(sqlite_with_test_config_atlases, configuration):
    with pytest.raises(traitlets.traitlets.TraitError, match="Parameter analysis_number cannot be negative."):
        AnalysisIdentifiers(
            source_atlas_unique_id=configuration.workflows[0].rt_alignment.atlas.unique_id,
            experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
            analysis_number=-9,
            google_folder="0B-ZDcHbPi-aqZzE5V3hOZFc0dms",
            project_directory="/foo/bar",
            configuration=configuration,
            workflow="Test-HILIC",
            analysis="EMA-POS",
        )


def test_analysis_identifiers07(sqlite_with_test_config_atlases, configuration):
    with pytest.raises(
        traitlets.traitlets.TraitError,
        match="Parameter 'experiment' should contains 9 fields when split on '_', but has",
    ):
        AnalysisIdentifiers(
            source_atlas_unique_id=configuration.workflows[0].rt_alignment.atlas.unique_id,
            experiment="experiment_name_not_valid",
            analysis_number=0,
            google_folder="0B-ZDcHbPi-aqZzE5V3hOZFc0dms",
            project_directory="/foo/bar",
            configuration=configuration,
            workflow="Test-HILIC",
            analysis="EMA-POS",
        )


def test_analysis_identifiers08(sqlite_with_test_config_atlases, caplog, mocker, lcmsrun, configuration):
    mocker.patch("metatlas.plots.dill2plots.get_metatlas_files", return_value=[lcmsrun])
    AnalysisIdentifiers(
        source_atlas_unique_id=configuration.workflows[0].rt_alignment.atlas.unique_id,
        experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_EXTRA-FIELD",
        analysis_number=0,
        google_folder="0B-ZDcHbPi-aqZzE5V3hOZFc0dms",
        project_directory=str(os.getcwd()),
        configuration=configuration,
        workflow="Test-HILIC",
        analysis="EMA-POS",
    )
    assert "Parameter 'experiment' should contains 9 fields when split on '_', but has 10." in caplog.text


def test_analysis_identifiers09(sqlite_with_test_config_atlases, mocker, lcmsrun, configuration):
    mocker.patch("metatlas.plots.dill2plots.get_metatlas_files", return_value=[lcmsrun])
    experiment = "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583"
    analysis_number = 0
    google_folder = "0B-ZDcHbPi-aqZzE5V3hOZFc0dms"
    AnalysisIdentifiers(
        project_directory=str(os.getcwd()),
        experiment=experiment,
        analysis_number=analysis_number,
        google_folder=google_folder,
        source_atlas_unique_id=configuration.workflows[0].rt_alignment.atlas.unique_id,
        configuration=configuration,
        workflow="Test-HILIC",
        analysis="EMA-POS",
    )


def test_analysis_identifiers_atlas01(analysis_ids, username):
    assert (
        analysis_ids.atlas == f"505892_OakGall_final_HILICz150_ANT20190824_TPL_EMA_Unlab_POS_{username}_0_0"
    )


def test_analysis_identifiers_atlas02(analysis_ids, username, sqlite_with_test_config_atlases):
    # call .atlas twice to get cached value
    analysis_ids.atlas  # pylint: disable=pointless-statement
    assert (
        analysis_ids.atlas == f"505892_OakGall_final_HILICz150_ANT20190824_TPL_EMA_Unlab_POS_{username}_0_0"
    )


def test_set_output_state01(analysis, analysis_ids):
    with pytest.raises(AssertionError):
        analysis_ids.set_output_state(analysis.parameters, "NOT VALID STATE")


def test_set_output_state02(analysis, analysis_ids):
    analysis_ids.set_output_state(analysis.parameters, "ids_spreadsheet")
    assert analysis_ids.exclude_lcmsruns == ["QC"]


def test_exclude_lcmsruns01(analysis, analysis_ids):
    analysis_ids.exclude_lcmsruns = ["Run34"]
    runs = analysis_ids.lcmsruns
    assert len(runs) == 0


def test_include_lcmsruns01(analysis, analysis_ids):
    analysis_ids.include_lcmsruns = ["Run34"]
    runs = analysis_ids.lcmsruns
    assert len(runs) == 1
