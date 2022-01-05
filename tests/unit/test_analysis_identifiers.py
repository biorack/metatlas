""" tests for AnalysisIds """
# pylint: disable=missing-function-docstring,unused-argument

import os
import pytest
import traitlets

from metatlas.datastructures.analysis_identifiers import AnalysisIdentifiers


def test_analysis_identifiers01(sqlite):
    with pytest.raises(traitlets.traitlets.TraitError, match=r"Database does not contain an atlas.*"):
        AnalysisIdentifiers(
            source_atlas="Not_A_Real_Atlas_Name",
            experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
            output_type="FinalEMA-HILIC",
            polarity="positive",
            analysis_number=1,
            google_folder="0B-ZDcHbPi-aqZzE5V3hOZFc0dms",
            project_directory="/foo/bar",
        )


def test_analysis_identifiers02(sqlite_with_atlas, username):
    with pytest.raises(
        traitlets.traitlets.TraitError,
        match="Parameter output_type must be one of ISTDsEtc, FinalEMA-HILIC, data_QC",
    ):
        AnalysisIdentifiers(
            source_atlas=f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}0",
            experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
            output_type="output_type_not_valid",
            polarity="positive",
            analysis_number=1,
            google_folder="0B-ZDcHbPi-aqZzE5V3hOZFc0dms",
            project_directory="/foo/bar",
        )


def test_analysis_identifiers03(username, sqlite_with_atlas):
    with pytest.raises(
        traitlets.traitlets.TraitError,
        match="Parameter polarity must be one of positive, negative, fast-polarity-switching",
    ):
        AnalysisIdentifiers(
            source_atlas=f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}0",
            experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
            output_type="FinalEMA-HILIC",
            polarity="not a polarity value",
            analysis_number=1,
            google_folder="0B-ZDcHbPi-aqZzE5V3hOZFc0dms",
            project_directory="/foo/bar",
        )


def test_analysis_identifiers04(username, sqlite_with_atlas):
    with pytest.raises(
        traitlets.traitlets.TraitError,
        match="The 'analysis_number' trait of an AnalysisIdentifiers instance expected an int, not",
    ):
        AnalysisIdentifiers(
            source_atlas=f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}0",
            experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
            output_type="FinalEMA-HILIC",
            polarity="positive",
            analysis_number="this is a string",
            google_folder="0B-ZDcHbPi-aqZzE5V3hOZFc0dms",
            project_directory="/foo/bar",
        )


def test_analysis_identifiers05(username, sqlite_with_atlas):
    with pytest.raises(
        traitlets.traitlets.TraitError,
        match="The 'analysis_number' trait of an AnalysisIdentifiers instance expected an int, not",
    ):
        AnalysisIdentifiers(
            source_atlas=f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}0",
            experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
            output_type="FinalEMA-HILIC",
            polarity="positive",
            analysis_number="1",
            google_folder="0B-ZDcHbPi-aqZzE5V3hOZFc0dms",
            project_directory="/foo/bar",
        )


def test_analysis_identifiers06(username, sqlite_with_atlas):
    with pytest.raises(traitlets.traitlets.TraitError, match="Parameter analysis_number cannot be negative."):
        AnalysisIdentifiers(
            source_atlas=f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}0",
            experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
            output_type="FinalEMA-HILIC",
            polarity="positive",
            analysis_number=-9,
            google_folder="0B-ZDcHbPi-aqZzE5V3hOZFc0dms",
            project_directory="/foo/bar",
        )


def test_analysis_identifiers07(username, sqlite_with_atlas):
    with pytest.raises(
        traitlets.traitlets.TraitError,
        match="Parameter 'experiment' should contains 9 fields when split on '_', but has",
    ):
        AnalysisIdentifiers(
            source_atlas=f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}0",
            experiment="experiment_name_not_valid",
            output_type="FinalEMA-HILIC",
            polarity="positive",
            analysis_number=0,
            google_folder="0B-ZDcHbPi-aqZzE5V3hOZFc0dms",
            project_directory="/foo/bar",
        )


def test_analysis_identifiers08(username, sqlite_with_atlas, caplog, mocker, lcmsrun):
    mocker.patch("metatlas.plots.dill2plots.get_metatlas_files", return_value=[lcmsrun])
    AnalysisIdentifiers(
        source_atlas=f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}0",
        experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_EXTRA-FIELD",
        output_type="FinalEMA-HILIC",
        polarity="positive",
        analysis_number=0,
        google_folder="0B-ZDcHbPi-aqZzE5V3hOZFc0dms",
        project_directory=str(os.getcwd()),
    )
    assert "Parameter 'experiment' should contains 9 fields when split on '_', but has 10." in caplog.text


def test_analysis_identifiers09(sqlite_with_atlas, username, mocker, lcmsrun, groups_controlled_vocab):
    mocker.patch("metatlas.plots.dill2plots.get_metatlas_files", return_value=[lcmsrun])
    experiment = "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583"
    output_type = "FinalEMA-HILIC"
    polarity = "positive"
    analysis_number = 0
    google_folder = "0B-ZDcHbPi-aqZzE5V3hOZFc0dms"
    AnalysisIdentifiers(
        str(os.getcwd()),
        experiment,
        output_type,
        polarity,
        analysis_number,
        google_folder,
        source_atlas=f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}0",
        groups_controlled_vocab=groups_controlled_vocab,
    )


def test_analysis_identifiers_atlas01(analysis_ids, username):
    assert analysis_ids.atlas == f"505892_OakGall_final_FinalEMA-HILIC_POS_{username}0"


def test_analysis_identifiers_atlas02(analysis_ids, username):
    # call .atlas twice to get cached value
    analysis_ids.atlas  # pylint: disable=pointless-statement
    assert analysis_ids.atlas == f"505892_OakGall_final_FinalEMA-HILIC_POS_{username}0"
