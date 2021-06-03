# pylint: disable=missing-function-docstring, missing-module-docstring

import pytest

from metatlas.datastructures import metatlas_dataset as mads


@pytest.fixture(name="analysis_ids")
def fixture_analysis_ids(tmp_path):
    return mads.AnalysisIdentifiers(
        "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
        "FinalEMA-HILIC",
        "positive",
        0,
        str(tmp_path),
        atlas="HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_root0",
        username="root",
    )