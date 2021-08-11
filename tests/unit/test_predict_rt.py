""" unit testing of predict_rt functions """
# pylint: disable=missing-function-docstring

import os

from metatlas.datastructures import metatlas_dataset as mads
from metatlas.tools import predict_rt


def test_get_rts01(metatlas_dataset):
    rts_df = predict_rt.get_rts(metatlas_dataset, include_atlas_rt_peak=False)
    assert f"{rts_df.iloc[0]['min']:0.5f}" == "2.29224"


def test_get_rts02(
    mocker, df_container, analysis_ids, qc_lcmsruns, sqlite_with_atlas, username, groups_controlled_vocab
):
    mocker.patch(
        "metatlas.io.metatlas_get_data_helper_fun.df_container_from_metatlas_file", return_value=df_container
    )
    mocker.patch("metatlas.plots.dill2plots.get_metatlas_files", return_value=qc_lcmsruns)
    ids = mads.AnalysisIdentifiers(
        source_atlas=f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}0",
        experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
        output_type="FinalEMA-HILIC",
        polarity="positive",
        analysis_number=0,
        project_directory=str(os.getcwd()),
        groups_controlled_vocab=groups_controlled_vocab,
    )
    metatlas_dataset = mads.MetatlasDataset(ids=ids, save_metadata=False)
    rts_df = predict_rt.get_rts(metatlas_dataset, include_atlas_rt_peak=False)
    assert (
        rts_df.to_json()
        == """{"20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_FPS_MS1_0_QC_Pre_Rg70to1050-CE102040--QC_Run6.h5":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":2.2922415733},"20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_0_QC_Pre_Rg70to1050-CE102040--QC_Run7.h5":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":2.2922415733},"20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_FPS_MS1_0_QC_Post_Rg70to1050-CE102040--QC_Run307.h5":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":2.2922415733},"20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_0_QC_Post_Rg70to1050-CE102040--QC_Run308.h5":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":2.2922415733},"mean":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":2.2922415733},"median":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":2.2922415733},"min":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":2.2922415733},"max":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":2.2922415733},"standard deviation":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":0.0},"standard error":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":0.0},"#NaNs":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":0}}"""
    )
