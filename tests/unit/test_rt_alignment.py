""" unit testing of rt_alignment functions """
# pylint: disable=missing-function-docstring

import pandas as pd

from metatlas.targeted import rt_alignment


def test_get_rts01(metatlas_dataset):
    # pylint: disable=line-too-long
    rts_df = rt_alignment.get_rts(metatlas_dataset, include_atlas_rt_peak=False)
    assert f"{rts_df.iloc[0]['min']:0.5f}" == "2.29224"
    assert (
        rts_df.to_json()
        == """{"20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":2.2922415733},"mean":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":2.2922415733},"median":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":2.2922415733},"min":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":2.2922415733},"max":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":2.2922415733},"standard deviation":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":null},"standard error":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":null},"#NaNs":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":0}}"""
    )


def test_plot_actual_vs_aligned_rts01(model):
    arrays = [[]]
    rts_df = pd.DataFrame(data={"1": [], "2": [], "3": [], "4": [], "5": [], "6": []})
    rt_alignment.plot_actual_vs_aligned_rts(arrays, arrays, rts_df, "file_name", model, model)


def test_align_atlas(atlas_with_2_cids, model):
    out = rt_alignment.align_atlas(atlas_with_2_cids, model, 0)
    assert out.name == "HILICz150_ANT20190824_PRD_EMA_Unlab_POS"
    assert atlas_with_2_cids.compound_identifications[0].rt_references[0].rt_peak != -0.8035359946292822
    assert out.compound_identifications[0].rt_references[0].rt_peak == -0.8035359946292822
    assert atlas_with_2_cids.compound_identifications[1].rt_references[0].rt_peak != 0.02331840799266649
    assert out.compound_identifications[1].rt_references[0].rt_peak == 0.02331840799266649
    assert len(out.compound_identifications) == 2
