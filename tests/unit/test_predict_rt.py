""" unit testing of predict_rt functions """
# pylint: disable=missing-function-docstring

from metatlas.tools import predict_rt


def test_get_rts01(metatlas_dataset):
    # pylint: disable=line-too-long
    rts_df = predict_rt.get_rts(metatlas_dataset, include_atlas_rt_peak=False)
    assert f"{rts_df.iloc[0]['min']:0.5f}" == "2.29224"
    assert (
        rts_df.to_json()
        == """{"20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":2.2922415733},"mean":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":2.2922415733},"median":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":2.2922415733},"min":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":2.2922415733},"max":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":2.2922415733},"standard deviation":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":null},"standard error":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":null},"#NaNs":{"0000_2deoxyadenosine_positive_M+H252p1091_2p20":0}}"""
    )


def test_adjust_atlas_rt_range01(atlas):
    orig_rt_min = atlas.compound_identifications[0].rt_references[0].rt_min
    mod_atlas = predict_rt.adjust_atlas_rt_range(atlas, -0.1, 0.1)
    mod_rt_min = mod_atlas.compound_identifications[0].rt_references[0].rt_min
    assert orig_rt_min != mod_rt_min


def test_adjust_atlas_rt_range02(atlas):
    orig_rt_max = atlas.compound_identifications[0].rt_references[0].rt_max
    mod_atlas = predict_rt.adjust_atlas_rt_range(atlas, -0.1, 0.1)
    mod_rt_max = mod_atlas.compound_identifications[0].rt_references[0].rt_max
    assert orig_rt_max != mod_rt_max
