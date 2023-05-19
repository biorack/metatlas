""" tests for metatlas.targeted.process"""
# pylint: disable=missing-function-docstring,protected-access,unused-argument


from metatlas.targeted.process import annotation_gui


def test_annotation_gui01(metatlas_dataset, hits, mocker, instructions):
    mocker.patch("metatlas.plots.dill2plots.get_msms_hits", return_value=hits)
    mocker.patch("pandas.read_csv", return_value=instructions)
    agui = annotation_gui(metatlas_dataset)
    agui.compound_idx = 0
    agui.set_msms_flag("1, co-isolated precursor but all reference ions are in sample spectrum")
    agui.set_peak_flag("remove")
    agui.data.set_rt(0, "rt_min", 2.1245)
    agui.data.set_rt(0, "rt_max", 2.4439)
    assert metatlas_dataset.rts[0].rt_min == 2.1245
    assert metatlas_dataset.rts[0].rt_max == 2.4439
    assert metatlas_dataset.data[0][0]["identification"].ms1_notes == "remove"
    assert (
        metatlas_dataset.data[0][0]["identification"].ms2_notes
        == "1, co-isolated precursor but all reference ions are in sample spectrum"
    )


def test_annotation_gui02(metatlas_dataset, hits, mocker, instructions):
    metatlas_dataset[0][0]["identification"].compound = []
    mocker.patch("metatlas.plots.dill2plots.get_msms_hits", return_value=hits)
    mocker.patch("pandas.read_csv", return_value=instructions)
    annotation_gui(metatlas_dataset)
