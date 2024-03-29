""" tests for metatlas.targeted.process"""
# pylint: disable=missing-function-docstring,protected-access,unused-argument,too-many-arguments


from metatlas.io.targeted_output import generate_standard_outputs


def test_generate_standard_outputs01(metatlas_dataset, hits, mocker, workflow, analysis):
    mocker.patch("metatlas.plots.dill2plots.get_msms_hits", return_value=hits)
    generate_standard_outputs(metatlas_dataset, workflow, analysis)
    assert len(list(metatlas_dataset.ids.output_dir.glob("*"))) == 7
    assert len(list(metatlas_dataset.ids.output_dir.glob("*/*"))) == 19
