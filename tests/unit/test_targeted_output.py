""" tests for metatlas.targeted.process"""
# pylint: disable=missing-function-docstring,protected-access,unused-argument,too-many-arguments

import glob

from metatlas.io.targeted_output import generate_all_outputs


def test_generate_all_outputs01(metatlas_dataset, hits, mocker, configuration, workflow, analysis):
    mocker.patch("metatlas.plots.dill2plots.get_msms_hits", return_value=hits)
    generate_all_outputs(metatlas_dataset, configuration, workflow, analysis)
    assert len(glob.glob(metatlas_dataset.ids.output_dir + "/*")) == 15
    assert len(glob.glob(metatlas_dataset.ids.output_dir + "/*/*")) == 19
