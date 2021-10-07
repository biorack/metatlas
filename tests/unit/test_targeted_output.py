""" unit testing of targeted_output functions """
# pylint: disable=missing-function-docstring

from metatlas.io import targeted_output


def test_write_msms_fragment_ions01(metatlas_dataset):
    out = targeted_output.write_msms_fragment_ions(metatlas_dataset, min_mz=100, max_mz_offset=0.5)
    assert out.loc[0, "spectrum"] == "[[252.11, 252.16], [100000, 7912]]"
