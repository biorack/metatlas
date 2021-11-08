"""Test of EIC plotting functions"""
# pylint: disable=missing-function-docstring

import os

import numpy as np

from metatlas.plots import compound_eic


def test_save_compound_eic_pdf(metatlas_dataset):
    file_name = "foo.pdf"
    compound_eic.save_compound_eic_pdf(metatlas_dataset, 0, file_name, sharey=True, overwrite=True)
    assert os.path.isfile(file_name)


def test_insert_in_sorted_array():
    arr = [1.0, 2.0, 3.0]
    assert np.array_equal(compound_eic.insert_in_sorted_array(arr, [0, 99]), [0, 1, 2, 3, 99])
    assert np.array_equal(compound_eic.insert_in_sorted_array(arr, [1.5]), [1, 1.5, 2, 3])
    assert np.array_equal(compound_eic.insert_in_sorted_array(arr, []), arr)
    assert np.array_equal(compound_eic.insert_in_sorted_array(arr, [1]), [1, 1, 2, 3])


def test_add_interp_at():
    assert compound_eic.add_interp_at([1.0, 2.0], [3.0, 4.0], [1.5]) == ([1, 1.5, 2], [3, 3.5, 4])
