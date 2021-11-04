"""Test of EIC plotting functions"""
# pylint: disable=missing-function-docstring

import os

import numpy as np

from metatlas.plots import compound_eic


def test_save_compound_eic_pdf(metatlas_dataset):
    file_name = "foo.pdf"
    compound_eic.save_compound_eic_pdf(metatlas_dataset, 0, file_name, sharey=True, overwrite=True)
    assert os.path.isfile(file_name)


def test_is_in_range():
    assert compound_eic.is_in_range([0, 1, 2], 1, 1) == [False, True, False]
    assert compound_eic.is_in_range([], 1, 1) == []
    assert compound_eic.is_in_range([5], 1, 1) == [False]


def test_subplot_dimensions():
    assert compound_eic.subplot_dimensions(1) == (1, 1)
    assert compound_eic.subplot_dimensions(3) == (3, 1)
    assert compound_eic.subplot_dimensions(4) == (2, 2)
    assert compound_eic.subplot_dimensions(5) == (2, 3)
    assert compound_eic.subplot_dimensions(8) == (3, 3)


def test_insert_in_sorted_array():
    arr = [1.0, 2.0, 3.0]
    assert np.array_equal(compound_eic.insert_in_sorted_array(arr, [0, 99]), [0, 1, 2, 3, 99])
    assert np.array_equal(compound_eic.insert_in_sorted_array(arr, [1.5]), [1, 1.5, 2, 3])
    assert np.array_equal(compound_eic.insert_in_sorted_array(arr, []), arr)
    assert np.array_equal(compound_eic.insert_in_sorted_array(arr, [1]), [1, 1, 2, 3])


def test_add_interp_at():
    assert compound_eic.add_interp_at([1.0, 2.0], [3.0, 4.0], [1.5]) == ([1, 1.5, 2], [3, 3.5, 4])


def test_colors_generator():
    gen = compound_eic.colors()
    assert next(gen) == compound_eic.BACKGROUND_COLORS[0]
    assert next(gen) == compound_eic.BACKGROUND_COLORS[1]
