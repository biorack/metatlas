""" unit testing of spectralprocessing functions """
# pylint: disable=missing-function-docstring

import pytest

import numpy as np
import pandas as pd

import metatlas.tools.spectralprocessing as sp


def test_calc_data_to_ref_frag_ratio01():
    num_matches = 5
    spectrum = np.array([[1.1, 1.2, 1.3, 1.4, 1.5], [10, 100, 1000, 100, 10]])
    row = pd.Series({'num_matches': num_matches, 'spectrum': spectrum})
    out = sp.calc_data_to_ref_frag_ratio(row)
    assert out == 1.0

def test_calc_data_to_ref_frag_ratio02():
    num_matches = 1
    spectrum = np.array([[1.1, 1.2, 1.3, 1.4, 1.5], [10, 100, 1000, 100, 10]])
    row = pd.Series({'num_matches': num_matches, 'spectrum': spectrum})
    out = sp.calc_data_to_ref_frag_ratio(row)
    assert out == 0.2

def test_calc_data_to_ref_frag_ratio03():
    num_matches = 0
    spectrum = np.array([[1.1, 1.2, 1.3, 1.4, 1.5], [10, 100, 1000, 100, 10]])
    row = pd.Series({'num_matches': num_matches, 'spectrum': spectrum})
    out = sp.calc_data_to_ref_frag_ratio(row)
    assert out == 0.0

def test_calc_data_to_ref_frag_ratio04():
    num_matches = 0
    spectrum = np.array([])
    row = pd.Series({'num_matches': num_matches, 'spectrum': spectrum})
    out = sp.calc_data_to_ref_frag_ratio(row)
    assert out == 0.0

def test_calc_jaccard_of_spectra01():
    query_spectrum = np.array([[1.1, 1.2, 1.3, 1.4, 1.5], [10, 100, 1000, 100, 10]])
    spectrum = np.array([[1.1, 1.2, 1.3, 1.4, 1.5], [10, 100, 1000, 100, 10]])
    row = pd.Series({'query_spectrum': query_spectrum, 'spectrum': spectrum})
    out = sp.calc_jaccard_of_spectra(row)
    assert out == 1.0

def test_calc_jaccard_of_spectra02():
    query_spectrum = np.array([[1.1], [10]])
    spectrum = np.array([[1.1, 1.2, 1.3, 1.4, 1.5], [10, 100, 1000, 100, 10]])
    row = pd.Series({'query_spectrum': query_spectrum, 'spectrum': spectrum})
    out = sp.calc_jaccard_of_spectra(row)
    assert out == 0.2

def test_calc_jaccard_of_spectra03():
    query_spectrum = np.array([[1.1, 1.3, 1.7, 1.9, 2.1], [10, 100, 1000, 100, 10]])
    spectrum = np.array([[1.1, 1.2, 1.3, 1.4, 1.5], [10, 100, 1000, 100, 10]])
    row = pd.Series({'query_spectrum': query_spectrum, 'spectrum': spectrum})
    out = sp.calc_jaccard_of_spectra(row)
    assert out == 0.25

def test_calc_jaccard_of_spectra04():
    query_spectrum = np.array([[1.1, 1.2, 1.3, 1.4, 1.5], [10, 100, 1000, 100, 10]])
    spectrum = np.array([])
    row = pd.Series({'query_spectrum': query_spectrum, 'spectrum': spectrum})
    out = sp.calc_jaccard_of_spectra(row)
    assert out == 0.0

def test_calc_jaccard_of_spectra05():
    query_spectrum = np.array([])
    spectrum = np.array([[1.1, 1.2, 1.3, 1.4, 1.5], [10, 100, 1000, 100, 10]])
    row = pd.Series({'query_spectrum': query_spectrum, 'spectrum': spectrum})
    out = sp.calc_jaccard_of_spectra(row)
    assert out == 0.0

def array_equal(arr_a, arr_b, equal_nan=False):
    # https://numpy.org/neps/nep-0034-infer-dtype-is-object.html#usage-and-impact
    # required for dealing with "ragged" arrays
    return np.array_equal(np.array(arr_a, dtype=object), np.array(arr_b, dtype=object), equal_nan=equal_nan)


def test_find_all_ms_matches01():
    msv_1 = np.array([[1, 2, 3], [4, 5, 6]], np.float64)
    msv_2 = np.array([[1.2, 2.04, 3.09], [4, 5, 6]], np.float64)
    out = sp.find_all_ms_matches(msv_1, msv_2, mz_tolerance=0.1)
    assert array_equal(out[0], [[], [1], [2]])
    assert array_equal(out[1], [1.0, 2.0])


def test_find_all_ms_matches02():
    msv_1 = np.array([[1], [1]], np.float64)
    msv_2 = np.array([[], []], np.float64)
    out = sp.find_all_ms_matches(msv_1, msv_2, mz_tolerance=0.1)
    assert np.array_equal(out[0], [[]])
    assert np.array_equal(out[1], [])


def test_match_ms_vectors01():
    msv_1 = np.array([[1, 2, 3], [4, 5, 6]], np.float64)
    msv_2 = np.array([[1.2, 2.04, 3.09], [4, 5, 6]], np.float64)
    out = sp.match_ms_vectors(msv_1, msv_2, 0.1, "shape")
    assert np.array_equal(out[0], [np.NaN, 1.0, 2.0], equal_nan=True)
    assert np.array_equal(out[1], [np.NaN, 1.0, 2.0], equal_nan=True)


def test_match_ms_vectors02():
    msv_1 = np.array([[1, 2, 3, 3.1], [4, 5, 6, 4.4]], np.float64)
    msv_2 = np.array([[1.2, 2.04, 3.09], [4, 5, 3]], np.float64)
    out = sp.match_ms_vectors(msv_1, msv_2, 0.1, "intensity")
    assert np.array_equal(out[0], [np.NaN, 1.0, 2.0, np.NaN], equal_nan=True)
    assert np.array_equal(out[1], [np.NaN, 1.0, 2.0], equal_nan=True)


def test_match_ms_vectors03():
    msv_1 = np.array([[1, 2, 3, 3.1], [4, 5, 6, 4.4]], np.float64)
    msv_2 = np.array([[1.2, 2.04, 3.09], [4, 5, 3]], np.float64)
    out = sp.match_ms_vectors(msv_1, msv_2, 0.1, "distance")
    assert np.array_equal(out[0], [np.NaN, 1.0, np.NaN, 2.0], equal_nan=True)
    assert np.array_equal(out[1], [np.NaN, 1.0, 3.0], equal_nan=True)


def test_match_ms_vectors04():
    msv_1 = np.array([[1, 2, 3, 3.1], [4, 5, 6, 4.4]], np.float64)
    msv_2 = np.array([[1.2, 2.04, 3.09], [4, 5, 3]], np.float64)
    with pytest.raises(AssertionError):
        sp.match_ms_vectors(msv_1, msv_2, 0.1, "foobar")


def test_match_ms_vectors05():
    msv_1 = np.array([[1], [1]], np.float64)
    msv_2 = np.array([[], []], np.float64)
    out = sp.match_ms_vectors(msv_1, msv_2, 0.1, "distance")
    assert np.array_equal(out[0], [np.NaN], equal_nan=True)
    assert np.array_equal(out[1], [])
