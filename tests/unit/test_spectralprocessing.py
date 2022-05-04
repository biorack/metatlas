""" unit testing of spectralprocessing functions """
# pylint: disable=missing-function-docstring

import pytest

import numpy as np

import metatlas.tools.spectralprocessing as sp


def test_find_all_ms_matches01():
    msv_1 = np.array([[1, 2, 3], [4, 5, 6]], np.float64)
    msv_2 = np.array([[1.2, 2.04, 3.09], [4, 5, 6]], np.float64)
    out = sp.find_all_ms_matches(msv_1, msv_2, mz_tolerance=0.1)
    assert np.array_equal(out[0], [[], [1], [2]])
    assert np.array_equal(out[1], [1.0, 2.0])


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
