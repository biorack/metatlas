""" unit testing of spectralprocessing functions """
# pylint: disable=missing-function-docstring

import pytest

import numpy as np
import pandas as pd
import metatlas.tools.spectralprocessing as sp


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


## Test MSMS hits sorting methods here

def test_sort_msms_hits_cosine_score():
    unsorted_msms_data_df = pd.DataFrame({
        "name": ["compound1", "compound2", "compound3", "compound4"],
        "score": np.array([0.1, 0.4, 0.3, 0.2]),
        "num_matches": np.array([4, 3, 2, 1]),
        "ref_frags": np.array([7, 6, 5, 4]),
        "data_frags": np.array([10, 9, 8, 7])
    })
    true_sorted_msms_data_df = pd.DataFrame({
        "name": ["compound2", "compound3", "compound4", "compound1"],
        "score": np.array([0.4, 0.3, 0.2, 0.1]),
        "num_matches": np.array([3, 2, 1, 4]),
        "ref_frags": np.array([6, 5, 4, 7]),
        "data_frags": np.array([9, 8, 7, 10])
    })
    true_sorting_method = 'cosine_score'
    function_sorted_msms_data_df, function_sorting_method = sp.sort_msms_hits(unsorted_msms_data_df, sorting_method=true_sorting_method)
    assert np.array_equal(true_sorted_msms_data_df.to_numpy(), function_sorted_msms_data_df.to_numpy())
    assert true_sorting_method == function_sorting_method

def test_sort_msms_hits_sums():
    unsorted_msms_data_df = pd.DataFrame({
        "name": ["compound1", "compound2", "compound3", "compound4"],
        "score": np.array([0.1, 0.4, 0.3, 0.2]),
        "num_matches": np.array([4, 3, 2, 1]),
        "ref_frags": np.array([7, 6, 5, 4]),
        "data_frags": np.array([10, 9, 8, 7])
    })
    true_sorted_msms_data_df = pd.DataFrame({
        "name": ["compound2", "compound1", "compound3", "compound4"],
        "score": np.array([0.4, 0.1, 0.3, 0.2]),
        "num_matches": np.array([3, 4, 2, 1]),
        "ref_frags": np.array([6, 7, 5, 4]),
        "data_frags": np.array([9, 10, 8, 7]),
        "summed_ratios_and_score": np.array([1.2333333333333334, 1.0714285714285714, 0.95, 0.5928571428571429])
    })
    true_sorting_method = 'sums'
    function_sorted_msms_data_df, function_sorting_method = sp.sort_msms_hits(unsorted_msms_data_df, sorting_method=true_sorting_method)
    assert np.array_equal(true_sorted_msms_data_df.to_numpy(), function_sorted_msms_data_df.to_numpy())
    assert true_sorting_method == function_sorting_method

def test_sort_msms_hits_weighted():
    unsorted_msms_data_df = pd.DataFrame({
        "name": ["compound1", "compound2", "compound3", "compound4"],
        "score": np.array([0.1, 0.4, 0.3, 0.2]),
        "num_matches": np.array([4, 3, 2, 1]),
        "ref_frags": np.array([7, 6, 5, 4]),
        "data_frags": np.array([10, 9, 8, 7])
    })
    true_sorted_msms_data_df = pd.DataFrame({
        "name": ["compound2", "compound3", "compound1", "compound4"],
        "score": np.array([0.4, 0.3, 0.1, 0.2]),
        "num_matches": np.array([3, 2, 4, 1]),
        "ref_frags": np.array([6, 5, 7, 4]),
        "data_frags": np.array([9, 8, 10, 7]),
        "weighted_score": np.array([0.4083333333333333, 0.3125, 0.29285714285714287, 0.19821428571428573])
    })
    true_sorting_method = 'weighted'
    function_sorted_msms_data_df, function_sorting_method = sp.sort_msms_hits(unsorted_msms_data_df, sorting_method=true_sorting_method)
    assert np.array_equal(true_sorted_msms_data_df.to_numpy(), function_sorted_msms_data_df.to_numpy())
    assert true_sorting_method == function_sorting_method

def test_sort_msms_hits_numeric_hierarchy():
    unsorted_msms_data_df = pd.DataFrame({
        "name": ["compound1", "compound2", "compound3", "compound4"],
        "score": np.array([0.1, 0.4, 0.3, 0.2]),
        "num_matches": np.array([4, 3, 2, 1]),
        "ref_frags": np.array([7, 6, 5, 4]),
        "data_frags": np.array([10, 9, 8, 7])
    })
    true_sorted_msms_data_df = pd.DataFrame({
        "name": ["compound2", "compound3", "compound4", "compound1"],
        "score": np.array([0.4, 0.3, 0.2, 0.1]),
        "num_matches": np.array([3, 2, 1, 4]),
        "ref_frags": np.array([6, 5, 4, 7]),
        "data_frags": np.array([9, 8, 7, 10]),
        "score_bin": ['0.4-0.5', '0.3-0.4', '0.2-0.3', '0.1-0.2']
    })
    true_sorting_method = 'numeric_hierarchy'
    function_sorted_msms_data_df, function_sorting_method = sp.sort_msms_hits(unsorted_msms_data_df, sorting_method=true_sorting_method)
    assert np.array_equal(true_sorted_msms_data_df.to_numpy(), function_sorted_msms_data_df.to_numpy())
    assert true_sorting_method == function_sorting_method

def test_sort_msms_hits_quantile_hierarchy():
    unsorted_msms_data_df = pd.DataFrame({
        "name": ["compound1", "compound2", "compound3", "compound4"],
        "score": np.array([0.1, 0.4, 0.3, 0.2]),
        "num_matches": np.array([4, 3, 2, 1]),
        "ref_frags": np.array([7, 6, 5, 4]),
        "data_frags": np.array([10, 9, 8, 7])
    })
    true_sorted_msms_data_df = pd.DataFrame({
        "name": ["compound2", "compound3", "compound4", "compound1"],
        "score": np.array([0.4, 0.3, 0.2, 0.1]),
        "num_matches": np.array([3, 2, 1, 4]),
        "ref_frags": np.array([6, 5, 4, 7]),
        "data_frags": np.array([9, 8, 7, 10]),
        "score_bin": [
                        pd.Interval(0.34, 0.4, closed='right'),
                        pd.Interval(0.28, 0.34, closed='right'),
                        pd.Interval(0.16, 0.22, closed='right'),
                        pd.Interval(0.099, 0.16, closed='right')
                    ]
    })
    true_sorting_method = 'quantile_hierarchy'
    function_sorted_msms_data_df, function_sorting_method = sp.sort_msms_hits(unsorted_msms_data_df, sorting_method=true_sorting_method)
    assert np.array_equal(true_sorted_msms_data_df.to_numpy(), function_sorted_msms_data_df.to_numpy())
    assert true_sorting_method == function_sorting_method
