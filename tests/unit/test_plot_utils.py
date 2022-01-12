"""Test of plotting utility functions"""
# pylint: disable=missing-function-docstring

from metatlas.plots import utils


def test_is_in_range():
    assert utils.is_in_range([0, 1, 2], 1, 1) == [False, True, False]
    assert utils.is_in_range([], 1, 1) == []
    assert utils.is_in_range([5], 1, 1) == [False]


def test_subplot_dimensions():
    assert utils.subplot_dimensions(1) == (1, 1)
    assert utils.subplot_dimensions(3) == (3, 1)
    assert utils.subplot_dimensions(4) == (2, 2)
    assert utils.subplot_dimensions(5) == (2, 3)
    assert utils.subplot_dimensions(8) == (3, 3)


def test_colors_generator():
    gen = utils.colors()
    assert next(gen) == utils.BACKGROUND_COLORS[0]
    assert next(gen) == utils.BACKGROUND_COLORS[1]
