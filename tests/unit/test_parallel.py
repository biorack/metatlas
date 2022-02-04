""" tests of parallel processing utilities """
# pylint: disable=missing-function-docstring
import pytest

from metatlas.tools import parallel


def times_two(num):
    return num * 2


def error_on_zero(num):
    if num == 0:
        raise ValueError("Zero is illegal here!")
    return num


def test_parallel_process01():
    data = [0, 1, 2]
    assert sorted(parallel.parallel_process(times_two, data, len(data))) == [0, 2, 4]


def test_parallel_process02():
    data = [0, 1, 2]
    with pytest.raises(ValueError):
        parallel.parallel_process(error_on_zero, data, len(data))
