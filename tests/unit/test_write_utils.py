""" unit testing of write_utils functions """
# pylint: disable=missing-function-docstring

import pytest

from pathlib import Path

import pandas
from metatlas.io import write_utils


def test_make_dir_for01(mocker):
    mocker.patch("pathlib.Path.mkdir")
    write_utils.make_dir_for(Path("foo/bar"))
    Path.mkdir.assert_called_with(parents=True, exist_ok=True)  # pylint: disable=no-member


def test_make_dir_for02(mocker):
    mocker.patch("pathlib.Path.mkdir")
    # doesn't have a parent directory, so no mkdir required
    write_utils.make_dir_for(Path("bar"))
    assert not Path.mkdir.called  # pylint: disable=no-member


def test_check_existing_file01(mocker):
    mocker.patch("pathlib.Path.exists", return_value=True)
    with pytest.raises(FileExistsError):
        write_utils.check_existing_file(Path("exists_file.txt"))


def test_check_existing_file02(mocker):
    mocker.patch("pathlib.Path.exists", return_value=False)
    write_utils.check_existing_file(Path("does_not_exist_file.txt"))
    # Should not raise an error. No assert needed.


def test_export_dataframe01(mocker):
    mocker.patch("pandas.DataFrame.to_csv")
    mocker.patch("pathlib.Path.exists", return_value=False)
    dataframe = pandas.DataFrame({1: [10], 2: [20]})
    write_utils.export_dataframe(dataframe, Path("foo/bar"), "test")
    assert pandas.DataFrame.to_csv.called  # pylint: disable=no-member


def test_raise_on_diff01(mocker):
    mocker.patch("pathlib.Path.exists", return_value=False)
    dataframe = pandas.DataFrame({1: [10], 2: [20]})
    write_utils.raise_on_diff(dataframe, Path("foo/bar"), "test")
    # Should not raise an error. No assert needed.


def test_raise_on_diff02(mocker):
    mocker.patch("pathlib.Path.exists", return_value=True)
    dataframe = pandas.DataFrame({1: [10], 2: [20]})
    mocker.patch("filecmp.cmp", return_value=True)
    write_utils.raise_on_diff(dataframe, Path("foo/bar"), "test")
    # Should not raise an error. No assert needed.


def test_raise_on_diff03(mocker):
    mocker.patch("pathlib.Path.exists", return_value=True)
    mocker.patch("filecmp.cmp", return_value=False)
    to_write = pandas.DataFrame({1: [10], 2: [99]})
    with pytest.raises(ValueError):
        write_utils.raise_on_diff(to_write, Path("foo/bar"), "test")


def test_export_dataframe_die_on_diff01():
    dataframe = pandas.DataFrame({1: [10], 2: [20]})
    write_utils.export_dataframe_die_on_diff(dataframe, Path("foo/bar"), "test")
    # Should not raise an error. No assert needed.


def test_export_dataframe_die_on_diff02():
    dataframe = pandas.DataFrame({1: [10], 2: [20]})
    write_utils.export_dataframe(dataframe, Path("foo/bar"), "test")
    write_utils.export_dataframe_die_on_diff(dataframe, Path("foo/bar"), "test")
    # Should not raise an error. No assert needed.


def test_export_dataframe_die_on_diff03():
    existing = pandas.DataFrame({1: [10], 2: [20]})
    write_utils.export_dataframe(existing, Path("foo/bar"), "test")
    to_write = pandas.DataFrame({1: [10], 2: [99]})
    with pytest.raises(ValueError):
        write_utils.export_dataframe_die_on_diff(to_write, Path("foo/bar"), "test")
