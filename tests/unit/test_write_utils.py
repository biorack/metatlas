""" unit testing of write_utils functions """
# pylint: disable=missing-function-docstring
import os
import pytest

import pandas
from metatlas.io import write_utils


@pytest.fixture(scope="function", autouse=True)
def change_test_dir(request):
    os.chdir(request.fspath.dirname)
    yield
    os.chdir(request.config.invocation_dir)


def test_make_dir_for01(mocker):
    mocker.patch("os.makedirs")
    write_utils.make_dir_for("foo/bar")
    os.makedirs.assert_called_with("foo", exist_ok=True)  # pylint: disable=no-member


def test_make_dir_for02(mocker):
    mocker.patch("os.makedirs")
    write_utils.make_dir_for("bar")
    assert not os.makedirs.called  # pylint: disable=no-member


def test_check_existing_file01(mocker):
    mocker.patch("os.path.exists", return_value=True)
    with pytest.raises(FileExistsError):
        write_utils.check_existing_file("exists_file.txt")


def test_check_existing_file02(mocker):
    mocker.patch("os.path.exists", return_value=False)
    write_utils.check_existing_file("does_not_exist_file.txt")
    # Should not raise an error. No assert needed.


def test_export_dataframe01(mocker):
    mocker.patch("pandas.DataFrame.to_csv")
    mocker.patch("os.path.exists", return_value=False)
    dataframe = pandas.DataFrame({1: [10], 2: [20]})
    write_utils.export_dataframe(dataframe, "foo/bar", "test")
    assert pandas.DataFrame.to_csv.called  # pylint: disable=no-member


def test_raise_on_diff01(mocker):
    mocker.patch("os.path.exists", return_value=False)
    dataframe = pandas.DataFrame({1: [10], 2: [20]})
    write_utils.raise_on_diff(dataframe, "foo/bar", "test")
    # Should not raise an error. No assert needed.


def test_raise_on_diff02(mocker):
    mocker.patch("os.path.exists", return_value=True)
    dataframe = pandas.DataFrame({1: [10], 2: [20]})
    mocker.patch("filecmp.cmp", return_value=True)
    write_utils.raise_on_diff(dataframe, "foo/bar", "test")
    # Should not raise an error. No assert needed.


def test_raise_on_diff03(mocker):
    mocker.patch("os.path.exists", return_value=True)
    mocker.patch("filecmp.cmp", return_value=False)
    to_write = pandas.DataFrame({1: [10], 2: [99]})
    with pytest.raises(ValueError):
        write_utils.raise_on_diff(to_write, "foo/bar", "test")


def test_export_dataframe_die_on_diff01(mocker):
    mocker.patch("os.path.exists", return_value=False)
    dataframe = pandas.DataFrame({1: [10], 2: [20]})
    write_utils.export_dataframe_die_on_diff(dataframe, "foo/bar", "test")
    # Should not raise an error. No assert needed.


def test_export_dataframe_die_on_diff02(mocker):
    mocker.patch("os.path.exists", return_value=True)
    dataframe = pandas.DataFrame({1: [10], 2: [20]})
    mocker.patch("pandas.read_csv", return_value=dataframe)
    write_utils.export_dataframe_die_on_diff(dataframe, "foo/bar", "test")
    # Should not raise an error. No assert needed.


def test_export_dataframe_die_on_diff03(mocker):
    mocker.patch("os.path.exists", return_value=True)
    existing = pandas.DataFrame({1: [10], 2: [20]})
    mocker.patch("pandas.read_csv", return_value=existing)
    to_write = pandas.DataFrame({1: [10], 2: [99]})
    with pytest.raises(ValueError):
        write_utils.export_dataframe_die_on_diff(to_write, "foo/bar", "test")
