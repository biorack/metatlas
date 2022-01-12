""" Tests of RClone """
# pylint: disable=missing-function-docstring

import os
import subprocess

import pytest

from metatlas.io import rclone


def has_rclone():
    return os.system("rclone") == 256


def rclone_path():
    result = subprocess.run(["which", "rclone"], stdout=subprocess.PIPE, check=True)
    return result.stdout.decode("utf-8").rstrip()


def test_config_file01():
    rci = rclone.RClone("/bin/foobarz")
    assert rci.config_file() is None


@pytest.mark.skipif(not has_rclone(), reason="rclone not in PATH")
def test_config_file02():
    rci = rclone.RClone(rclone_path())
    config_path = os.path.join(os.environ["HOME"], ".config", "rclone", "rclone.conf")
    assert rci.config_file() == config_path


def test_parse_path01():
    assert ("drive", ["foo", "bar"]) == rclone.parse_path("drive:foo/bar")
