# pylint: disable=missing-function-docstring, missing-module-docstring, duplicate-code

import os
import subprocess
import unicodedata

from typing import Any, Dict

import pandas as pd

PAPERMILL_ENV = {"PAPERMILL_EXECUTION": "True"}


def num_files_in(path: os.PathLike) -> int:
    """Returns number of files in path. Does not count directories"""
    return int(subprocess.check_output(f"find {str(path)} -type f | wc -l", shell=True, text=True).strip())


def compare_strs(s_1: str, s_2: str) -> bool:
    """String comparision with unicode normalization"""

    def norm_str(in_str: str) -> str:
        """Unicode string normalization"""
        return unicodedata.normalize("NFD", in_str)

    return norm_str(s_1) == norm_str(s_2)


def assert_files_match(expected: Dict[os.PathLike, str]) -> None:
    """
    Throw assertion error if expected does not contain the same data as files on disk
    inputs:
        expected: dict with Path objects as keys and strings representing file contents as values
    returns None
    """
    for path, contents in expected.items():
        with open(path, "r", encoding="utf8") as handle:
            expected_lines = contents.split("\n")
            num = None
            for num, line in enumerate(handle.readlines()):
                assert num < len(expected_lines), "File has more lines than expected."
                clean_line = line.rstrip("\n")
                assert compare_strs(expected_lines[num], clean_line), (
                    "Expected line differss from actual:\n"
                    f'Expected: "{expected_lines[num]}"\n'
                    f'Actual:   "{clean_line}"'
                )
            if num is None and contents == "":
                continue
            assert len(expected_lines) == num + 1, "File has fewer lines than expected."


def assert_dfs_match(expected: Dict[os.PathLike, Dict[str, Dict[int, Any]]]) -> None:
    """
    Throw assertion error if expected does not contain the same data as dataframes on disk
    inputs:
        expected: dict with Path objects as keys and dicts as values
    returns None
    """
    for path, data_dict in expected.items():
        disk_df = pd.read_csv(path)
        expected_df = pd.DataFrame(data_dict)
        pd.testing.assert_frame_equal(disk_df, expected_df)


def exec_docker(image: str, command: str, out_path: os.PathLike, env: Dict) -> None:
    """execute command in image with out_path mounted at /out"""
    env_flags = [x for key, value in env.items() for x in ("--env", f"{key}={value}")]
    subprocess.run(
        [
            "docker",
            "run",
            "--rm",
            f"--user={os.getuid()}:{os.getgid()}",
            "-v",
            f"{os.getcwd()}:/src",
            "-v",
            f"{out_path}:/out",
            *env_flags,
            image,
            "/bin/bash",
            "-c",
            command,
        ],
        check=True,
    )
