# pylint: disable=missing-function-docstring, missing-module-docstring, duplicate-code

import os
import subprocess
import unicodedata


def num_files_in(path) -> int:
    """Returns number of files in path. Does not count directories"""
    return int(subprocess.check_output(f"find {str(path)} -type f | wc -l", shell=True, text=True).strip())


def compare_strs(s_1: str, s_2: str) -> bool:
    """String comparision with unicode normalization"""
    def norm_str(in_str: str) -> str:
        """Unicode string normalization"""
        return unicodedata.normalize('NFD', in_str)

    return norm_str(s_1) == norm_str(s_2)


def assert_files_match(expected) -> None:
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
                clean_line = line.rstrip("\n")
                if not compare_strs(expected_lines[num], clean_line):
                    print('Expected line differss from actual:')
                    print(f'Expected: "{expected_lines[num]}"')
                    print(f'Actual:   "{expected_lines[num]}"')
                assert compare_strs(expected_lines[num], clean_line)
            if num is None and contents == "":
                continue
            assert len(expected_lines) == num + 1


def exec_docker(image, command, out_path) -> None:
    """execute command in image with out_path mounted at /out"""
    subprocess.run(
        [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{os.getcwd()}:/src",
            "-v",
            f"{out_path}:/out",
            image,
            "/bin/bash",
            "-c",
            command,
        ],
        check=True,
    )
