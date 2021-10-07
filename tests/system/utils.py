# pylint: disable=missing-function-docstring, missing-module-docstring, duplicate-code

import os
import subprocess


def num_files_in(path):
    """Returns number of files in path. Does not count directories"""
    return int(subprocess.check_output(f"find {str(path)} -type f | wc -l", shell=True, text=True).strip())


def assert_files_match(expected):
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
                if expected_lines[num] != clean_line:
                    print('Expected line differss from actual:')
                    print('Expected: "{expected_lines[num]}"')
                    print('Actual:   "{expected_lines[num]}"')
                assert expected_lines[num] == clean_line
            if num is None and contents == "":
                continue
            assert len(expected_lines) == num + 1


def exec_docker(image, command, out_path):
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
