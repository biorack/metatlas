"""Defines top-level execution of test suite"""
# pylint: disable=missing-function-docstring

import os
from pathlib import Path

import nox

py_versions = ["3.11"]

nox.options.sessions = [
    "flake8",
    "black",
    "pylint-3.11",
    "mypy-3.11",
    "unit_tests-3.11",
    "flake8_nb",
    "black_nb",
    "pylint_nb-3.11",
    "system_tests-3.11",
    "update_git_hooks",
]

# files we can run all the checks on, as they don't contain legacy code that
# has not yet been updated to pass all checks.
more_checks = [
    "metatlas/interfaces/compounds/populate.py",
    "metatlas/io/gdrive.py",
    "metatlas/io/file_converter.py",
    "metatlas/io/rclone.py",
    "metatlas/io/targeted_output.py",
    "metatlas/io/system_utils.py",
    "metatlas/io/write_utils.py",
    "metatlas/io/feature_tools.py",
    "metatlas/datastructures/atlas.py",
    "metatlas/datastructures/analysis_identifiers.py",
    "metatlas/datastructures/id_types.py",
    "metatlas/datastructures/groups.py",
    "metatlas/datastructures/metatlas_dataset.py",
    "metatlas/datastructures/spectrum.py",
    "metatlas/datastructures/utils.py",
    "metatlas/plots/compound_eic.py",
    "metatlas/plots/mirror.py",
    "metatlas/plots/plot_set.py",
    "metatlas/plots/tic.py",
    "metatlas/plots/utils.py",
    "metatlas/scripts/copy_ms_to_new_names.py",
    "metatlas/scripts/yaml_validation.py",
    "metatlas/targeted/rt_alignment.py",
    "metatlas/targeted/process.py",
    "metatlas/tools/add_msms_ref.py",
    "metatlas/tools/cheminfo.py",
    "metatlas/tools/config.py",
    "metatlas/tools/environment.py",
    "metatlas/tools/logging.py",
    "metatlas/tools/notebook.py",
    "metatlas/tools/parallel.py",
    "metatlas/tools/util.py",
    "metatlas/tools/validate_filenames.py",
    "noxfile.py",
    "utils/mzml_to_h5.py",
    "utils/ts.py",
    "utils/validate_config.py",
    "utils/validate_file_name.py",
    "tests",
]

notebooks = [
    "notebooks/reference/",
]

pytest_deps = [
    "coverage==7.0.5",
    "iniconfig==2.0.0",
    "pandas==2.1.2",
    "pluggy==1.0.0",
    "py==1.11.0",
    "pytest==7.2.1",
    "pytest-cov==3.0.0",
    "pytest-mock==3.10.0",
    "toml==0.10.2",
]

mypy_deps = [
    "mypy==0.991",
    "types-PyYAML",
    "types-requests",
    "types-simplejson",
    "types-six",
    "types-tabulate",
]

pylint_deps = [
    "nox==2022.11.21",
    "pylint==2.15.10",
    "pytest==7.2.1",  # so "import pytest" doesn't get reported
]

nbqa_deps = [
    "nbqa==1.6.1",
    "tokenize-rt==5.0.0",
    "importlib-metadata==4.12.0",
    "astroid==2.13.2",
    "wrapt==1.14.1",
    "lazy_object_proxy==1.9.0",
    "isort==5.11.4",
]

flake8_deps = [
    "flake8==6.0.0",
    "flake8-bugbear==23.1.17",
    "flake8-builtins==2.1.0",
    "flake8-comprehensions==3.10.1",
]

pytest_flags = ["-vv", f"--basetemp={Path.home() / '.pytest_tmp'}"]
nbqa_flags = ["--nbqa-dont-skip-bad-cells"]

nox.options.error_on_external_run = True
REUSE_LARGE_VENV = True
# the NB_LINE_LEN value also appears in pre-commit-config.yaml and must be
# manually updated when this value is changed
NB_LINE_LEN = 140

# parallel testings doesn't work well on NERSC login nodes
if "NERSC_HOST" not in os.environ:
    pytest_deps.append("pytest-xdist[psutil]==3.1.0")
    pytest_flags.extend(["--numprocesses", "auto"])


@nox.session(python=py_versions[0])
def flake8_all(session):
    session.install(*flake8_deps)
    session.run("flake8", "metatlas", "tests")


@nox.session(python=py_versions[0])
def flake8(session):
    session.install(*flake8_deps)
    session.run("flake8", *more_checks)


@nox.session(python=py_versions[0])
def black_all(session):
    session.install("black")
    session.run("black", "--check", "--diff", "--color", "metatlas", "tests")


@nox.session(python=py_versions[0])
def black(session):
    session.install("black")
    session.run("black", "--check", "--diff", "--color", *more_checks)


@nox.session(python=py_versions[0])
def blacken(session):
    """this modifies the files to meet black's requirements"""
    session.install("black")
    session.run("black", *more_checks)


@nox.session(python=py_versions, reuse_venv=REUSE_LARGE_VENV)
def mypy(session):
    session.install("-r", "docker/requirements.txt", *mypy_deps)
    session.run("mypy", "metatlas")


@nox.session(python=py_versions, reuse_venv=REUSE_LARGE_VENV)
def pylint(session):
    session.install("-r", "docker/requirements.txt", *pylint_deps)
    session.run("pylint", *more_checks)


@nox.session(python=py_versions, reuse_venv=REUSE_LARGE_VENV)
def pylint_nb(session):
    session.install("-r", "docker/requirements.txt", *nbqa_deps, *pylint_deps)
    # dupliate code cannot be disabled on per-cell level https://github.com/PyCQA/pylint/issues/214
    # Some duplicate code is required to setup the notebook and do error handling.
    # So turn off duplicate code for whole session -- not ideal.
    session.run(
        "nbqa",
        "pylint",
        "--disable=duplicate-code",
        f"--max-line-length={NB_LINE_LEN}",
        *nbqa_flags,
        *notebooks,
    )


@nox.session(python=py_versions[0])
def flake8_nb(session):
    session.install(*nbqa_deps, *flake8_deps)
    session.run("nbqa", "flake8", *nbqa_flags, *notebooks)


@nox.session(python=py_versions[0])
def black_nb(session):
    session.install("black", *nbqa_deps)
    session.run("nbqa", "black", f"--line-length={NB_LINE_LEN}", "--check", *nbqa_flags, *notebooks)


@nox.session(python=py_versions[0])
def blacken_nb(session):
    """this modifies notebook files to meet black's requirements"""
    session.install("black", *nbqa_deps)
    session.run("nbqa", "black", f"--line-length={NB_LINE_LEN}", *nbqa_flags, *notebooks)


@nox.session(python=py_versions, reuse_venv=REUSE_LARGE_VENV)
def unit_tests(session):
    session.install("-r", "docker/requirements.txt", *pytest_deps)
    session.run(
        "pytest",
        *pytest_flags,
        *session.posargs,
        "--cov",
        "metatlas",
        "--cov-report",
        "term:skip-covered",
        "--cov-branch",
        "tests/unit/",
        env={"METATLAS_LOCAL": "TRUE"},
    )


@nox.session(python=py_versions[0], reuse_venv=REUSE_LARGE_VENV)
def cov_report(session):
    session.install("-r", "docker/requirements.txt", *pytest_deps)
    session.run(
        "pytest",
        *pytest_flags,
        *session.posargs,
        "--cov",
        "metatlas",
        "--cov-report",
        "term-missing:skip-covered",
        "--cov-branch",
        "tests/unit/",
        env={"METATLAS_LOCAL": "TRUE"},
    )


@nox.session(python=py_versions[0])
def system_tests(session):
    session.install(*pytest_deps)
    session.run("pytest", *pytest_flags, *session.posargs, "tests/system/")


@nox.session(python=py_versions[0])
def install_git_hooks(session):
    session.install("pre-commit")
    session.run("pre-commit", "install")


@nox.session(python=py_versions[0])
def update_git_hooks(session):
    session.install("pre-commit")
    session.run("pre-commit", "autoupdate")
