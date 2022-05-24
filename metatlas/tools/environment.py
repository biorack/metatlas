"""
Environment setup functions

Try to keep the imports in this file to the python standard libraries.
Otherwise some of the metatlas_repo/kernel validation errors will
not correctly report problems with the notebook configuration
"""

import logging
import os
import shutil
import subprocess

from pathlib import Path

logger = logging.getLogger(__name__)


def repo_dir():
    """Returns a string with the path to the root of the Metatlas git repo"""
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def validate_data_dir(base_data_dir, experiment_id):
    """Raise FileNotFoundError if base_data_dir / experiment_id is not an existing directory"""
    experiment_dir = os.path.join(base_data_dir, experiment_id)
    try:
        if not os.path.isdir(experiment_dir):
            raise FileNotFoundError(f"Data directory does not exist at {experiment_dir}.")
    except FileNotFoundError as err:
        logger.exception(err)
        raise err


def get_repo_hash():
    """
    Returns the full hash for the current git commit or 'git not found, hash unknown'
    """
    try:
        result = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo_dir(), capture_output=True, check=True)
    except FileNotFoundError:
        return "git not found, hash unknown"
    return result.stdout.strip()


def set_git_head(source_code_version_id: str) -> None:
    """Performs a git checkout"""
    cmd = ["git", "checkout", source_code_version_id]
    subprocess.run(cmd, cwd=repo_dir(), check=True, capture_output=True)


def get_commit_date() -> str:
    """Returns a string describing when the HEAD commit was created"""
    cmd = ["git", "show", "-s", "--format=%ci -- %cr", "HEAD"]
    return subprocess.run(cmd, cwd=repo_dir(), check=True, capture_output=True, text=True).stdout.rstrip()


def install_metatlas_kernel() -> None:
    """Installs the Metatlas Targeted kernel into the user's home directory"""
    kernel_file_name = (
        Path.home() / ".local" / "share" / "jupyter" / "kernels" / "metatlas-targeted" / "kernel.json"
    )
    logger.debug("Installing jupyter kernel 'Metatlas Targeted' to %s", kernel_file_name)
    kernel_file_name.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(Path(repo_dir()) / "docker" / "shifter.kernel.json", kernel_file_name)
