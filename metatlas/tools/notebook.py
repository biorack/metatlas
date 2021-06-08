"""Jupyter notebook helper functions"""

import logging
import os
import shutil
import sys

from pathlib import Path
import pandas as pd
from IPython.core.display import display, HTML
from metatlas.tools.logging import activate_logging

logger = logging.getLogger(__name__)


def configure_environment(log_level):
    """
    Sets environment variables and configures logging
    inputs:
        log_level: one of 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
    """
    activate_logging(console_level=log_level)
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"


def validate_kernel():
    """
    Raise error if problem with kernel
    When on NERSC, this will install the correct kernel if needed
    """
    allowed_exe = [
        "/global/common/software/m2650/metatlas-targeted-20210521/bin/python",
    ]
    error_msg = "Invalid kernel setting in Jupyter Notebook."
    on_nersc = "METATLAS_LOCAL" not in os.environ
    if on_nersc and sys.executable not in allowed_exe:
        install_kernel()
        logger.critical('Please check that the kernel is set to "Metatlas Targeted".')
        raise ValueError(error_msg)
    try:
        # pylint: disable=import-outside-toplevel,unused-import
        import dataset  # noqa: F401
    except ModuleNotFoundError as module_error:
        logger.critical(
            'Could not find dataset module. Please check that the kernel is set to "Metatlas Targeted".'
        )
        raise ModuleNotFoundError from module_error


def install_kernel():
    """
    Copies kernel.json from repo to active location under home directory.
    Only for use on NERC!
    """
    logger.info('Installing kernel.json for "Metatlas Targeted".')
    repo_path = Path(__file__).resolve().parent.parent.parent
    source = repo_path / "notebooks" / "kernels" / "metatlas-targeted.kernel.json"
    dest_dir = Path.home() / ".local" / "share" / "jupyter" / "kernels" / "metatlas-targeted"
    os.makedirs(dest_dir, exist_ok=True)
    shutil.copyfile(source, dest_dir / "kernel.json")
    logger.info(
        (
            'Reload the page and then change kernel to "Metatlas Targeted". '
            "On the menu bar at the top of this page select 'Kernel'>'Change Kernel..' "
            "then find 'Metatlas Targeted' in the drop down list."
        )
    )


def configure_pandas_display(max_rows=5000, max_columns=500, max_colwidth=100):
    """Set pandas display options"""
    pd.set_option("display.max_rows", max_rows)
    pd.set_option("display.max_columns", max_columns)
    pd.set_option("display.max_colwidth", max_colwidth)


def configure_notebook_display():
    """Configure output from Jupyter"""
    # set notebook to have minimal side margins
    display(HTML("<style>.container { width:100% !important; }</style>"))


def setup(log_level):
    """High level function to prepare the metatlas notebook"""
    configure_environment(log_level)
    validate_kernel()
    configure_notebook_display()
    configure_pandas_display()
