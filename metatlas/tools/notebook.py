"""Jupyter notebook helper functions"""

import logging
import os
import sys

import pandas as pd
from IPython.core.display import display, HTML
from metatlas.tools.logging import activate_logging
from metatlas.tools.logging import activate_module_logging
from metatlas.tools.environment import install_kernel


logger = logging.getLogger(__name__)


def configure_environment(log_level):
    """
    Sets environment variables and configures logging
    inputs:
        log_level: one of 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
    """
    activate_logging(console_level=log_level)
    logger.debug('Running import and environment setup block of notebook.')
    logger.debug('Configuring notebook environment with console log level of %s.', log_level)
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
    logger.debug('Kernel validation passed. Using python from %s.', sys.executable)


def configure_pandas_display(max_rows=5000, max_columns=500, max_colwidth=100):
    """Set pandas display options"""
    logger.debug('Settings pandas display options')
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


def activate_sql_logging(console_level="INFO", console_format=None, file_level="DEBUG", filename=None):
    """
    Turns on logging from sqlalchemy.
    Level 'INFO' gets SQL statements and 'DEBUG' gets SQL statements and results.
    inputs:
        console_level: one of 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
        console_format: input to logging.setFormatter
        file_level: one of 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
        filename: logging destination

    """
    logger.debug('Activaing SQL logging with console_level=%s and file_level=%s.', console_level, file_level)
    activate_module_logging("sqlalchemy.engine", console_level, console_format, file_level, filename)
