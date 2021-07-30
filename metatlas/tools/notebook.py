"""Jupyter notebook helper functions"""

import json
import logging
import os

import pandas as pd
from IPython.core.display import display, HTML
from metatlas.tools.logging import activate_logging
from metatlas.tools.logging import activate_module_logging
from metatlas.tools.environment import validate_kernel
from metatlas.tools.environment import get_repo_hash


logger = logging.getLogger(__name__)


def configure_environment(log_level):
    """
    Sets environment variables and configures logging
    inputs:
        log_level: one of 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
    """
    activate_logging(console_level=log_level)
    logger.debug("Running import and environment setup block of notebook.")
    logger.debug("Configuring notebook environment with console log level of %s.", log_level)
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    logger.info('Running on git commit: %s', get_repo_hash())


def configure_pandas_display(max_rows=5000, max_columns=500, max_colwidth=100):
    """Set pandas display options"""
    logger.debug("Settings pandas display options")
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
    logger.debug("Activaing SQL logging with console_level=%s and file_level=%s.", console_level, file_level)
    activate_module_logging("sqlalchemy.engine", console_level, console_format, file_level, filename)


def create_notebook(input_file_name, output_file_name, parameters, injection_cell=2):
    """
    Copies from input_file_name to output_file_name and then places the parameters into a
    cell of the output notebook.
    inputs:
        input_file_name: source notebook
        output_file_name: destination notebook
        parameters: dict where keys are LHS of assignment and values are RHS of assignment
        injection_cell: zero-indexed number of cell to overwrite with the parameters
    """
    with open(input_file_name, "r") as in_fh:
        notebook = json.load(in_fh)
    notebook["cells"][injection_cell]["source"] = [assignment_string(k, v) for k, v in parameters.items()]
    with open(output_file_name, "w", encoding="utf-8") as out_fh:
        json.dump(notebook, out_fh, ensure_ascii=False, indent=4)
    logger.info("Created jupyter notebook %s", output_file_name)


def assignment_string(lhs, rhs):
    """
    inputs:
        lhs: name of variable to be assigned value
        rhs: python object that will be assigned
    returns a string
    """
    if isinstance(rhs, bool):
        rhs_str = "True" if rhs else "False"
    else:
        rhs_str = json.dumps(rhs)
    return f"{lhs} = {rhs_str}\n"
