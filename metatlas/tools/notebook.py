"""Jupyter notebook helper functions"""

import json
import logging
import os
import re

from os import PathLike
from typing import Any, List, Optional, Sequence

import pandas as pd
from IPython.core.display import display, HTML
from metatlas.tools.logging import activate_logging
from metatlas.tools.logging import activate_module_logging
from metatlas.tools.environment import get_commit_date
from metatlas.tools.environment import get_repo_hash
from metatlas.tools.environment import set_git_head
from metatlas.tools.environment import install_metatlas_kernel

logger = logging.getLogger(__name__)

PAPERMILL_EXEC_ENV = "PAPERMILL_EXECUTION"


def configure_environment(log_level: str) -> None:
    """
    Sets environment variables and configures logging
    inputs:
        log_level: one of 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
    """
    activate_logging(console_level=log_level)
    logger.debug("Running import and environment setup block of notebook.")
    logger.debug("Configuring notebook environment with console log level of %s.", log_level)
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    install_metatlas_kernel()


def configure_pandas_display(max_rows: int = 5000, max_columns: int = 500, max_colwidth: int = 100) -> None:
    """Set pandas display options"""
    logger.debug("Settings pandas display options")
    pd.set_option("display.max_rows", max_rows)
    pd.set_option("display.max_columns", max_columns)
    pd.set_option("display.max_colwidth", max_colwidth)


def configure_notebook_display() -> None:
    """Configure output from Jupyter"""
    # set notebook to have minimal side margins
    display(HTML("<style>.container { width:100% !important; }</style>"))


def setup(log_level: str, source_code_version_id: Optional[str] = None) -> None:
    """High level function to prepare the metatlas notebook"""
    configure_environment(log_level)
    if source_code_version_id is not None:
        set_git_head(source_code_version_id)
    logger.info("Running on git commit: %s from %s", get_repo_hash(), get_commit_date())
    configure_notebook_display()
    configure_pandas_display()


def activate_sql_logging(
    console_level: str = "INFO",
    console_format: Optional[logging.Formatter] = None,
    file_level: str = "DEBUG",
    filename: Optional[PathLike] = None,
) -> None:
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


def cells_matching_tags(data: dict, tags: Sequence[str]) -> List[int]:
    """
    For a jupyter notebook represented by data, return the list of cells with
    one or more tags that are within the tags input
    """
    return [i for i, cell in enumerate(data["cells"]) if has_intersection(get_metadata_tags(cell), tags)]


def has_intersection(one: Sequence, two: Sequence) -> bool:
    """True if the set intersection of one and two is non-empty"""
    return len(set.intersection(set(one), set(two))) > 0


def get_metadata_tags(cell: dict) -> List[str]:
    """Return a list of metadata tags for the input cell"""
    try:
        return cell["metadata"]["tags"]
    except KeyError:
        return []


def create_notebook(
    input_file_name: PathLike, output_file_name: PathLike, parameters: dict, injection_cell: int = 2
) -> None:
    """
    Copies from input_file_name to output_file_name and then places the parameters into a
    cell of the output notebook.
    inputs:
        input_file_name: source notebook
        output_file_name: destination notebook
        parameters: dict where keys are LHS of assignment and values are RHS of assignment
        injection_cell: zero-indexed number of cell to overwrite with the parameters
    """
    with open(input_file_name, "r", encoding="utf8") as in_fh:
        notebook = json.load(in_fh)
    logger.debug('Creating notebook with parameters: %s', str(parameters))
    notebook["cells"][injection_cell]["source"] = [assignment_string(k, v) for k, v in parameters.items()]
    with open(output_file_name, "w", encoding="utf-8") as out_fh:
        json.dump(notebook, out_fh, ensure_ascii=False, indent=4)
    logger.info("Created jupyter notebook %s", output_file_name)


def replace_parameters(source: Sequence[str], parameters: dict) -> List[str]:
    """Update parameter values in a list of strings and return a new list of strings"""
    eq_pat = re.compile(r"^([^#= ]+)\s*=.+$")
    out = []
    updated = []
    for line in source:
        re_match = eq_pat.match(line)
        if re_match:
            param_name = re_match.group(1)
            if param_name in parameters:
                new_value = parameters[param_name]
                out_value = f"'{new_value}'" if isinstance(new_value, str) else new_value
                out.append(f"{param_name} = {out_value}\n")
                updated.append(param_name)
                continue
        out.append(line)
    unused = set(parameters.keys()) - set(updated)
    if len(unused) > 0:
        raise ValueError(f"The following parameters could not be found in the source notebook: {unused}")
    return out


def assignment_string(lhs: str, rhs: Any) -> str:
    """
    inputs:
        lhs: name of variable to be assigned value
        rhs: python object that will be assigned
    returns a string
    """
    if rhs is None:
        rhs_str = "None"
    elif isinstance(rhs, bool):
        rhs_str = "True" if rhs else "False"
    else:
        rhs_str = json.dumps(rhs, default=vars, sort_keys=True, indent=4)
    return f"{lhs} = {rhs_str}\n"


def in_papermill():
    """
    Returns True if notebook is running in papermill.
    Requires setting of envionmental variable before running papermill!
    """
    return PAPERMILL_EXEC_ENV in os.environ
