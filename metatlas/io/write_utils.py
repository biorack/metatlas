""" Utility functions used in writing files"""

import logging
import os

import pandas as pd

logger = logging.getLogger(__name__)


def make_dir_for(file_path):
    """makes directories for file_path if they don't already exist"""
    directory = os.path.dirname(file_path)
    if directory != "":
        os.makedirs(directory, exist_ok=True)


def check_existing_file(file_path, overwrite=False):
    """Creates directories as needed and throws an error if file exists and overwrite is False"""
    make_dir_for(file_path)
    try:
        if not overwrite and os.path.exists(file_path):
            raise FileExistsError(f"Not overwriting {file_path}.")
    except FileExistsError as err:
        logger.exception(err)
        raise


def export_dataframe(dataframe, file_path, description, overwrite=False, **kwargs):
    """
    inputs:
        dataframe: pandas DataFrame to save
        file_path: string with path of file to create
        description: free string for logging
        overwrite: if False, raise error if file already exists
        remaining arguments are passed through to to_csv()
    """
    check_existing_file(file_path, overwrite)
    dataframe.to_csv(file_path, **kwargs)
    logger.info("Exported %s to %s.", description, file_path)


def raise_on_diff(dataframe, file_path, description, **kwargs):
    """
    inputs:
        dataframe: pandas DataFrame to save
        file_path: string with path of file to compare against
        description: free string for logging
        remaining arguments are passed through to read_csv()

    If file_path exists and does not match dataframe then raise ValueError
    """
    if not os.path.exists(file_path):
        return
    existing_df = pd.read_csv(file_path, **kwargs)
    if dataframe.equals(existing_df):
        logging.info("Data in %s is the same as %s.", description, file_path)
    else:
        try:
            raise ValueError("Data in %s is not the same as %s." % (description, file_path))
        except ValueError as err:
            logger.exception(err)
            raise


def export_dataframe_die_on_diff(dataframe, file_path, description, **kwargs):
    """
    inputs:
        dataframe: pandas DataFrame to save
        file_path: string with path of file to create
        description: free string for logging
        remaining arguments are passed through to to_csv()

    If file_path does not exist then save the dataframe there
    If file_path exists and matches data in dataframe then do nothing
    If file_path exists and does not match dataframe then raise ValueError
    """
    raise_on_diff(dataframe, file_path, description, **kwargs)
    if not os.path.exists(file_path):
        export_dataframe(dataframe, file_path, description, **kwargs)
