""" stand alone utility functions """

from pathlib import Path
from typing import Optional, TypeVar

Generic = TypeVar("Generic")


def or_default(none_or_value: Optional[Generic], default: Generic) -> Generic:
    """
    inputs:
        none_or_value: variable to test
        default: value to return if none_or_value is None
    """
    return none_or_value if none_or_value is not None else default


def repo_path() -> Path:
    """returns Path object pointing to root directory of the git repo"""
    return Path(__file__).resolve().parent.parent.parent
