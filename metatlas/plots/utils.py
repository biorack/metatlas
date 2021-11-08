"""Utility functions used to create plots"""
# pylint: disable=invalid-name,too-many-arguments

import datetime
import logging
import math

from pathlib import Path
from typing import Any, Dict, Generator, List, Optional, Sequence, Tuple, Union

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.backends.backend_pdf import PdfPages

from metatlas.io import write_utils
from metatlas.datastructures import metatlas_objects as metob

logger = logging.getLogger(__name__)

BACKGROUND_COLORS = ["white", "lightyellow", "whitesmoke", "lavenderblush"]

Color = Union[str, Tuple[float, float, float, float]]

MetatlasDataset = List[List[Any]]  # avoiding a circular import


def colors() -> Generator:
    """Infinite generator that outputs a color string from BACKGROUND_COLORS"""
    num = 0
    num_colors = len(BACKGROUND_COLORS)
    while True:
        yield BACKGROUND_COLORS[num % num_colors]
        num += 1


def subplot_dimensions(num: int) -> Tuple[int, int]:
    """Returns (num_rows, num_columns)"""
    if num < 4:
        return (num, 1)
    num_cols = math.ceil(math.sqrt(num))
    num_rows = math.ceil(num / num_cols)
    return (num_rows, num_cols)


def wrap_subplots(num: int, **kwargs) -> Tuple[matplotlib.figure.Figure, List[matplotlib.axes.Axes]]:
    """Gets a figure with a grid of subplots, but internally deals with the grid dimensions"""
    nrows, ncols = subplot_dimensions(num)
    fig, axs = plt.subplots(nrows, ncols, squeeze=False, **kwargs)
    needed_axs = axs.flatten()[:num]
    blank_axs = axs.flatten()[num:]
    for ax in blank_axs:
        ax.set_axis_off()
    return (fig, needed_axs)


def is_in_range(a: List[float], start: float, stop: float) -> List[bool]:
    """Determine which members of a list fall between two values (inclusive)"""
    return [start <= x <= stop for x in a]


def fill_under(
    ax: matplotlib.axes.Axes,
    x: List[float],
    y: List[float],
    between: Optional[Tuple[float, float]] = None,
    **kwargs,
) -> None:
    """Fill under a curve with fill limited by x-range in between"""
    where = None if between is None else is_in_range(x, between[0], between[1])
    ax.fill_between(x, y, [0] * len(x), where=where, **kwargs)
