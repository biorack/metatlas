"""Functions for generating mirror plots of MSMS spectra"""
import json

from typing import List, Sequence, Tuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from numpy.typing import ArrayLike

import metatlas.plots.dill2plots as dp
import metatlas.tools.spectralprocessing as sp


def mirror_plot(
    axis: matplotlib.axes.Axes,
    top: Tuple[ArrayLike, ArrayLike],
    bottom: Tuple[ArrayLike, ArrayLike],
    mz_tolerance: float = 1e-6,
    resolve_by: str = "shape",
) -> None:
    """
    Draw a mirror plot on axis

    inputs:
      axis: where to draw
      top: tuple containing an array of mz values and intensity values
      bottom: tuple containing an array of mz values and intensity values
      mz_tolerance: tolerance used for matching mz values
      resolve_by: one of 'distance', 'shape', or 'intensity' use to break ties between matching mzs
    """
    if len(top[0]) > 0 and len(bottom[0]) > 0:
        aligned_top, aligned_bottom = sp.pairwise_align_ms_vectors(top, bottom, mz_tolerance, resolve_by)
        dp.plot_msms_comparison(1, 0.0, axis, aligned_top, aligned_bottom)


def get_msms_data(
    data: pd.DataFrame, idx: int, mz_column_name: str, i_column_name: str
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extracts the msms spectra from a row in a DataFrame

    inputs:
      data: the DataFrame
      idx: row index to extract from
      mz_column_name: name of the column containing an list of mz values as a json string
      i_column_name: name of the column containing a list of intensity values as json string
    returns a tuple of numpy arrays containing mz values and intensity values
    """
    row = data.loc[idx]
    try:
        pairs = zip(json.loads(row[mz_column_name]), json.loads(row[i_column_name]))
    except TypeError:
        return np.array([]), np.array([])
    ordered = sorted(pairs, key=lambda x: x[0])
    out = list(zip(*ordered))  # inverse of zip
    return np.array(out[0]), np.array(out[1])


# pylint: disable-next=too-many-arguments
def multiple_mirror_plot(
    data: pd.DataFrame,
    top_idx: int,
    bottom_idx: int,
    mz_column_names: List[str],
    i_column_names: List[str],
    mz_tolerance: float = 1e-3,
    resolve_by: str = "shape",
    title: bool = True,
    width: float = 20,
    height: float = 12,
) -> Tuple[matplotlib.figure.Figure, Sequence[Sequence[matplotlib.axes.Axes]]]:
    """
    Generates a figure containing multiple mirror plots from two rows of a DataFrame

    inputs:
      data: the DataFrame
      top_idx: row index from which to extract values for the top half of the mirror plots
      bottom_idx: row index from which to extract values for the bottom half of the mirror plots
      mz_column_names: List of names of the columns containing an list of mz values as a json string
      i_column_name: List of names of the columns containing a list of intensity values as json string
      mz_tolerance: tolerance used for matching mz values
      resolve_by: one of 'distance', 'shape', or 'intensity' use to break ties between matching mzs
      title: display a title containing the columns names if True
      width: figure widgth
      height: figure height
    returns a tuple figure, axes where axes is always a 2D array
    """
    assert isinstance(mz_column_names, list)
    assert isinstance(i_column_names, list)
    assert len(mz_column_names) == len(i_column_names)
    fig, axes = plt.subplots(1, len(mz_column_names), squeeze=False)
    fig.set(figwidth=width, figheight=height)
    data = [
        [get_msms_data(data, idx, mz_col, i_col) for idx in [top_idx, bottom_idx]]
        for mz_col, i_col in zip(mz_column_names, i_column_names)
    ]
    for i, axis in enumerate(axes[0]):
        mirror_plot(axis, data[i][0], data[i][1], mz_tolerance, resolve_by)
        if title:
            axis.set_title(f"{mz_column_names[i]}, {i_column_names[i]}")
    return fig, axes
