"""Plot set of EIC plots for a compound. One sample per sub-plot"""
# pylint: disable=invalid-name,too-many-arguments,too-few-public-methods

import logging

from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import matplotlib
import numpy as np

from metatlas.datastructures import metatlas_objects as metob
from metatlas.plots import plot_set
from metatlas.plots import utils

logger = logging.getLogger(__name__)

MetatlasDataset = List[List[Any]]  # avoiding a circular import


class CompoundEic(plot_set.Plot):
    """EIC for one compound within a single sample"""

    def __init__(self, title: str, group_name: str, compound: Dict[str, Any], rt_buffer: float = 0.5):
        """
        compound: Compound instance
        rt_buffer: amount of time in minutes to show to each side of rt_min/rt_max/rt_peak
        """
        super().__init__(title, group_name)
        self.compound = compound
        rt_ref: metob.RtReference = compound["identification"].rt_references[0]
        self.rt_range: Tuple[float, float] = (rt_ref.rt_min, rt_ref.rt_max)
        self.rt_peak: float = rt_ref.rt_peak
        self.rt_buffer = rt_buffer

    def plot(self, ax: matplotlib.axes.Axes, back_color: utils.Color = "white") -> None:
        """Draw plot of EIC on ax"""
        super().plot(ax, back_color)
        self._draw_curve(ax, self.rt_buffer)
        self._draw_rt_ref_lines(ax)

    def _draw_rt_ref_lines(self, ax: matplotlib.axes.Axes) -> None:
        """Draw vertical lines for RT min, RT max, and RT peak"""
        for x_pos in self.rt_range:
            ax.axvline(x_pos, color="black")
        ax.axvline(self.rt_peak, color="red")

    def _draw_curve(self, ax: matplotlib.axes.Axes, rt_buffer: float) -> None:
        """Draw the EIC data and fill under the curve betweeen RT min and RT max"""
        eic = self.compound["data"]["eic"]
        if len(eic["rt"]) > 0:
            # fill_between requires a data point at each end of range, so add points via interpolation
            x, y = add_interp_at(eic["rt"], eic["intensity"], self.rt_range)
            ax.plot(x, y)
            utils.fill_under(ax, x, y, between=self.rt_range, color="c", alpha=0.3)
        x_min = min(self.rt_range[0], self.rt_peak) - rt_buffer
        x_max = max(self.rt_range[1], self.rt_peak) + rt_buffer
        ax.set_xlim(x_min, x_max)


def add_interp_at(
    x: List[float], y: List[float], x_inserts: Sequence[float]
) -> Tuple[List[float], List[float]]:
    """Assumes x and new_x are sorted"""
    x_array = np.array(x)
    x_inserts_array = np.array(x_inserts)
    full_x_array = insert_in_sorted_array(x_array, x_inserts_array)
    return full_x_array.tolist(), np.interp(full_x_array, x_array, np.array(y)).tolist()


def insert_in_sorted_array(a1: np.ndarray, a2: np.ndarray) -> np.ndarray:
    """Combine two sorted arrays into one sorrted array"""
    insert_indices = np.searchsorted(a1, a2)
    return np.insert(a1, insert_indices, a2)


def save_compound_eic_pdf(
    data: MetatlasDataset,
    compound_idx: int,
    file_name: str,
    overwrite: bool = False,
    sharey: bool = True,
    max_plots_per_page: int = 30,
) -> None:
    """Generate a PDF of EIC plots for compound_idx within data"""
    file_order_eics = []
    for sample in data:
        compound = sample[compound_idx]
        fields = Path(compound["lcmsrun"].name).stem.split("_")
        title = f"{fields[0]}_{fields[11]}_{fields[13]}_{fields[15]}"
        group_name = compound["group"].short_name
        file_order_eics.append(CompoundEic(title, group_name, compound))
    eics = sorted(file_order_eics, key=lambda x: (x.group_name, x.title))
    compound_name = data[0][compound_idx]["identification"].name
    with plot_set.PlotSet(eics, max_plots_per_page, sharey) as plots:
        plots.save_pdf(file_name, title=f"EICs for {compound_name}", overwrite=overwrite)
