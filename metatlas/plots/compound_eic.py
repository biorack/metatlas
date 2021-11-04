"""Plot set of EIC plots for a compound. One sample per sub-plot"""
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


class CompoundEicPlotSet:
    """All the EIC plots for one compound"""

    # pylint: disable=too-few-public-methods,too-many-locals,too-many-arguments

    def __init__(
        self,
        data: MetatlasDataset,
        compound_idx: int,
        max_plots_per_fig: int = 30,
        sharey: bool = True,
        rt_buffer: float = 0.5,
        font_scale: 2,
    ):
        file_order_eics = [CompoundEic(sample[compound_idx]) for sample in data]
        self.eics = sorted(file_order_eics, key=lambda x: (x.short_group_name, x.short_run_id))
        num_plots = len(self.eics)
        num_pages = math.ceil(num_plots / max_plots_per_fig)
        self.plots_per_page = math.ceil(num_plots / num_pages)
        self.figures = []
        color_generator = colors()
        current_group = ""
        eic_idx = y_max = 0
        scale_factor = font_scale/num_plots**0.5
        matplotlib.rcParams.update({'font.size': 10*scale_factor})
        for _ in range(num_pages):
            plots_remaining = num_plots - eic_idx
            num_plots_this_page = min(self.plots_per_page, plots_remaining)
            cur_fig, axs = wrap_subplots(num_plots_this_page, sharey=sharey, sharex=True,
                                         constrained_layout=True)
            self.figures.append(cur_fig)
            for ax in axs:
                eic = self.eics[eic_idx]
                if eic.short_group_name != current_group:
                    current_color = next(color_generator)
                    current_group = eic.short_group_name
                eic.plot(ax, back_color=current_color, rt_buffer=rt_buffer)
                y_max = max(y_max, ax.get_ylim()[1])
                eic_idx += 1
        for fig in self.figures:
            for ax in fig.axes:
                ax.set_ylim(bottom=0, top=y_max if sharey else None)
        matplotlib.rcParams.update({'font.size': 10})


class CompoundEic:
    """EIC for one compound within a single sample"""

    def __init__(self, compound: Dict[str, Any]):
        self.compound = compound
        rt_ref: metob.RtReference = compound["identification"].rt_references[0]
        self.rt_range: Tuple[float, float] = (rt_ref.rt_min, rt_ref.rt_max)
        self.rt_peak: float = rt_ref.rt_peak

    def plot(self, ax: matplotlib.axes.Axes, back_color: Color = "white", rt_buffer=0.5) -> None:
        """
        Draw plot of EIC on ax
        back_color: background color for plot
        rt_buffer: amount of time in minutes to show to each side of rt_min/rt_max/rt_peak
        """
        ax.ticklabel_format(axis="y", scilimits=[0, 0])
        ax.set_facecolor(back_color)
        self._draw_curve(ax, rt_buffer)
        self._draw_rt_ref_lines(ax)
        self._draw_title(ax)

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
            fill_under(ax, x, y, between=self.rt_range, color="c", alpha=0.3)
        x_min = min(self.rt_range[0], self.rt_peak) - rt_buffer
        x_max = max(self.rt_range[1], self.rt_peak) + rt_buffer
        ax.set_xlim(x_min, x_max)

    def _draw_title(self, ax: matplotlib.axes.Axes) -> None:
        """Add title to plot"""
        title = f"{self.short_run_id}\n{self.short_group_name}"
        ax.set_title(title, fontdict={"fontsize": "x-small"})

    @property
    def short_run_id(self) -> str:
        """Output is RunDate_NorthenLabSampleNum_ReplicateNum_SequenceInjectionNum"""
        fields = Path(self.compound["lcmsrun"].name).stem.split("_")
        return f"{fields[0]}_{fields[11]}_{fields[13]}_{fields[15]}"

    @property
    def short_group_name(self) -> str:
        """Short name for a related group of samples"""
        return self.compound["group"].short_name


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


def save_compound_eic_pdf(
    data: MetatlasDataset,
    compound_idx: int,
    file_name: str,
    overwrite: bool = False,
    sharey: bool = True,
    max_plots_per_page: int = 30,
) -> None:
    """Create a PDF file containing the set of EIC plots for a single compound"""
    write_utils.check_existing_file(file_name, overwrite)
    compound_name = data[0][compound_idx]["identification"].compound[0].name
    plt.ioff()  # don't display the plots
    with PdfPages(file_name) as pdf:
        for fig in CompoundEicPlotSet(data, compound_idx, max_plots_per_page, sharey).figures:
            pdf.savefig(fig)
        metadata = pdf.infodict()
        metadata["Title"] = f"EICs for {compound_name}"
        metadata["Author"] = "Joint Genome Institute"
        metadata["CreationDate"] = datetime.datetime.today()
    logger.debug("Exported EIC chromatograms for %s to %s.", compound_name, file_name)
