"""Set of plots"""
# pylint: disable=invalid-name,too-many-arguments,too-few-public-methods

import datetime
import logging
import math

from abc import ABC
from typing import Optional, Sequence, Tuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

from matplotlib.backends.backend_pdf import PdfPages

from metatlas.io import write_utils
from metatlas.plots import utils

logger = logging.getLogger(__name__)


def restore_matplotlib_settings(func):
    """Stop function from permanently modifiying rcParams"""

    def wrapper(*args, **kwargs):
        original = matplotlib.rcParams
        out = func(*args, **kwargs)
        matplotlib.rcParams = original
        return out

    return wrapper


class Plot(ABC):
    """A single plot"""

    def __init__(self, title: str, group_name: str):
        self.title = title
        self.group_name = group_name

    def plot(self, ax: matplotlib.axes.Axes, back_color: utils.Color = "white") -> None:
        """
        Draw plot on ax
        back_color: background color for plot
        """
        ax.ticklabel_format(axis="y", scilimits=[0, 0])
        ax.set_facecolor(back_color)
        ax.set_title(f"{self.title}\n{self.group_name}")


class PlotSet(ABC):
    """A Set of related plots"""

    # pylint: disable=too-few-public-methods,too-many-locals,too-many-arguments

    @restore_matplotlib_settings
    def __init__(
        self,
        plots: Sequence[Plot],
        max_plots_per_fig: int = 30,
        x_min: Optional[float] = None,
        x_max: Optional[float] = None,
        y_min: Optional[float] = None,
        y_max: Optional[float] = None,
        sharey: bool = True,
        font_scale: float = 2,
    ):
        num_plots = len(plots)
        num_pages = math.ceil(num_plots / max_plots_per_fig)
        plots_per_page = math.ceil(num_plots / num_pages)
        self.figures = []
        color_generator = utils.colors()
        current_group = None
        plot_idx = 0
        for _ in range(num_pages):
            plots_remaining = num_plots - plot_idx
            num_plots_this_page = min(plots_per_page, plots_remaining)
            _set_font_size(num_plots_this_page, font_scale)
            cur_fig, axs = utils.wrap_subplots(
                num_plots_this_page, sharey=sharey, sharex=True, constrained_layout=True
            )
            self.figures.append(cur_fig)
            for ax in axs:
                if current_group is None or plots[plot_idx].group_name != current_group:
                    current_color = next(color_generator)
                    current_group = plots[plot_idx].group_name
                plots[plot_idx].plot(ax, back_color=current_color)
                plot_idx += 1
        if sharey:
            self.sharey_between_pages()
        self.limit_axes(x_min, x_max, y_min, y_max)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        self.close()

    def sharey_between_pages(self):
        min_y = max_y = None
        for fig in self.figures:
            for ax in fig.axes:
                ax_min, ax_max = ax.get_ylim()
                if min_y is None or ax_min < min_y:
                    min_y = ax_min
                if max_y is None or ax_max > max_y:
                    max_y = ax_max
        if min_y is not None and max_y is not None:
            for fig in self.figures:
                for ax in fig.axes:
                    ax.set_ylim(bottom=min_y, top=max_y)

    def limit_axes(
        self,
        x_min: Optional[float],
        x_max: Optional[float],
        y_min: Optional[float],
        y_max: Optional[float],
    ) -> None:
        """Update all plots to have the desired axes limits"""
        if x_min is None and x_max is None and y_min is None and y_max is None:
            return
        sides = None
        if x_min is not None or x_max is not None:
            if y_min is None and y_max is None:
                sides = "both"
            elif y_min is None and y_max is not None:
                sides = "min"
            elif y_min is not None and y_max is None:
                sides = "max"
        for fig in self.figures:
            for ax in fig.axes:
                ax.set_xlim(left=x_min, right=x_max)
                ax.set_ylim(bottom=y_min, top=y_max)
                if sides is not None:
                    _autoscale(ax, axis="y", sides=sides)

    def save_pdf(self, file_name: str, title: str, overwrite: bool = False) -> None:
        """Create a PDF file containing one figure per page"""
        write_utils.check_existing_file(file_name, overwrite)
        plt.ioff()  # don't display the plots
        matplotlib.use("agg")  # mitigates a memory leak to not use backend_nbagg
        with PdfPages(file_name) as pdf:
            for fig in self.figures:
                pdf.savefig(fig)
            metadata = pdf.infodict()
            metadata["Title"] = title
            metadata["Author"] = "Joint Genome Institute"
            metadata["CreationDate"] = datetime.datetime.today()
        logger.debug("Exported PDF containing %s to %s.", title, file_name)

    def close(self):
        """Close all plots and free their memory"""
        for fig in self.figures:
            plt.close(fig)


# adapted from https://stackoverflow.com/questions/51323505/how-to-make-relim-and-autoscale-in-a-scatter-plot
def _get_xy(artist: matplotlib.artist.Artist) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Gets the xy coordinates of a given artist"""
    if "Collection" in str(artist):
        return artist.get_offsets().T
    if "Line" in str(artist):
        return artist.get_xdata(), artist.get_ydata()
    raise ValueError("This type of object isn't implemented yet")


# adapted from https://stackoverflow.com/questions/51323505/how-to-make-relim-and-autoscale-in-a-scatter-plot
def _autoscale(ax: matplotlib.axes.Axes, axis: str = "y", sides: str = "both", margin: float = 0.1) -> None:
    """Autoscales the x or y axis of a given matplotlib ax object
    to fit the margins set by manually limits of the other axis,
    with margins in fraction of the width of the plot
    if sides is 'max' or 'min' then only adjust the limit on that side of axis"""
    assert axis in ["x", "y"]
    assert sides in ["both", "min", "max"]
    low, high = np.inf, -np.inf
    for artist in ax.collections + ax.lines:
        if axis == "y":
            set_lim = ax.set_ylim
            get_lim = ax.get_ylim
            cur_fixed_limit = ax.get_xlim()
            fixed, dependent = _get_xy(artist)
        else:
            set_lim = ax.set_xlim
            get_lim = ax.get_xlim
            cur_fixed_limit = ax.get_ylim()
            dependent, fixed = _get_xy(artist)
        low, high = _update_limts(low, high, fixed, dependent, cur_fixed_limit)
    margin = margin * (high - low)
    if low == np.inf and high == -np.inf:
        return
    assert low != np.inf and high != -np.inf
    new_min = (low - margin) if sides in ["both", "min"] else get_lim()[0]
    new_max = (high + margin) if sides in ["both", "max"] else get_lim()[1]
    set_lim(new_min, new_max)


def _update_limts(
    low: float,
    high: float,
    fixed: npt.NDArray[np.float64],
    dependent: npt.NDArray[np.float64],
    fixed_limits: Tuple[float, float],
) -> Tuple[float, float]:
    """Current limits low and high are updated to include data with ranges
    in the lists fixed and dependent subject to fixed_limits
    """
    local_low, local_high = _calculate_new_limit(fixed, dependent, fixed_limits)
    return min(local_low, low), max(local_high, high)


# adapted from https://stackoverflow.com/questions/51323505/how-to-make-relim-and-autoscale-in-a-scatter-plot
def _calculate_new_limit(
    fixed: npt.NDArray[np.float64], dependent: npt.NDArray[np.float64], fixed_limits: Tuple[float, float]
) -> Tuple[float, float]:
    """Calculates the min/max of the dependent axis given a fixed axis with limits"""
    if len(fixed) > 2:
        mask = (fixed > fixed_limits[0]) & (fixed < fixed_limits[1])
        window = dependent[mask]
        if len(window) == 0:
            return np.inf, -np.inf
        return window.min(), window.max()
    low = dependent[0]
    high = dependent[-1]
    if low == 0.0 and high == 1.0:
        # This is a axhline in the autoscale direction
        return np.inf, -np.inf
    return low, high


def _set_font_size(num_plots: int, font_scale: float) -> None:
    """Scales the font size down based on the number of plots"""
    nrows, ncols = utils.subplot_dimensions(num_plots)
    scale_factor = font_scale / (nrows * ncols)**0.5
    matplotlib.rcParams.update({"font.size": 10 * scale_factor})
