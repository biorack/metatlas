""" For minipulating msms_refs files """

import io
import json
import logging

from typing import cast, List, Sequence, Tuple, TypedDict

import ipywidgets as widgets
import numpy as np
import pandas as pd
import traitlets

from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from traitlets import Float, HasTraits, TraitError, TraitType, validate

logger = logging.getLogger(__name__)


class Proposal(TypedDict):
    """for use with traitlets.validate"""

    owner: HasTraits
    value: object
    trait: TraitType


class Spectrum(HasTraits):
    # pylint: disable=too-few-public-methods
    """List of intensities with list of corresponding MZ values"""
    intensities: List[float] = traitlets.List(trait=Float())
    mzs: List[float] = traitlets.List(trait=Float())

    def __init__(self, mzs: Sequence[float], intensities: Sequence[float], **kwargs) -> None:
        """required fields are inputs"""
        with self.hold_trait_notifications():
            super().__init__(**kwargs)
            self.intensities = list(intensities)
            self.mzs = list(mzs)

    def __repr__(self) -> str:
        """Return representation of data"""
        nested_list_form = [[f"{m:.5f}" for m in self.mzs], [f"{x:.1f}" for x in self.intensities]]
        return str(nested_list_form).replace("'", "")

    def __str__(self) -> str:
        """Return string representation of data"""
        return self.__repr__()

    def __len__(self) -> int:
        return len(self.mzs)

    def toJSON(self) -> str:  # pylint: disable=invalid-name
        """returns JSON str representation of Spectrum"""
        return json.dumps([self.mzs, self.intensities])

    def plot(self) -> Figure:
        """returns a matplotlib Figure object showing the mz vs intensity"""
        fig, axis = plt.subplots()
        axis.vlines(self.mzs, [0] * len(self), self.intensities, colors="b", linewidth=2)
        return fig

    def widget(self) -> widgets.Image:
        """Returns a Jupyter Image widget of the mz vs intensity plot"""
        mem_fh = io.BytesIO()
        plt.ioff()
        fig = self.plot()
        fig.savefig(mem_fh, format="svg")
        plt.close(fig)
        return widgets.Image(value=mem_fh.getvalue(), format="svg+xml")

    def mz_by_intensity(self) -> List[float]:
        """m/z values ordered by intensity"""
        pair_lists = sort_mzs_intensities(self.mzs, self.intensities)
        return list(reversed(pair_lists[0]))  # pylint: disable=unsubscriptable-object

    @validate("intensities")
    def _valid_intensities(self, proposal: Proposal) -> List[float]:
        """validate positive values, not empty, and same length as mzs list"""
        value = cast(List[float], proposal["value"])
        if len(value) != len(self.mzs):
            raise TraitError("length of intensities and mzs must be equal")
        if any(x <= 0 for x in value):
            raise TraitError("intensities must be positive")
        return value

    @validate("mzs")
    def _valid_mzs(self, proposal: Proposal) -> List[float]:
        """validate positive values, not empty, and same length as intensities list"""
        value = cast(List[float], proposal["value"])
        if len(value) != len(self.intensities):
            raise TraitError("length of intensities and mzs must be equal")
        if not pd.Series(value, dtype=np.float64).is_monotonic_increasing:
            raise TraitError("mzs values must be monotonically increasing")
        if any(x <= 0 for x in value):
            raise TraitError("mzs values must be positive")
        return value


def str_to_spectrum(in_str: str) -> Spectrum:
    """Converts a spectrum string into a Spectrum class instance"""
    try:
        decoded = json.loads(in_str)
        mzs = decoded[0]
        intensities = decoded[1]
    except (TypeError, json.JSONDecodeError, IndexError):
        logger.error("Cannot convert '%s' to a Spectrum object, setting to empty spectrum", in_str)
        return Spectrum(mzs=[], intensities=[])
    if len(decoded) > 2:
        logger.error("Invalid spectrum '%s'. Truncating elements after first two lists.", in_str)
    if len(mzs) > len(intensities):
        logger.error("Invalid spectrum '%s'. Truncating mzs list as intensities list is shorter.", in_str)
        mzs = mzs[: len(intensities)]
    elif len(mzs) < len(intensities):
        logger.error("Invalid spectrum '%s'. Truncating intensities list as mzs list is shorter.", in_str)
        intensities = intensities[: len(mzs)]
    if not pd.Series(mzs, dtype=np.float64).is_monotonic_increasing:
        logger.error("Invalid spectrum '%s'. mzs values must be monotonically increasing. Sorting.", in_str)
        mzs, intensities = sort_mzs_intensities(mzs, intensities)
    try:
        return Spectrum(mzs=mzs, intensities=intensities)
    except TraitError as err:
        logger.exception("%s -- setting to empty Spectrum", err)
        return Spectrum(mzs=[], intensities=[])


def sort_mzs_intensities(mzs: List[float], intensities: List[float]) -> Tuple[List[float], List[float]]:
    """sort paired lists of mz and intensity values by intensity"""
    pairs = zip(mzs, intensities)
    sort_pairs = sorted(pairs, key=lambda x: x[0])
    return cast(Tuple[List[float], List[float]], list(zip(*sort_pairs)))
