"""Test of plot_set functions"""
# pylint: disable=missing-function-docstring,protected-access

import matplotlib.pyplot as plt

from metatlas.plots import plot_set


def test_autoscale():
    _, axes = plt.subplots()
    axes.plot([], [])
    plot_set._autoscale(axes)
