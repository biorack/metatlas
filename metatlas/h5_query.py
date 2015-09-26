"""
Module used to define H5 queries.
"""
from __future__ import print_function

import tables
import numpy as np


def plot_XIC(x, y, title='XIC for Sample', **kwargs):
    """
    Plots an XIC.

    Parameters
    ----------
    x : array-like
        X values.
    y : array-like
        Y values.
    title : str, optional
        Title of plot.
    **kwargs
        Keyword arguments for ``plt.plot``.
    """
    import matplotlib.pyplot as plt
    plt.plot(x, y, **kwargs)
    plt.xlabel('Time (min)')
    plt.ylabel('Intensity')
    plt.title(title)


def plot_spectrogram(x, y, title='Spectrogram for Sample', **kwargs):
    """
    Plots a spectrogram.

    Parameters
    ----------
    x : array-like
        X values.
    y : array-like
        Y values.
    title : str, optional
        Title of plot.
    **kwargs
        Keyword arguments for ``plt.plot``.
    """
    import matplotlib.pyplot as plt
    plt.plot(x, y, **kwargs)
    plt.xlabel('Mass (m/z)')
    plt.ylabel('Intensity')
    plt.title(title)


def plot_heatmap(arr, rt_bins, mz_bins, title='HeatMap for Sample', **kwargs):
    """
    Plots the given numpy array in pyplot on a log scale (in the given format)
    and returns the image.

    Parameters
    ----------
    arr : array
        Array returned by get_HeatMapRTMZ (typically a slice).
    rt_bins : array-like
        Selected bins on the rt axis (typically a slice).
    mz_bins : array-like
        Selected bins on the mz axis (typically a slice).
    title : str, optional
        Title of plot.
    **kwargs
        Keyword arguments for ``plt.imshow``.
    """
    import matplotlib.pyplot as plt
    kwargs.setdefault('interpolation', 'nearest')
    kwargs.setdefault('aspect', 'auto')
    kwargs.setdefault('cmap', 'YlGnBu_r')
    kwargs.setdefault('origin', 'lower')

    kwargs['extent'] = [rt_bins[0], rt_bins[-1], mz_bins[0], mz_bins[-1]]

    plt.imshow(arr, **kwargs)
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)

    plt.xlabel('Time (min)')
    plt.ylabel('Mass (m/z)')
    plt.title(title)
    plt.colorbar()


def get_data(h5file, ms_level, polarity, **kwargs):
    """
    Get raw data from an h5file meeting given criteria.

    Parameters
    ----------
    h5file: string or open pytables file
        The path to the h5file or open file handle.
    ms_level : int
        MS Level.
    polarity : int
        Plus proton (1) or Minus proton (0).
    **kwargs
        Optional search modifiers.  (e.g. min_r=0, max_mz=100, precursor_MZ=1).
        Use verbose=True for displaying query messages.

    Returns
    -------
    out : dictionary
        Dictionary with arrays for 'i', 'mz', and 'rt' values meeting criteria.
    """
    if not isinstance(h5file, tables.File):
        h5file = tables.open_file(h5file)

    if ms_level == 1:
        if not polarity:
            table = h5file.root.ms1_neg
        else:
            table = h5file.root.ms1_pos
    elif not polarity:
        table = h5file.root.ms2_neg
    else:
        table = h5file.root.ms2_pos

    if not table.nrows:
        raise ValueError('No data in chosen table: %s' % table._v_name)

    query = ''

    for name in ['rt', 'mz', 'precursor_MZ', 'precursor_intensity',
                 'collision_energy']:
        if 'min_%s' % name in kwargs:
            query += ' & (%s >= %s)' % (name, kwargs['min_%s' % name])
        if 'max_%s' % name in kwargs:
            query += ' & (%s <= %s)' % (name, kwargs['max_%s' % name])
        if name in kwargs:
            query += ' & (%s == %s)' % (name, kwargs[name])

    if not query:
        data = table.read()
        return data

    # chop off the initial ' & '
    query = query[3:]
    if kwargs.get('verbose', None):
        print('Querying: %s from %s' % (query, table._v_name))

    data = table.read_where(query)

    if not data.size:
        raise ValueError('No data found matching criteria')

    if kwargs.get('verbose', None):
        print('Query complete')

    return data


def get_XIC(h5file, min_mz, max_mz, ms_level, polarity, bins=None, **kwargs):
    """
    Get Extracted-ion chromatogram (XIC) data - RT vs. cum sum of intensities

    Parameters
    ----------
    h5file: string or open pytables file
        The path to the h5file or open file handle.
    min_mz : float
        Minimum m/z value.
    max_mz : float
        Maximum m/z value.
    ms_level : int
        MS Level.
    polarity: int
        Plus proton (1) or Minus proton (0).
    bins : int or array-like, optional.
        Desired bins to use for the histogram, defaults to unique retention
        times.
    **kwargs
        Optional search modifiers.  (e.g. precursor_MZ=1,
            min_collision_energy=4)

    Returns
    -------
    out : tuple of arrays
        (rt_vals, i_vals) arrays in the desired range.
    """
    data = get_data(h5file, ms_level, polarity, min_mz=min_mz,
                    max_mz=max_mz, **kwargs)

    if not bins:
        # the last bin edget is inclusive, so we have to add another bin
        bins = np.unique(data['rt'])
        delta = bins[1] - bins[0]
        bins = np.hstack((bins, bins[-1] + delta))

    i, rt = np.histogram(data['rt'], bins=bins, weights=data['i'])

    return rt[:-1], i[:-1]


def get_heatmap(h5file, mz_bins, rt_bins, ms_level, polarity, **kwargs):
    """
    Get a HeatMap of RT vs MZ.

    Parameters
    ----------
    h5file: string or open pytables file
        The path to the h5file or open file handle.
    mz_steps : int or array-like
        Bins to use for the mz axis.
    rt_steps : int or array-like
        Bins to use for the rt axis.
    ms_level : int
        MS Level.
    polarity: int
        Plus proton (1) or Minus proton (0).
    **kwargs
        Optional search modifiers.  (e.g. precursor_MZ=1,
            min_collision_energy=4)

    Returns
    -------
    out : dict
        Dictionary containing: 'arr', 'rt_bins', 'mz_bins'.  We return a log
        scaled output array.
    """
    data = get_data(h5file, ms_level, polarity, **kwargs)

    arr, mz_bins, rt_bins = np.histogram2d(data['mz'], data['rt'],
                                           weights=data['i'],
                                           bins=(mz_bins, rt_bins))
    # center the bins
    mz_bins = (mz_bins[:-1] + mz_bins[1:]) / 2
    rt_bins = (rt_bins[:-1] + rt_bins[1:]) / 2

    mz_centroid = (np.sum(np.multiply(np.sum(arr, axis=1), mz_bins))
                   / np.sum(arr))

    return dict(arr=arr, rt_bins=rt_bins, mz_bins=mz_bins,
                mz_centroid=mz_centroid)


def get_spectrogram(h5file, min_rt, max_rt, ms_level, polarity,
                    bins=2000, **kwargs):
    """
    Get cumulative I vs MZ in RT Range (spectrogram)

    Parameters
    ----------
    h5file : table file handle
        Handle to an open tables file.
    min_rt : float
        Minimum retention time.
    max_rt : float
        Maximum retention time.
    ms_level : int
        MS Level.
    polarity: int
        Plus proton (1) or Minus proton (0).
    bins : int or array-like
        Desired bins for the histogram.
    **kwargs
        Optional search modifiers.  (e.g. precursor_MZ=1,
            min_collision_energy=4)

    Returns
    -------
    out : tuple of arrays
        (mz_vals, i_vals) arrays in the desired range.
    """
    data = get_data(h5file, ms_level, polarity, min_rt=min_rt,
                    max_rt=max_rt, **kwargs)

    i, mz = np.histogram(data['mz'], bins=bins, weights=data['i'])
    # center the bins
    mz = (mz[:-1] + mz[1:]) / 2

    return mz, i


def get_info(h5file):
    """Get info about an LCMS HDF file

    Parameters
    ----------
    h5file: string or open pytables file
        The path to the h5file or open file handle.

    Returns
    -------
    out : dict
        Number of rows for all of the tables in the file.
    """
    if not isinstance(h5file, tables.File):
        h5file = tables.open_file(h5file)

    info = dict()
    for table_name in ['ms1_neg', 'ms1_pos', 'ms2_neg', 'ms2_pos']:
        table = h5file.get_node('/%s' % table_name)
        data = dict()
        data['nrows'] = table.nrows
        if not table.nrows:
            info[table_name] = data
            continue
        data['min_mz'] = table.col('mz').min()
        data['max_mz'] = table.col('mz').max()
        data['min_rt'] = table.col('rt').min()
        data['max_rt'] = table.col('rt').max()
        info[table_name] = data

    return info


if __name__ == '__main__':  # pragma: no cover
    import argparse
    import os
    import matplotlib.pyplot as plt
    plt.iof()

    desc = "Query and plot MZML data from HDF files"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-x", "--xic", action="store_true",
                        help="Get and plot XIC")
    parser.add_argument("-s", "--spectrogram", action="store_true",
                        help="Get and plot Spectrogram")
    parser.add_argument("--heatmap", action="store_true",
                        help="Get and plot Heatmap")
    parser.add_argument('input_file', help="Input HDF file",
                        action='store', required=True)

    args = parser.parse_args()

    fname = args.input_file
    fid = tables.open_file(fname)
    basename = os.path.splitext(fname)[0]

    if args.xic:
        x, y = get_XIC(fid, 0, 100000, 1, 0)
        plot_XIC(x, y, title=basename)

    if args.spectrogram:
        x, y = get_spectrogram(fid, 1, 5, 1, 0)
        plot_spectrogram(x, y)

    if args.heatmap:
        data = get_heatmap(fid, 1000, 1000, 1, 0)
        plot_heatmap(data['arr'], data['rt_bins'], data['mz_bins'])

    plt.show()
