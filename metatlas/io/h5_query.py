"""
Module used to define H5 queries.
"""
from __future__ import print_function

import tables
import numpy as np


def get_data(h5file, **kwargs):
    """
    Get raw data from an h5file meeting given criteria.

    Parameters
    ----------
    h5file: string or open pytables file
        The path to the h5file or open file handle.
    **kwargs
        Optional search modifiers.  (e.g. min_r=0, max_mz=100, precursor_MZ=1).
        Use verbose=True for displaying query messages.

    Returns
    -------
    out : dictionary
        Dictionary with arrays for 'i', 'mz', and 'rt' values meeting criteria.
    """
    is_h5_file = False
    if not isinstance(h5file, tables.File):
        is_h5_file = True
        h5file = tables.open_file(h5file)


    # Select the ms_level
    ms_level = kwargs.get('ms_level', None)
    if ms_level is None:
        ms1 = h5file.root.ms1_neg.nrows + h5file.root.ms1_pos.nrows
        ms2 = h5file.root.ms2_neg.nrows + h5file.root.ms2_pos.nrows
        ms_level = 1 if ms1 > ms2 else 2

    # Select the polarity
    polarity = kwargs.get('polarity', None)
    if polarity is None:
        if ms_level == 1:
            polarity = h5file.root.ms1_pos.nrows > h5file.root.ms1_neg.nrows
        else:
            polarity = h5file.root.ms2_pos.nrows > h5file.root.ms2_neg.nrows

    # Select the table
    if ms_level == 1:
        name = 'ms1_pos' if polarity else 'ms1_neg'
    else:
        name = 'ms2_pos' if polarity else 'ms2_neg'
    if 'rt' in kwargs or 'min_rt' in kwargs or 'max_rt' in kwargs:
        data_table = h5file.get_node('/' + name)
    else:
        try:
            data_table = h5file.get_node('/' + name + '_mz')
        except Exception:
            # Fall back on original node.
            data_table = h5file.get_node('/' + name)

    # Get the selected entries
    queries = []
    for name in ['rt', 'mz', 'precursor_MZ', 'precursor_intensity',
                 'collision_energy']:
        if 'min_%s' % name in kwargs:
            queries.append('(%s >= %s)' % (name, kwargs['min_%s' % name]))
        if 'max_%s' % name in kwargs:
            queries.append('(%s <= %s)' % (name, kwargs['max_%s' % name]))
        if name in kwargs:
            queries.append('(%s == %s)' % (name, kwargs[name]))
    query = ' & '.join(queries)

    if kwargs.get('verbose', None):
        print('Querying: %s from %s' % (query, data_table._v_name))

    if not query:
        data = data_table.read()
    else:
        data = data_table.read_where(query)

    if not data.size:
        raise ValueError('No data found matching criteria')

    if kwargs.get('verbose', None):
        print('Query complete')

    # close the file if opened
    if is_h5_file:
        h5file.close()

    return data


def get_chromatogram(h5file, min_mz, max_mz, aggregator=np.sum, **kwargs):
    """
    Get Chromatogram data - RT vs. intensity aggregation

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
    aggregator: function
        Function to aggregate the intensity data.  Defaults to np.sum,
        producing an XIC. For Base Peak Chromatogram, use np.max.
    **kwargs
        Optional search modifiers.  (e.g. precursor_MZ=1,
            min_collision_energy=4)

    Returns
    -------
    out : tuple of arrays
        (rt_vals, i_vals) arrays in the desired range.
    """
    data = get_data(h5file, min_mz=min_mz, max_mz=max_mz, **kwargs)
    if data is None:
        return [], []

    rt = np.unique(data['rt'])
    if len(rt)==1:
        return [rt], [np.sum(data['i'])]
    if aggregator == np.sum:
        d = np.diff(rt) / 2
        edges = np.hstack([rt[0] - d[0], rt[0:-1] + d, rt[-1] + d[-1]])
        i, _ = np.histogram(data['rt'], bins=edges, weights=data['i'])
    else:
        i = []
        for val in rt:
            indices = np.argwhere(data['rt'] == val)
            i.append(aggregator(np.take(data['i'], indices)))
        i = np.array(i)
    return rt, i


def get_heatmap(h5file, mz_bins, **kwargs):
    """
    Get a HeatMap of RT vs MZ.

    Parameters
    ----------
    h5file: string or open pytables file
        The path to the h5file or open file handle.
    mz_steps : int or array-like
        Bins to use for the mz axis.
    **kwargs
        Optional search modifiers.  (e.g. precursor_MZ=1,
            min_collision_energy=4)

    Returns
    -------
    out : dict
        Dictionary containing: 'arr', 'rt_bins', 'mz_bins'.
    """
    data = get_data(h5file, **kwargs)
    if data is None:
        return None

    rt_values = np.unique(data['rt'])
    rt_bins = np.hstack((rt_values, rt_values[-1] + 1))
    arr, mz_bins, _ = np.histogram2d(data['mz'], data['rt'],
                                     weights=data['i'],
                                     bins=(mz_bins, rt_bins))

    mz_centroid = (np.sum(np.multiply(np.sum(arr, axis=1), mz_bins[:-1]))
                   / np.sum(arr))

    return dict(arr=arr, rt_bins=rt_values, mz_bins=mz_bins,
                mz_centroid=mz_centroid)


def get_spectrogram(h5file, min_rt, max_rt, bins=2000, **kwargs):
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
    data = get_data(h5file, min_rt=min_rt, max_rt=max_rt, **kwargs)
    if data is None:
        return [], []

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
    is_h5_file = False
    if not isinstance(h5file, tables.File):
        is_h5_file = True
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

    # close the file if opened
    if is_h5_file:
        h5file.close()


    return info


if __name__ == '__main__':  # pragma: no cover
    from os.path import expandvars
    import sys

    # sys.path.insert(0,'/global/homes/b/bpb/metatlas/metatlas/helpers/' )
    sys.path.insert(0, expandvars('$HOME/dev/metatlas'))

    import argparse
    import os
    import matplotlib.pyplot as plt
    from metatlas import plot_chromatogram, plot_spectrogram, plot_heatmap

    desc = "Query and plot MZML data from HDF files"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-x", "--xic", action="store_true",
                        help="Get and plot XIC")
    parser.add_argument("-s", "--spectrogram", action="store_true",
                        help="Get and plot Spectrogram")
    parser.add_argument("--heatmap", action="store_true",
                        help="Get and plot Heatmap")
    parser.add_argument('input_file', help="Input HDF file",
                        action='store')

    args = parser.parse_args()

    fname = args.input_file
    fid = tables.open_file(fname)
    basename = os.path.splitext(fname)[0]

    if args.xic:
        x, y = get_chromatogram(fid, 0, 1000)
        plot_chromatogram(x, y, title=basename)

    if args.spectrogram:
        x, y = get_spectrogram(fid, 1, 5)
        plot_spectrogram(x, y)

    if args.heatmap:
        data = get_heatmap(fid, 1000)
        plot_heatmap(data['arr'], data['rt_bins'], data['mz_bins'])

    # close the file
    fid.close()

    plt.show()
