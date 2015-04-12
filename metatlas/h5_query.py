"""
Module used to define H5 queries.
"""
from __future__ import print_function
import io

import sys
import tables
import numpy as np

try:
    #import matplotlib
    #matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


def plot(x, y, xlabel, ylabel, title='', image_type="png", **kwargs):
    """
    Plots the given numpy array in pyplot (in the given format) and
    returns the image

    Parameters
    ----------
    x : array-like
        X values.
    y : array-like
        Y values.
    image_type : str
        Format of image to be returned.
    **kwargs
        Keyword arguments for ``plt.plot``.

    Returns
    -------
    out : io.BytesIO
        Image in the format of image_type.
    """
    if plt is None:
        raise ValueError('Please install matplotlib')

    plt.plot(x, y, **kwargs)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    buf = io.BytesIO()
    plt.savefig(buf, format=image_type)
    buf.seek(0)
    plt.close()
    return buf


def plot_heatmap(data, vmin, vmax, hmin, hmax, image_type='png', **kwargs):
    """
    Plots the given numpy array in pyplot on a log scale (in the given format)
    and returns the image.

    Parameters
    ----------
    data : dict
        Dictionary returned by get_HeatMapRTMZ.
    vmin : int
        Minimum index along the rt axis.
    vmax : int
        Maximum index along the rt axis.
    hmin : int
        Minimum index along the mz axis.
    hmax : int
        Maximum index along the mz axis.
    image_type : str
        Format of image to be returned.
    **kwargs
        Keyword arguments for ``plt.imshow``.

    Returns
    -------
    out : io.BytesIO
        Image in the format of image_type.
    """
    if plt is None:
        raise ValueError('Please install matplotlib')

    kwargs.setdefault('interpolation', 'nearest')
    kwargs.setdefault('aspect', 'auto')
    kwargs.setdefault('cmap', 'cubehelix')

    arr = data['data'][vmin: vmax, hmin: hmax]
    hmin *= data['mz_step']
    hmax *= data['mz_step']
    vmin *= data['rt_step']
    vmax *= data['rt_step']

    kwargs['extent'] = [hmin, hmax, vmin, vmax]
    plt.imshow(np.log10(arr + 1), **kwargs)
    plt.ylabel('Time (min)')
    plt.xlabel('M/Z')
    plt.title('HeatMap for %s' % data['name'])

    buf = io.BytesIO()
    plt.savefig(buf, format=image_type)
    buf.seek(0)
    plt.close()
    return buf


def get_data(h5file, ms_level, polarity, **kwargs):
    """
    Get data from an h5file meeting given criteria.

    Parameters
    ----------
    ms_level : int
        MS Level.
    polarity : int
        Plus proton (1) or Minus proton (0).
    **kwargs
        Optional search modifiers.  (e.g. min_r=0, max_mz=100)

    Returns
    -------
    out : dictionary
        Dictionary with arrays for 'i', 'mz', and 'rt' values meeting criteria.
        Also returns the 'name' of the file.
    """
    query = '(ms_level == %s) & (polarity == %s)' % (ms_level, polarity)
    if 'min_rt' in kwargs:
        query += ' & (scan_time > %s)' % kwargs['min_rt']
    if 'max_rt' in kwargs:
        query += ' & (scan_time < %s)' % kwargs['max_rt']

    for name in ['mz', 'precursor_MZ', 'precursor_intensity',
                 'collision_energy']:
        if 'min_%s' % name in kwargs:
            query += ' & (%s > %s)' % (name, kwargs['min_%s' % name])
        if 'max_%s' % name in kwargs:
            query += ' & (%s < %s)' % (name, kwargs['max_%s' % name])

    print('Querying: %s' % query)

    table = h5file.root.spectra
    i = np.array([x['i'] for x in table.where(query)])
    rt = np.array([x['scan_time'] for x in table.where(query)])
    mz = np.array([x['mz'] for x in table.where(query)])

    return dict(i=i, rt=rt, mz=mz, name=h5file.filename.replace('.h5', ''))


def get_XICof(h5file, min_mz, max_mz, ms_level, polarity):
    """
    Get XICof data - RT vs. cumulative sum of intensities (Chromatagram)

    Parameters
    ----------
    h5file : table file handle
        Handle to an open tables file.
    min_mz : float
        Minimum m/z value.
    max_mz : float
        Maximum m/z value.
    ms_level : int
        MS Level.
    polarity: int
        Plus proton (1) or Minus proton (0).

    Returns
    -------
    out : tuple of arrays
        (rvals, ivals) arrays in the desired range.
    """
    data = get_data(h5file, ms_level, polarity, min_mz=min_mz,
                    max_mz=max_mz)

    jumps = np.nonzero(np.diff(data['rt']) > 0)[0]
    jumps = np.hstack((0, jumps, data['rt'].size - 1))

    isum = np.cumsum(data['i'])

    ivals = np.diff(np.take(isum, jumps))
    rvals = np.take(data['rt'], jumps)[1:]

    return rvals, ivals


def get_HeatMapRTMZ(h5file, mz_steps, rt_steps, ms_level, polarity):
    """
    Get a HeatMap of RT vs MZ.

    Parameters
    ----------
    h5file : table file handle
        Handle to an open tables file.
    mz_steps : int
        Sample rate in m/z.
    rt_steps : int
        Sample rate in retention time.
    ms_level : int
        MS Level.
    polarity: int
        Plus proton (1) or Minus proton (0).

    """
    data = get_data(h5file, ms_level, polarity)

    mvals = data['mz']
    ivals = data['i']
    rvals = data['rt']
    name = data['name']

    minds = np.linspace(0, mvals.size - 1, mz_steps + 1).astype(int)
    rinds = np.linspace(0, rvals.size - 1, rt_steps + 1).astype(int)

    morder = np.argsort(mvals)

    print('Building array', end='')
    sys.stdout.flush()
    arr = np.zeros((mz_steps, rt_steps))
    for ir in range(rt_steps):
        # get the intensities in this rt bin
        row = ivals[rinds[ir]: rinds[ir + 1]]
        # get the mz indices for this rt bin
        mrow = morder[rinds[ir]: rinds[ir + 1]]
        # sum the intensities within each mz bin
        for im in range(rt_steps):
            vals = row[(mrow > minds[im]) & (mrow < minds[im + 1])]
            arr[ir, im] = np.sum(vals)

        if not (ir % int(rt_steps / 10)):
            print('.', end='')
            sys.stdout.flush()

    rt_step = (rvals[-1]- rvals[0]) / rt_steps
    mz_step = (mvals.max() - mvals.min()) / mz_steps
    return dict(data=arr, rt_step=rt_step, mz_step=mz_step, name=name)


def get_IvsMZinRTRange(h5file, min_rt, max_rt, ms_level, polarity,
                       nsteps=2000):
    """
    Get cumulative I vs MZ in RT Range (Spectragram)

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
    nsteps : int
        Desired number of steps in the output array.

    """
    data = get_data(h5file, ms_level, polarity, min_rt=min_rt,
                    max_rt=max_rt)

    mvals = data['mz']
    ivals = data['i']

    order = np.argsort(mvals)
    ivals = np.take(ivals, order)
    mvals = np.take(mvals, order)

    jumps = np.linspace(0, mvals.size - 1, nsteps + 1).astype(int)
    isum = np.cumsum(ivals)

    ivals = np.diff(isum[jumps])
    mvals = mvals[jumps][:-1]

    return (mvals, ivals)


if __name__ == '__main__':
    fid = tables.open_file('140808_1_RCH2_neg.h5')
    x, y = get_XICof(fid, 1, 1000, 1, 0)
    np.save('xicof_new.npy', np.vstack((x, y)).T)
    #x, y = get_IvsMZinRTRange(fid, 1, 5, 1, 0)
    #np.save('ivsmz_new.npy', np.vstack((x, y)).T)

    #plot(x, y, 'Sum(I)', 'M/Z', 'Spectragram of %s' % fid.name)

    #data = get_HeatMapRTMZ(fid, 1000, 1000, 1, 0)
    #img = plot_heatmap(data, 0, 450, 280, 890)
    #np.save('heatmap_new.npy', (data['data'][0:450,280:890]+1)**0.1)

    #with open('test.png', 'wb') as fid:
    #    fid.write(img.read())
