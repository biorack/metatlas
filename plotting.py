
def plot_chromatogram(x, y, title='XIC for Sample', **kwargs):
    """
    Plots a Chromatogram.

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
