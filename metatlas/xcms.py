
from rpy2 import robjects
from rpy2.robjects import pandas2ri


library = robjects.r['library']
library('xcms')
pandas2ri.activate()


def get_xmcs_set(files, *args, **kwargs):
    """Get an xmcs set for a group of files

    Parameters
    -----------------
    files : list of strings
        mzML files for extraction.
    *args
        Specification search terms for file selection.
    **kwargs
        Keyword arguments to xcmsSet command, such as:
        method='centWave' ppm=10, peakwidth=(5,30), snthresh=6,
        mzdiff=0.01,prefilter=(3,100)

    Returns
    -----------
    out : xcmsSet
        R xcms set.
    """
    grepl = robjects.r['grepl']
    for arg in args:
        hits = grepl(arg, files)
        files = [f for (f, h) in zip(files, hits) if h]
    paste = robjects.r['paste']
    files = paste('', files, sep='')
    xmcs_set = robjects.r['xcmsSet']
    return xmcs_set(files, **kwargs)


def group(xcms_set, **kwargs):
    """Group an xmcs set using a given set of criteria.

    Parameters
    ----------------
    xcms_set : xcmsSet
        R xcms set.
    **kwargs
        Keyword arguments to the group command, such as:
        bw=5,mzwid=0.015,minfrac=0.5,minsamp=1

    Returns
    -----------
    out : xcmsSet
        R xcms set.
    """
    grp = robjects.r['group']
    return grp(xcms_set, **kwargs)


def fill_peaks(xmcs_set):
    """Fill peaks in an xmcs_set"""
    fill = robjects.r['fillPeaks']
    return fill(xmcs_set)


def ret_cor(xmcs_set, **kwargs):
    """RETENTION CORRECTION CAN MAKE THINGS WORSE
    DONT USE IT UNLESS YOU NEED IT
    DONT USE IT UNLESS YOU KNOW WHAT YOU ARE DOING

    Example keyword parameters:
    family="symmetric", plottype="mdevden"
    """
    retcor = robjects.r['retcor']
    return retcor(xmcs_set, **kwargs)


def run_xcms(files, *args):
    """Convenience function to run xcms using default settings

    Parameters
    -----------------
    files : list of strings
        mzML files for extraction.
    *args
        Specification search terms for file selection.
    **kwargs
        Keyword arguments to xcmsSet command, such as:
        method='centWave' ppm=10, peakwidth=(5,30), snthresh=6,
        mzdiff=0.01,prefilter=(3,100)

    Returns
    -----------
    out : dataFrame
        xcms peak dataFrame.
    """
    xmcs_set = get_xmcs_set(files, *args)
    xmcs_set = group(xmcs_set)
    xmcs_set = fill_peaks(xmcs_set)
    return peak_table(xmcs_set)


def peak_table(xmcs_set, filebase='peakList'):
    """Export the global peak table

    Parameters
    ----------------
    xcms_set : xcmsSet
        R xcms set.
    filebase : str
        Type of filebase to use.

    Returns
    -----------
    out : dataFrame
        xcms peak dataFrame.
    """
    peak = robjects.r['peakTable']
    tab = peak(xmcs_set, filebase)
    df = pandas2ri.ri2py_dataframe(tab)
    df.columns = tab.colnames
    return df


if __name__ == '__main__':
    xset = get_xmcs_set(['test_basic.mzML'])
    xset = group(xset)
    df = peak_table(xset)
    print(df.head())
