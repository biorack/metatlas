
from rpy2 import robjects
from rpy2.robjects import pandas2ri


library = robjects.r['library']
library('xcms')
pandas2ri.activate()


def get_xmcs_set(files, *args, **kwargs):
    """Get an xmcs set for a group of files"""
    grepl = robjects.r['grepl']
    for arg in args:
        hits = grepl(arg, files)
        files = [f for (f, h) in zip(files, hits) if h]
    paste = robjects.r['paste']
    files = paste('', files, sep='')
    xmcs_set = robjects.r['xcmsSet']
    return xmcs_set(files, **kwargs)


def group(xmcs_set, **kwargs):
    """Group an xmcs set using a given set of criteria"""
    grp = robjects.r['group']
    return grp(xmcs_set, **kwargs)


def fill_peaks(xmcs_set):
    """Fill peaks in an xmcs_set"""
    fill = robjects.r['fillPeaks']
    return fill(xmcs_set)


def ret_cor(xmcs_set, **kwargs):
    """RETENTION CORRECTION CAN MAKE THINGS WORSE
    DONT USE IT UNLESS YOU NEED IT
    DONT USE IT UNLESS YOU KNOW WHAT YOU ARE DOING
    """
    retcor = robjects.r['retcor']
    return retcor(xmcs_set, **kwargs)


def run_xcms(files, *args):
    """Convenience function to run xcms using default settings"""
    xmcs_set = get_xmcs_set(files, *args)
    xmcs_set = group(xmcs_set)
    xmcs_set = fill_peaks(xmcs_set)
    return peak_table(xmcs_set)


def peak_table(xmcs_set, filebase='peakList'):
    """Export the global peak table"""
    peak = robjects.r['peakTable']
    tab = peak(xmcs_set, filebase)
    df = pandas2ri.ri2py_dataframe(tab)
    df.columns = tab.colnames
    return df


if __name__ == '__main__':
    xset = get_xmcs_set(['021715_QC_6_neg.mzML'], 'neg')
    xset = group(xset)
    df = peak_table(xset)
    print(df.head())
    import ipdb; ipdb.set_trace()
