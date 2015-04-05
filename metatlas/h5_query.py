"""
Module used to define H5 queries.
"""
import io

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


def plot(self, np_arr, image_type="png", **kwargs):
    """
    Plots the given numpy array in pyplot (in the given format) and
    returns the image

    Parameters
    ----------
    np_arr : array-like
        Numpy array to be graphed in the form of [(x,y),...].
    image_type : str
        Format of image to be returned.

    Returns
    -------
    out : io.BytesIO
        Image in the format of image_type.
    """
    plt.plot(np_arr[:, 0], np_arr[:, 1])
    buf = io.BytesIO()
    plt.savefig(buf, format=image_type)
    buf.seek(0)
    plt.close()
    return buf


def plot_log(self, np_arr, image_type, **kwargs):
    plt.imshow(np.log10(np_arr+1))
    buf = io.BytesIO()
    plt.savefig(buf, format=image_type)
    buf.seek(0)
    plt.close()
    return buf


def get_XICof(h5file, min_mz, max_mz, level, polarity):
    """
    Get XICof data

    Parameters
    ----------
    h5file : table file handle
        Handle to an open tables file.
    min_mz : float
        Minimum m/z value.
    max_mz : float
        Maximum m/z value.
    level : int
        Level.
    polarity: int
        Plus proton (1) or Minus proton (0).

    """
    # Ensures all parameters are in the expected format
    L = int(level)
    P = int(polarity)
    max_mz = float(max_mz)
    min_mz = float(min_mz)

    mindlogmz = int(np.log10(min_mz)*1000000)
    maxdlogmz = int(np.log10(max_mz)*1000000)

    # Takes a chunk of data between MIN_MZ and MAX_MZ
    btw = "between(%(ARR)s, %(L)d, %(P)d, %(FILE)s, null, null, %(MIN_DLOGMZ)d,   %(L)d, %(P)d, %(FILE)s, null, null, %(MAX_DLOGMZ)d)" % {
                'ARR': arrayname,
                'FILE': fileid,
                'MIN_DLOGMZ': mindlogmz,
                'MAX_DLOGMZ': maxdlogmz,
                'L': L,
                'P': P}

    # Calculates nsteps
    aggr1 = "aggregate(%s, sum(I) as sumI, D_RT)" % btw
    count = "count(%s)" % aggr1
    res_count = sdb._execute_query(count, response=True, n=0)
    nsteps = (eval(res_count)[0])

    # Get values as [RT, I]
    redim_arr = "<D_RT: int64, sumI: double NULL DEFAULT null> [x=0:*,2000000,0]"
    redim = "redimension(%s, %s)" % (aggr1, redim_arr)
    app = "apply(%s, RT, double(D_RT/1000000.0))" % redim
    proj = "project(%s, RT, sumI)" % app
    querystr = proj

    res = sdb._execute_query(querystr, response=True, n=0, fmt="csv")
    c = StringIO(res)
    nparr = np.loadtxt(c, delimiter=",", skiprows=1)
    return nparr


def get_XICof_mf(h5file, min_mz, max_mz, min_rt, max_rt, level, polarity):
    """
    Get XICof_mf data

    Parameters
    ----------
    h5file : table file handle
        Handle to an open tables file.
    min_mz : float
        Minimum m/z value.
    max_mz : float
        Maximum m/z value.
    min_rt : float
        Minimum RT value minutes
    max_rt : float
        Maximum RT value minutes
    level : int
        Level.
    polarity: int
        Plus proton (1) or Minus proton (0).

    """
    # Ensures all parameters are in the expected format
    P = int(polarity)
    L = int(level)
    max_mz = float(max_mz)
    min_mz = float(min_mz)
    max_rt = float(max_rt)
    min_rt = float(min_rt)

    fileidlist=str(fileidlist)

    arrayname = unicode(arrayname)

    mindlogmz = int(np.log10(min_mz)*1000000)
    maxdlogmz = int(np.log10(max_mz)*1000000)

    minrt=int(min_rt*1000000)
    maxrt=int(max_rt*1000000)

    sdb = self.connect()

    #STARTSTART
    #now create a small array with fileIDs to cross_join with:
    filesdb=sdb.from_array(np.loadtxt(StringIO(fileidlist),delimiter=',',dtype=np.int64))

    #Hard coding dimension info for now, can query the array to get dynamically

    #using a string dictionary to build the final query
    sdict={'FL_ARR':filesdb.name,
           'ARR':arrayname,
           'MIN_DLOGMZ': mindlogmz,
           'MAX_DLOGMZ': maxdlogmz,
           'MIN_RT': minrt,
           'MAX_RT': maxrt,
           'L': L,
           'P': P}
    sdict['FL_APPLY']="apply(%(FL_ARR)s,D_file,f0)"%sdict
    sdict['FL_REDIM']="redimension(%(FL_APPLY)s,<f0:int64>[D_file=0:100000,10,0])"%sdict
    sdict['CROSSJOIN']="cross_join(%(ARR)s as l, %(FL_REDIM)s as r,l.D_file, r.D_file)"%sdict

    # Takes a chunk of data between MIN_MZ and MAX_MZ

    sdict['btw'] = "between(%(CROSSJOIN)s, %(L)d, %(P)d, null, %(MIN_RT)d, null, %(MIN_DLOGMZ)d, %(L)d, %(P)d, null, %(MAX_RT)d, null, %(MAX_DLOGMZ)d)" % sdict

    # Calculates nsteps
    sdict['aggr1'] = "aggregate(%(btw)s, sum(I) as sumI, D_RT,D_file)" % sdict
    sdict['count'] = "count(%(aggr1)s)" % sdict

#res_count = sdb._execute_query(sdict['count'], response=True, n=0)
#nsteps = (eval(res_count)[0])

# Get values as [RT, I]
    sdict['redim_arr'] = "<D_RT: int64, D_file: int64, sumI: double NULL DEFAULT null> [x=0:*,2000000,0]"
    sdict['redim'] = "redimension(%(aggr1)s, %(redim_arr)s)" % sdict
    sdict['app'] = "apply(%(redim)s, RT, double(D_RT/1000000.0))" % sdict
    sdict['proj'] = "project(%(app)s, RT, sumI,D_file)" % sdict
    querystr = "count(%s)"%sdict['proj']
    querystr = sdict['proj']

    res = sdb._execute_query(querystr, response=True, n=0, fmt="csv")
    c = StringIO(res)
    nparr = np.loadtxt(c, delimiter=",", skiprows=1)

    #ENDEND

    return nparr


def get_HeatMapRTMZ(h5file, mz_steps, rt_steps, level, polarity):
    """
    Get HeatMapRTMZ data

    Parameters
    ----------
    h5file : table file handle
        Handle to an open tables file.
    mz_steps : int
        Sample rate in m/z
    rt_steps : int
        Sample rate in retention time
    level : int
        Level.
    polarity: int
        Plus proton (1) or Minus proton (0).

    """
    # Make sure all parameters in the expected format
    P = int(polarity)
    L = int(level)
    mz_steps = int(mz_steps)
    rt_steps = int(rt_steps)
    fileid = int(fileid)
    arrayname = unicode(arrayname)

    sdb = self.connect()

    MZbins = (3500000-1000000)/mz_steps
    RTbins = (80000000)/rt_steps

    btw = "between(%(ARR)s, %(L)d, %(P)d, %(FILE)d, null, null, null, %(L)d, %(P)d, %(FILE)d, null, null, null)" % {
        'ARR':arrayname,
        'FILE':fileid,
        'L':L,
        'P':P}

    aggr1 = "aggregate(%s, sum(I) as sumI, D_RT, D_LOGMZ)" % btw
    regrid = "regrid(%s, %d, %d, sum(sumI) as sumI)" % (aggr1, RTbins, MZbins)

    zeros = "build(<sumI:double NULL DEFAULT null> [D_RT=0:%d,2000000,0, D_LOGMZ=1000000:%d,200000,0],0)" % (rt_steps, mz_steps+1000000)

    merge = "merge(%s,%s)" % (regrid,zeros)


    app = "apply(%s, MZ, double(pow(10,1.6+(D_LOGMZ-1000000.0)/%d*(3.3-1.6))), RT, double(D_RT*%d/10000.0))" % (merge, mz_steps, RTbins)

    proj2 = "project(%s, sumI)" % app

    # Port array to numpy
    schema = sdb._execute_query("show('%s', 'afl')" % proj2, n=0, response=True)
    datashape = SciDBDataShape.from_schema(eval(schema)[0])
    bytes_rep = sdb._execute_query(proj2, response=True, n=0, fmt=datashape.sdbtype.bytes_fmt)
    bytes_arr = np.atleast_1d(np.fromstring(bytes_rep, dtype=datashape.dtype))
    arr = bytes_arr.reshape(datashape.shape)

    return arr


def get_IvsMZinRTRange(h5file, min_rt, max_rt, level, polarity):
    """
    Get IvsMZinRTRange data

    Parameters
    ----------
    h5file : table file handle
        Handle to an open tables file.
    min_rt : float
        Minimum retention time.
    max_rt : float
        Maximum retention time.
    level : int
        Level.
    polarity: int
        Plus proton (1) or Minus proton (0).

    """
    sdb = self.connect()

    arrayname = unicode(arrayname)
    fileid = int(fileid)
    min_rt = float(min_rt)
    max_rt = float(max_rt)
    L = int(level)
    P = int(polarity)

    minDRT=int(min_rt*1000000)
    maxDRT=int(max_rt*1000000)

    # Takes a chunk of data from MIN_RT to MAX_RT
    btw="between(%(ARR)s, %(L)d, %(P)d, %(FILE)d, %(MIN_DRT)d, null, null, %(L)d, %(P)d, %(FILE)d, %(MAX_DRT)d, null, null)"%{
        'ARR':arrayname,
        'FILE':fileid,
        'MIN_DRT':minDRT,
        'MAX_DRT':maxDRT,
        'L':L,
        'P':P}

    # Counts available RTs inside region
    # nrt="apply(aggregate(%s, count(*) as ct, D_RT),RT, int64(D_RT/10000.))"%btw
    # res=sdb._execute_query(nrt,response=True, n=0, fmt="csv")

    # Sums the intensities across M/Z
    aggr1="aggregate(%s, sum(I) as sumI, D_LOGMZ)"%btw
    # Scales M/Z
    apply="apply(%s, MZ, pow(10,D_LOGMZ/1000000.))"%aggr1

    querystr = "project(%s, MZ, sumI)" % apply

    # Execute the query
    res=sdb._execute_query(querystr,response=True, n=0, fmt="csv")
    c = StringIO(res)
    nparr = np.loadtxt(c, delimiter=",", skiprows=1)
    return nparr
