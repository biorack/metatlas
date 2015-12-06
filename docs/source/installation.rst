Installation
************************

Log in to https://ipython.nersc.gov/.
Run the following in a notebook cell to get set up (note this only has to be done once for a given user):

::

    %%file ~/.ipython/profile_default/startup/ipython_config.py
    import sys
    sys.path.insert(0, '/global/project/projectdirs/metatlas/anaconda/lib/python2.7/site-packages')


You can then import and use all of `metatlas`.
Workspace objects will persist between notebooks in a database file.

You can press SHIFT+TAB to get documentation while accessing the functions, including
between arguments to the function.
