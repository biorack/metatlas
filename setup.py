"""Setup script for metatlas package.
"""
DISTNAME = 'metatlas'
DESCRIPTION = 'Metabolite Atlas'
LONG_DESCRIPTION = open('README.rst', 'rb').read().decode('utf-8')
MAINTAINER = 'Steven Silvester'
MAINTAINER_EMAIL = 'steven.silvester@ieee.org'
URL = 'http://github.com/metabolite-atlas/metatlas'
LICENSE = 'MIT'
REQUIRES = ["numpy", "pytables", "pymzml"]
PACKAGE_DATA = {DISTNAME: ['install/*.*']}
CLASSIFIERS = """\
Development Status :: 2 - Pre-Alpha
Intended Audience :: Developers
Intended Audience :: Science/Research
License :: OSI Approved :: BSD License
Operating System :: OS Independent
Programming Language :: Python
Programming Language :: Python :: 2.7
Programming Language :: Python :: 3.4
Topic :: Scientific/Engineering
Topic :: Software Development
"""
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('metatlas/__init__.py', 'rb') as fid:
    for line in fid:
        line = line.decode('utf-8')
        if line.startswith('__version__'):
            version = line.strip().split()[-1][1:-1]
            break


if __name__ == "__main__":

    setup(
        name=DISTNAME,
        version=version,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        package_data=PACKAGE_DATA,
        url=URL,
        download_url=URL,
        license=LICENSE,
        platforms=["Any"],
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        classifiers=list(filter(None, CLASSIFIERS.split('\n'))),
        install_requires=['pymzml'],
        requires=REQUIRES
     )

    import os
    os.system('cd install; bash install.sh')
