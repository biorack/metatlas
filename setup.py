"""Setup script for metatlas package.
"""
DISTNAME = 'metatlas'
DESCRIPTION = 'Metabolite Atlas'
LONG_DESCRIPTION = open('README.rst', 'rb').read().decode('utf-8')
MAINTAINER = 'Steven Silvester'
MAINTAINER_EMAIL = 'steven.silvester@ieee.org'
URL = 'http://github.com/biorack/metatlas'
LICENSE = 'MIT'
REQUIRES = ["numpy", "pytables", "pymzml==0.7.8", "simplejson", "rpy2", "pandas",
            "dataset", "ipython", "traitlets==4.1.0", "six", "tabulate", "dill",
            "gspread","pymysql", "qgrid", "pillow", 'oauth2client (== 1.5.2)']
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
import imp
import shutil
from setuptools import setup, find_packages
from setuptools.command.install import install


class custom_install(install):

    def run(self):
        install.run(self)
        # patch pymzml to use new obo file
        dirname = imp.find_module('pymzml')[1]
        shutil.copy('psi-ms-1.2.0.obo', '%s/obo' % dirname)
        with open('%s/obo.py' % dirname, 'r') as fid:
            lines = fid.readlines()
        with open('%s/obo.py' % dirname, 'w') as fid:
            for line in lines:
                if "version='1.1.0'" in line:
                    line = line.replace('1.1.0', '1.2.0')
                fid.write(line)


with open('metatlas/__init__.py') as fid:
    for line in fid:
        if line.startswith('__version__'):
            version = line.strip().split()[-1][1:-1]
            break


if __name__ == "__main__":

    setup(
        name=DISTNAME,
        version=version,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        url=URL,
        download_url=URL,
        license=LICENSE,
        platforms=["Any"],
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        classifiers=list(filter(None, CLASSIFIERS.split('\n'))),
        packages=find_packages(exclude=['doc']),
        include_package_data=True,
        zip_safe=False,  # the package can run out of an .egg file
        install_requires=['pymzml', 'simplejson', 'requests_toolbelt',
                          'dataset', 'ipython', 'traitlets', 'six',
                          'tabulate', 'dill', 'oauth2client==1.5.2', 'gspread',
                          'qgrid', 'pillow'],
        requires=REQUIRES,
        cmdclass={'install': custom_install},
     )
