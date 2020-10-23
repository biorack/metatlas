"""Setup script for metatlas package.
"""
DISTNAME = 'metatlas'
DESCRIPTION = 'Metabolite Atlas'
LONG_DESCRIPTION = open('README.rst', 'rb').read().decode('utf-8')
MAINTAINER = 'Ben Bowen'
MAINTAINER_EMAIL = 'ben.bowen@gmail.com'
URL = 'http://github.com/biorack/metatlas'
LICENSE = 'MIT'
# REQUIRES = ["pyyaml", "pytables", "pymzml", "simplejson", "rpy2",
#             "dataset", "traitlets", "six", "dill",
#             "gspread","pymysql", "qgrid", "pillow", 'oauth2client (== 1.5.2)']
CLASSIFIERS = """\
Development Status :: 2 - Pre-Alpha
Intended Audience :: Developers
Intended Audience :: Science/Research
License :: OSI Approved :: BSD License
Operating System :: OS Independent
Programming Language :: Python
Programming Language :: Python :: 3.8
Topic :: Scientific/Engineering
Topic :: Software Development
"""
import shutil
from importlib import util

from setuptools import setup
# , find_packages
from setuptools.command.install import install

#Make sure /global/common/software/m2650/nersc_python_install/metatlas/metatlas/helpers/isotope_dict.pkl /global/common/software/m2650/python3-cori/lib/python3.7/site-packages/metatlas/helpers/

# also copy /global/homes/b/bpb/repos/metatlas/metatlas/tools/isotope_dict.json

class custom_install(install):

    def run(self):
        install.run(self)
        # patch pymzml to use new obo file
#         dirname = imp.find_module('pymzml')[1]
        dirname = util.find_spec("pymzml").submodule_search_locations[0]
        shutil.copy('data/other/psi-ms-1.2.0.obo', '%s/obo' % dirname)
        with open('%s/obo.py' % dirname, 'r') as fid:
            lines = fid.readlines()
        with open('%s/obo.py' % dirname, 'w') as fid:
            for line in lines:
                if "version='1.1.0'" in line:
                    line = line.replace('1.1.0', '1.2.0')
                fid.write(line)




if __name__ == "__main__":
    with open('metatlas/__init__.py') as fid:
        for line in fid:
            if line.startswith('__version__'):
                version = line.strip().split()[-1][1:-1]
                break

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
#         packages=find_packages(exclude=['doc']),
#         include_package_data=True,
        zip_safe=False,  # the package can run out of an .egg file
#         install_requires=REQUIRES,
#         install_requires=['pymzml', 'pyyaml','simplejson', 'requests_toolbelt',
#                           'dataset', 'ipython', 'traitlets', 'six',
#                           'tabulate', 'dill', 'oauth2client==1.5.2', 'gspread',
#                           'qgrid', 'pillow'],
#         requires=REQUIRES,
        cmdclass={'install': custom_install},
     )
