# pylint: disable=missing-function-docstring

import nox

py_versions = ['3.8', '3.9']

nox.options.sessions = ["flake8_diff", "unit_tests-3.8", "update_git_hooks"]

single_py_args = {'venv_backend': 'conda', 'reuse_venv': False}
multi_py_args = {**single_py_args, **{'python': py_versions}}

nox.options.error_on_external_run = True

run_deps = [
    'alembic=1.5.8',
    'banal=1.0.6',
    'dill=0.3.3',
    'gspread=3.7.0',
    'hdf5=1.10.6',
    'ipywidgets=7.6.3',
    'jupyterlab=3.0.14',
    'matplotlib=3.4.1',
    'oauth2client=4.1.3',
    'pandas=1.2.4',
    'pip=21.1',
    'pymysql=1.0.2',
    'pymzml=2.4.7',
    'pytables=3.6.1',
    'pyyaml=5.4.1',
    'rdkit=2021.03.1',
    'simplejson=3.17.2',
    'scipy=1.6.3',
    'sqlalchemy=1.4.11',
    'tabulate=0.8.9',
    'xlsxwriter=1.4.0'
]

channel_names = ['conda-forge', 'bioconda']
channels = [f"--channel={name}" for name in channel_names]


@nox.session(python=py_versions[0])
def flake8_diff(session):
    session.install('flake8',
                    'flake8-bugbear',
                    'flake8-builtins',
                    'flake8-comprehensions')
    session.run('sh', '-c', 'git diff -U0 --staged HEAD | flake8 --diff', external=True)


@nox.session(python=py_versions[0])
def flake8(session):
    session.install('flake8',
                    'flake8-bugbear',
                    'flake8-builtins',
                    'flake8-comprehensions')
    session.run('flake8', 'metatlas', 'tests')


@nox.session(python=py_versions[0])
def black(session):
    session.install('black')
    session.run('black', '--check', '--diff', '--color', 'metatlas', 'tests')


@nox.session(python=py_versions[0])
def pylint(session):
    session.install('pylint')
    session.run('pylint', 'metatlas', 'tests')


@nox.session(venv_backend='conda', python=py_versions, reuse_venv=True)
def unit_tests(session):
    test_deps = ['pytest=6.2.3', 'pytest-cov=2.11.1']
    session.conda_install(*(channels+run_deps+test_deps))
    session.install('--no-deps', 'dataset==1.5.0')
    session.run('pytest', '-vv', '--cov', 'metatlas', 'tests/unit/', env={'METATLAS_LOCAL': 'TRUE'})


@nox.session(python=py_versions[0])
def system_tests(session):
    test_deps = ['pytest==6.2.3']
    session.install(*test_deps)
    session.run('pytest', 'tests/system/')


@nox.session(python=py_versions[0])
def install_git_hooks(session):
    session.install('pre-commit')
    session.run('pre-commit', 'install')


@nox.session(python=py_versions[0])
def update_git_hooks(session):
    session.install('pre-commit')
    session.run('pre-commit', 'autoupdate')
