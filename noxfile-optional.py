import nox

py_versions = ['3.8']

nox.options.error_on_external_run = True

single_py_args = {'venv_backend': 'conda', 'reuse_venv': False}
multi_py_args = {**single_py_args, **{'python': py_versions}}


@nox.session(reuse_venv=True)
def flake8(session):
    session.install('flake8',
                    'flake8-bugbear',
                    'flake8-builtins',
                    'flake8-comprehensions')
    session.run('flake8', 'metatlas', 'tests')


@nox.session(reuse_venv=True)
def black(session):
    session.install('black')
    session.run('black', '--check', '--diff', '--color', 'metatlas', 'tests')


@nox.session(**multi_py_args)
def pylint(session):
    session.conda_install('--channel=conda-forge', 'pylint')
    session.run('pylint', 'metatlas')


@nox.session(**multi_py_args)
def test(session):
    session.conda_install('--channel=conda-forge', 'pytest', 'pytest-cov')
    session.run('pytest', '--cov', 'metatlas')


@nox.session(reuse_venv=True)
def flake8_diff(session):
    session.install('flake8',
                    'flake8-bugbear',
                    'flake8-builtins',
                    'flake8-comprehensions')
    session.run('sh', '-c', 'git diff -U0 --staged HEAD | flake8 --diff', external=True)


@nox.session(reuse_venv=True)
def install_git_hooks(session):
    session.install('pre-commit')
    session.run('pre-commit', 'install')
