# pylint: disable=missing-function-docstring

import nox

py_versions = ['3.8', '3.9']

nox.options.sessions = [
        "flake8_diff",
        "flake8",
        "black",
        "pylint-3.8",
        "unit_tests-3.8",
        "system_tests-3.8",
        "update_git_hooks",
        ]

# files we can run all the checks on, as they don't contain legacy code that
# has not yet been updated to pass all checks.
more_checks = [
    'metatlas/datastructures/metatlas_dataset.py',
    'metatlas/tools/logging.py',
    'tests'
    ]

notebooks = [
    'notebooks/reference/Targeted.ipynb',
    ]

pytest_deps = [
        'attrs==21.2.0',
        'coverage==5.5',
        'iniconfig==1.1.1',
        'packaging==20.9',
        'pluggy==0.13.1',
        'py==1.10.0',
        'pyparsing==2.4.7',
        'pytest==6.2.4',
        'pytest-cov==2.11.1',
        'pytest-mock==3.6.1',
        'toml==0.10.2',
        ]

nbqa_deps = [
        'nbqa==0.8.1',
        'tokenize-rt==4.1.0',
        'importlib-metadata==4.0.1',
        'astroid==2.5.6',
        'wrapt==1.12.1',
        'lazy_object_proxy==1.6.0',
        'isort==5.8.0',
        ]

flake8_deps = [
    'flake8',
    'flake8-bugbear',
    'flake8-builtins',
    'flake8-comprehensions',
    ]

nox.options.error_on_external_run = True


@nox.session(python=py_versions[0])
def flake8_diff(session):
    session.install(*flake8_deps)
    session.run('sh', '-c', 'git diff -U0 -w --staged HEAD | flake8 --diff', external=True)


@nox.session(python=py_versions[0])
def flake8_all(session):
    session.install(*flake8_deps)
    session.run('flake8', 'metatlas', 'tests')


@nox.session(python=py_versions[0])
def flake8(session):
    session.install(*flake8_deps)
    session.run('flake8', *more_checks)


@nox.session(python=py_versions[0])
def black_all(session):
    session.install('black')
    session.run('black', '--check', '--diff', '--color', 'metatlas', 'tests')


@nox.session(python=py_versions[0])
def black(session):
    session.install('black')
    session.run('black', '--check', '--diff', '--color', *more_checks)


@nox.session(python=py_versions[0])
def blacken(session):
    """this modifies the files to meet black's requirements"""
    session.install('black')
    session.run('black', *more_checks)


@nox.session(venv_backend='conda', python=py_versions, reuse_venv=True)
def pylint(session):
    session.run('conda', 'env', 'update', '--prefix', session.virtualenv.location,
                '--file', 'docker/metatlas_env.yaml', silent=True)
    session.install('--no-deps', *pytest_deps)
    session.run('pylint', *more_checks)


@nox.session(venv_backend='conda', python=py_versions, reuse_venv=True)
def pylint_nb(session):
    session.run('conda', 'env', 'update', '--prefix', session.virtualenv.location,
                '--file', 'docker/metatlas_env.yaml', silent=True)
    session.install('--no-deps', *nbqa_deps, 'pylint')
    session.run('nbqa', 'pylint', *notebooks)


@nox.session(python=py_versions[0])
def flake8_nb(session):
    session.install(*nbqa_deps, *flake8_deps)
    session.run('nbqa', 'flake8', *notebooks)


@nox.session(python=py_versions[0])
def black_nb(session):
    session.install('black', *nbqa_deps)
    session.run('nbqa', 'black', '--check', *notebooks)


@nox.session(python=py_versions[0])
def blacken_nb(session):
    """this modifies notebook files to meet black's requirements"""
    session.install('black', *nbqa_deps)
    session.run('nbqa', 'black', '--nbqa-mutate', *notebooks)


@nox.session(venv_backend='conda', python=py_versions, reuse_venv=True)
def unit_tests(session):
    session.run('conda', 'env', 'update', '--prefix', session.virtualenv.location,
                '--file', 'docker/metatlas_env.yaml', silent=True)
    session.install('--no-deps', *pytest_deps)
    session.run('pytest', '-vv', *session.posargs, '--cov', 'metatlas', 'tests/unit/',
                env={'METATLAS_LOCAL': 'TRUE'})


@nox.session(python=py_versions[0])
def system_tests(session):
    session.install('--no-deps', *pytest_deps)
    session.run('pytest', '-vv', *session.posargs, 'tests/system/')


@nox.session(python=py_versions[0])
def install_git_hooks(session):
    session.install('pre-commit')
    session.run('pre-commit', 'install')


@nox.session(python=py_versions[0])
def update_git_hooks(session):
    session.install('pre-commit')
    session.run('pre-commit', 'autoupdate')
