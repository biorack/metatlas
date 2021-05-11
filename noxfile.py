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
    'tests'
    ]


nox.options.error_on_external_run = True


@nox.session(python=py_versions[0])
def flake8_diff(session):
    session.install('flake8',
                    'flake8-bugbear',
                    'flake8-builtins',
                    'flake8-comprehensions')
    session.run('sh', '-c', 'git diff -U0 -w --staged HEAD | flake8 --diff', external=True)


@nox.session(python=py_versions[0])
def flake8_all(session):
    session.install('flake8',
                    'flake8-bugbear',
                    'flake8-builtins',
                    'flake8-comprehensions')
    session.run('flake8', 'metatlas', 'tests')


@nox.session(python=py_versions[0])
def flake8(session):
    session.install('flake8',
                    'flake8-bugbear',
                    'flake8-builtins',
                    'flake8-comprehensions')
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
    session.install('--no-deps', 'pylint==2.7.2', 'pytest==6.2.3')
    session.run('pylint', *more_checks)


@nox.session(venv_backend='conda', python=py_versions, reuse_venv=True)
def unit_tests(session):
    session.run('conda', 'env', 'update', '--prefix', session.virtualenv.location,
                '--file', 'docker/metatlas_env.yaml', silent=True)
    session.install('--no-deps', 'pytest==6.2.3', 'pytest-cov==2.11.1')
    session.run('pytest', '-vv', *session.posargs, '--cov', 'metatlas', 'tests/unit/',
                env={'METATLAS_LOCAL': 'TRUE'})


@nox.session(python=py_versions[0])
def system_tests(session):
    session.install('pytest==6.2.3')
    session.run('pytest', '-vv', *session.posargs, 'tests/system/')


@nox.session(python=py_versions[0])
def install_git_hooks(session):
    session.install('pre-commit')
    session.run('pre-commit', 'install')


@nox.session(python=py_versions[0])
def update_git_hooks(session):
    session.install('pre-commit')
    session.run('pre-commit', 'autoupdate')
