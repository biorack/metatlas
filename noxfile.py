import nox

py_versions = ['3.8']

nox.options.error_on_external_run = True

full_env_args = {'venv_backend': 'conda', 'reuse_venv': False, 'python': py_versions}


@nox.session(**full_env_args)
def pylint(session):
    session.conda_install('--file', 'metatlas_env.yml')
    session.run('pylint', 'metatlas')


@nox.session(**full_env_args)
def test(session):
    session.conda_install('--file', 'metatlas_env.yml')
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
