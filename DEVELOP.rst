Developing For Metatlas
=======================

Setup
#####

1. git clone https://github.com/biorack/metatlas.git
2. Install Python 3.8+
3. Install Pip
4. Install a conda environment manager. You can get miniconda `here <https://docs.conda.io/en/latest/miniconda.html>`_.
5. Install `Docker <https://docs.docker.com/get-docker/>`_.
6. Install Nox with :code:`pip install --user --upgrade nox`
7. Install git pre-commit hooks with :code:`nox -s install_git_hooks`

Testing
#######

Running :code:`nox` will run the default tests. Using :code:`nox -s <session_name>` you can access more tests.

Metatlas makes use of `pre-commit <https://pre-commit.com/>`_ to manage pre-commit git hooks. The hooks are
automatically updated when the default nox test suite runs. These hooks check that the diff of the code to be
committed passes flake8 along with a few other minor tests (see .pre-commit-config.yaml).

The current tests are in ./tests and are fairly minimal. The "system_tests" uses Papermill to run
the workflow in ./notebooks/reference/Targeted.ipynb. More high-quality unit tests are needed and the
system_tests need to be expanded to cover more functionaltiy.

Merging into Main
#################

The default nox test suite must be passing before code can be merged into the main branch on GitHub.
This is enforced by Github.
