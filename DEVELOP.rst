Developing For Metatlas
=======================

Setup
#####

1. Install Python 3.8+
2. Install Pip
3. Install a conda environment manager. You can get miniconda `here <https://docs.conda.io/en/latest/miniconda.html>`_.
4. Install `Docker <https://docs.docker.com/get-docker/>`_.
5. Install Nox with :code:`pip install --user --upgrade nox`
6. :code:`git clone https://github.com/biorack/metatlas.git`
7. Install git pre-commit hooks with :code:`cd metatlas && nox -s install_git_hooks`

Local Development
#################

You can locally test within a Jupyter Notebook enviroment by running :code:`./docker/local_jupyter.sh`
and then pointing your web browser to `http://127.0.0.1:8888/ <http://127.0.0.1:8888/>`_.
This pulls in the same image used in the systems_tests that contains a small sqlite3 database at
/work/root_workspace.db. This database contains references to 4 LCMS run h5 files from experiment
20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583 along with a 6-compound atlas.
The working directory of the Jupyter notebook needs to be /work so that the database can be found.
You can either copy your notebook into /work or within your notebook do :code:`%cd /work`. Your local
copy of the metatlas git repo will be mounted at /src and $(pwd)/out on the host will be mounted at
/out in the container.

Testing
#######

Running :code:`nox` will run the default tests. Using :code:`nox -s <session_name>` you can access more tests. See
the `nox documentation <https://nox.thea.codes/>`_ for more information on running tests.

Metatlas makes use of `pre-commit <https://pre-commit.com/>`_ to manage pre-commit git hooks. The hooks are
automatically updated when the default nox test suite runs. These hooks check that the diff of the code to be
committed passes `flake8  <https://flake8.pycqa.org/>`_ along with a few other minor tests
(see .pre-commit-config.yaml). The flake8 parameters have been modified to allow for a maximum line length of
110 characters.

The current tests are in ./tests and are fairly minimal. The "system_tests" uses
`Papermill <https://papermill.readthedocs.io/>`_ to run
The workflow in `./notebooks/reference/Targeted.ipynb <blob/main/notebooks/reference/Targeted.ipynb>`_ in a
"headless", non-interactive mode.
More high-quality unit tests are needed and the system_tests need to be expanded to cover more functionaltiy.

Merging into Main
#################

The default nox test suite must be passing before code can be merged into the main branch on GitHub.
This is enforced by Github. This means you cannot commit directly to main. You'll need to make a
branch, commit to the branch, and then create a pull request.
