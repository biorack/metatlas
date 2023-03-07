Developing For Metatlas
=======================

Setup
#####

1. Install Python 3.8+ (`pyenv <https://github.com/pyenv/pyenv>`_ and `pyenv intstaller <https://github.com/pyenv/pyenv-installer>`_ can help here)
2. Install `Pip <https://pip.pypa.io/en/stable/installing/>`_
3. Install `Docker <https://docs.docker.com/get-docker/>`_.
4. Install Nox with :code:`pip install --user --upgrade nox`
5. Install Flake8 with :code:`pip install --user --upgrade flake8`
6. :code:`git clone https://github.com/biorack/metatlas.git`
7. Install git pre-commit hooks with :code:`cd metatlas && nox -s install_git_hooks`

Local Development
#################

You can locally test within a Jupyter Notebook enviroment by running :code:`./docker/local_jupyter.sh`
and then pointing your web browser to `http://127.0.0.1:8888/ <http://127.0.0.1:8888/>`_.
This pulls in the same image used in the systems_tests that contains a small sqlite3 database at
/work/root_workspace.db. This database contains references to 16 LCMS run h5 files from experiment
20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583 along with a 6-compound atlas.
The working directory of the Jupyter notebook needs to be /work so that the database can be found.
You can either copy your notebook into :code:`/work` or within your notebook do :code:`%cd /work`. Your local
copy of the metatlas git repo will be mounted at :code:`/src` and :code:`$(pwd)/out` on the host will be mounted at
:code:`/out` in the container.


The workflow in `./notebooks/reference/Targeted.ipynb </notebooks/reference/Targeted.ipynb>`_ can
easily be run locally through your web brower. In the second code block, which should be empty,
you'll want to put the following:

.. code-block:: python

  experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583"
  source_atlas_name="HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_root0"
  source_atlas_unique_id="4b05837a53494dd8b680e6b5059e1934"
  project_directory="/out"
  config_file_name="/src/test_config.yaml"
  workflow_name="Test-HILIC"
  analysis_name="EMA-POS"
  rt_alignment_number=0
  analysis_number=0
  exclude_lcmsruns['always'] = [
    '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_53_Cone-S1_5_Rg70to1050-CE102040-QlobataAkingi-S1_Run187',
    '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_54_Cone-S1_6_Rg70to1050-CE102040-QlobataAkingi-S1_Run221',
  ]

These parameters correspond to the system level test contained in :code:`./tests/system/test_targeted.py`.

You may also want to include a line for :code:`max_cpus` if your local docker configuration
makes fewer than 4 cores available.

Testing
#######

Running :code:`nox` will run the default tests. Using :code:`nox -s <session_name>` you can access more tests. See
the `nox documentation <https://nox.thea.codes/>`_ for more information on running tests.

Metatlas makes use of `pre-commit <https://pre-commit.com/>`_ to manage pre-commit git hooks. The hooks are
automatically updated when the default nox test suite runs. These hooks check that the diff of the code to be
committed passes `flake8  <https://flake8.pycqa.org/>`_ along with a few other minor tests
(see .pre-commit-config.yaml). The flake8 parameters have been modified to allow for a maximum line length of
110 characters.

The current tests are in ./tests and while the tests are substantial, there is still a lot of code
not being tested. The "system_tests" uses
`Papermill <https://papermill.readthedocs.io/>`_ to run
The workflow in `./notebooks/reference/Targeted.ipynb </notebooks/reference/Targeted.ipynb>`_ in a
"headless", non-interactive mode.

Merging into Main
#################

The unit and systems tests must be passing before code can be merged into the main branch on GitHub.
This is enforced by Github. This means you cannot commit directly to main. You'll need to make a
branch, commit to the branch, and then create a pull request.


Deployment
##########

Sucessuful pull requests to the main branch will build a new shifter image and push it to dockerhub.
A cronjob runs on cori21 under user pasteur every 5 minutes which runs shifterimg to import the image
to cori.
Another cronjob on cori21 under user pasteur pulls this git repo to
:code:`/global/common/software/m2650/metatlas-repo`.
An scrontab job runs on perlmutter under user msdata every 5 mintues which runs shifterimg to
import the image to perlmutter.
Another scrontab job on perlmutter under user msdata also pulls this git repo to
:code:`/global/common/software/m2650/metatlas-repo`.
