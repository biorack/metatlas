set -x

# Install a local version of the Kbase narrative
conda remove -n kbase --all -y
conda create -y -n kbase python=2.7 scipy pytables matplotlib readline
source activate kbase
pip install -e .
pip install pexpect

# based on instructions in narrative/docs/developer.md
rm -rf bootstrap
git clone https://github.com/kbase/bootstrap bootstrap
pip install -r bootstrap/kb_python_runtime/python-pip-list-narrative

rm -rf narrative
git clone https://github.com/kbase/narrative narrative
cd narrative
git submodule init
git submodule update

export VIRTUAL_ENV="kbase"
bash install.sh

# point to the CI services      
sed -i '' 's:kbase.us/services/:ci.kbase.us/services/:' src/config.json
