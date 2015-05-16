set -x

ipython_branch=1.x

# Install a local version of the Kbase narrative
conda remove -n kbase --all -y
conda create -y -n kbase python=2.7 scipy pytables matplotlib readline
source activate kbase
pip install -e .
pip install pexpect

rm -rf bootstrap
git clone https://github.com/kbase/bootstrap bootstrap
pip install -r bootstrap/kb_python_runtime/python-pip-list-narrative

rm -rf narrative
git clone https://github.com/kbase/narrative narrative
cd narrative

# install ipython from a specific git branch
git clone https://github.com/ipython/ipython.git -b ${ipython_branch}
cd ipython 
python setup.py install

git submodule init
git submodule update
cd ../src
pip install -r requirements.txt 
python setup.py install

# point to the CI services
sed -i '' 's:kbase.us/services/:ci.kbase.us/services/:' config.json

# Install IPython again (!)
cd ../ipython
python setup.py install
cd ..

# Set up the run_notebook script
# ------------------------------
echo "Installing scripts"
tgt="kbase-narrative"
i=0
while read s
do
    echo $s
    if [ $i = 0 ]
    then
        echo d=`pwd`
        echo e=$(dirname `which python`)
        i=1
    fi
done < run_notebook.tmpl > $tgt
d=$(dirname `which python`)
chmod 0755 $tgt
echo "Putting new $tgt command under $d"
/bin/mv $tgt $d
echo "Done installing scripts"

