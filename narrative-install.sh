set -x

# Install a local version of the Kbase narrative

conda create -y -n kbase python=2.7 scipy pytables matplotlib
source activate kbase
pip install -e .

git clone https://github.com/kbase/bootstrap ..
pip install -r ../bootstrap/kb_python_runtime/python-pip-list-narrative

git clone https://github.com/kbase/narrative ..
cd ../narrative
git submodule init
git submodule update
pip install -r src/requirements.txt 
cd src
python setup.py install

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
