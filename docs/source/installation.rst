Installation
************************

Log in to https://ipython.nersc.gov/.
Run the following in a notebook cell to get set up (note this only has to be done once for a given user):

```python
%%file ~/.ipython/profile_default/startup/ipython_config.py
import sys
sys.path.insert(0, '/project/projectdirs/metatlas/python_pkgs/')
import os
os.environ['R_LIBS_USER'] = '/project/projectdirs/metatlas/r_pkgs/'
```

You can then import and use all of `metatlas`.
Workspace objects will persist between notebooks in a database file.

You can press SHIFT+TAB to get documentation while accessing the functions, including 
between arguments to the function.

```python
from metatlas import metatlas_objects
e = metatlas_objects.Experiment(name='test6')
c = metatlas_objects.Compound(name='hiya')
a = metatlas_objects.Atlas(name='what', compounds=[c])
e.atlases.append(a)
e.save()
b = metatlas_objects.get_experiment('test6')
b.atlases[0].compounds
```
