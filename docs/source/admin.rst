Administration
--------------

NERSC
-----

The project folder at NERSC is ``/project/projectdirs/metatlas/``.

Raw data is stored in ``/project/projectdirs/metatlas/raw_data``.

Log in to NERSC via ``ssh username@corigrid.nersc.gov``.

There is a pseudo-super user account called ``pasteur``. To log in to
pasteur:

::

    module load globus
    myproxy-logon -s nerscca.nersc.gov
    NERSC PW
    gsissh localhost -p 2222 -l pasteur

MySQL
-----

We have a MySql database at NERSC.

The remote path (from python) is
``mysql+pymysql://meta_atlas_admin:<password>@scidb1.nersc.gov/meta_atlas``

The local path (at NERSC) is:
``mysql meta_atlas -u meta_atlas_admin -h scidb1.nersc.gov -p``.


You can use the database on a local machine or a machine not at NERC. Here are the steps

Clone the database on NERSC
::

    mysqldump meta_atlas -u meta_atlas_admin -h scidb1.nersc.gov  -p > dump.sql
    MYSQL PW

Bring the databsae to local machine

::

    scp username@corigrid.nersc.gov:dump.sql .


Add the database to your mysql server
::

    mysql -uroot -p
    LOCAL MYSQL ADMIN PW
    CREATE DATABASE metatlas;
    USE metatlas;
    SOURCE dump.sql;


Create a user
::

    CREATE USER new_user IDENTIFIED BY 'some_password';
    GRANT ALL ON my_project_copy.* TO 'new_user'@'localhost' IDENTIFIED BY 'some_password';
    FLUSH PRIVILEGES;

Alternatively, you can create a user using myphpadmin and give him/her full access to the database.

Add in index and specify type for inchi keys
::

    ALTER TABLE compounds MODIFY inchi_key VARCHAR(27);
    ALTER TABLE `compounds` ADD INDEX `inchi_key` (`inchi_key`);

CRON
----

We have a CRON job set up to run metatlas directory watcher every ten
minutes. The file that is run is
``/project/projectdirs/metatlas/run_dirwatch2.sh``. The CRON job
settings are in
``/project/projectdirs/metatlas/crontab.corigrid.pasteur``. You have to
log into pasteur to modify the CRON job.

Anaconda
--------

We have an Anaconda instance at NERSC:
``/project/projectdirs/metatlas/anaconda``. We have all the required
libraries stored there. See the ``installation`` instructions for how to
users are meant to use it. You can ``pip install`` and
``conda install``. The packages required for ``metatlas`` are listed the
``setup.py``, except for qgrid, which was installed directly from the
QGrid repo as ``pip install https://github.com/quantopian/qgrid/``.

R packages
----------

We are keeping R packages in
``/global/project/projectdirs/metatlas/r_pkgs/``, which are
automatically available because we set the ``R_LIBS`` variable in
``metatlas/__init__.py``. To load new packages, log in to NERSC, and
follow this example (replacing the source and lib name with your own):

::

    $ module load R
    $ R
    > source("http://bioconductor.org/biocLite.R")
    > # accept installation into a local directory (copy that directory)
    > biocLite("my_lib")
    > quit()
    $ cp -r ~/R/x86_64-unknown-linux-gnu-library/3.1/my_lib /global/project/projectdirs/metatlas/r_pkgs/my_lib

The currently installed packages are:

::

    biocLite("mzR")
    biocLite("xcms")



RDKIT Pacakge on Debian/Ubuntu
------------------------------
::

    conda install -c https://conda.anaconda.org/rdkit rdkit


Synchronizing the Repo
----------------------

There is a ``make deploy`` target in the top level ``Makefile`` in the
repo that will publish the latest changes to documentation and ``rsync``
the files to the anaconda environment at NERSC.


Administration
--------------

NERSC
=====

The project folder at NERSC is ``/project/projectdirs/metatlas/``.

Raw data is stored in ``/project/projectdirs/metatlas/raw_data``.

Log in to NERSC via ``ssh username@corigrid.nersc.gov``.

There is a pseudo-super user account called ``pasteur``. To log in to
pasteur:

::

    module load globus
    myproxy-logon -s nerscca.nersc.gov
    NERSC PW
    gsissh localhost -p 2222 -l pasteur

MySQL
=====

We have a MySql database at NERSC.

The remote path (from python) is
``mysql+pymysql://meta_atlas_admin:<password>@scidb1.nersc.gov/meta_atlas``

The local path (at NERSC) is:
``mysql meta_atlas -u meta_atlas_admin -h scidb1.nersc.gov -p``.


You can use the database on a local machine or a machine not at NERC. Here are the steps

Clone the database on NERSC
::

    mysqldump meta_atlas -u meta_atlas_admin -h scidb1.nersc.gov  -p > dump.sql
    MYSQL PW

Bring the databsae to local machine

::

    scp username@corigrid.nersc.gov:dump.sql .


Add the database to your mysql server
::

    mysql -uroot -p
    LOCAL MYSQL ADMIN PW
    CREATE DATABASE metatlas;
    USE metatlas;
    SOURCE dump.sql;


Create a user
::

    CREATE USER new_user IDENTIFIED BY 'some_password';
    GRANT ALL ON my_project_copy.* TO 'new_user'@'localhost' IDENTIFIED BY 'some_password';
    FLUSH PRIVILEGES;

Alternatively, you can create a user using myphpadmin and give him/her full access to the database.


CRON
====

We have a CRON job set up to run metatlas directory watcher every ten
minutes. The file that is run is
``/project/projectdirs/metatlas/run_dirwatch2.sh``. The CRON job
settings are in
``/project/projectdirs/metatlas/crontab.corigrid.pasteur``. You have to
log into pasteur to modify the CRON job.

Anaconda
========

We have an Anaconda instance at NERSC:
``/project/projectdirs/metatlas/anaconda``. We have all the required
libraries stored there. See the ``installation`` instructions for how to
users are meant to use it. You can ``pip install`` and
``conda install``. The packages required for ``metatlas`` are listed the
``setup.py``, except for 

qgrid, which was installed directly from the
QGrid repo as ``pip install https://github.com/quantopian/qgrid/``

and

RDKit, which was installed from the rdkit conda build as
``conda install -c https://conda.anaconda.org/rdkit rdkit``

R packages
==========

We are keeping R packages in
``/global/project/projectdirs/metatlas/r_pkgs/``, which are
automatically available because we set the ``R_LIBS`` variable in
``metatlas/__init__.py``. To load new packages, log in to NERSC, and
follow this example (replacing the source and lib name with your own):

::

    $ module load R
    $ R
    > source("http://bioconductor.org/biocLite.R")
    > # accept installation into a local directory (copy that directory)
    > biocLite("my_lib")
    > quit()
    $ cp -r ~/R/x86_64-unknown-linux-gnu-library/3.1/my_lib /global/project/projectdirs/metatlas/r_pkgs/my_lib

The currently installed packages are:

::

    biocLite("mzR")
    biocLite("xcms")



RDKIT Pacakge on Debian/Ubuntu
==============================
::

    conda install -c https://conda.anaconda.org/rdkit rdkit


Synchronizing the Repo
======================

There is a ``make deploy`` target in the top level ``Makefile`` in the
repo that will publish the latest changes to documentation and ``rsync``
the files to the anaconda environment at NERSC.


Synchronizing the Repo Between Github and NERSC
===============================================

Work in progress...
Also detailed here for reference:
https://github.com/biorack/metatlas/blob/master/update_nersc_packages.sh

