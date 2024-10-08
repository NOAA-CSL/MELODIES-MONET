Installation/Requirements
=========================

Required dependencies [#yaml]_
------------------------------

- Python 3.6+ (3.9 -- 3.11 recommended)
- ``pandas`` 1 (most of the system works with pandas 2 as well, but some of the :doc:`/cli` commands for downloading observational datasets using MONETIO might fail)
- ``pyyaml`` (to read control files)
- ``monet``, which brings `many dependencies <https://monet-arl.readthedocs.io/en/stable/installing.html>`__
- ``monetio``, which brings `a few dependencies <https://monetio.readthedocs.io/en/stable/installing.html>`__

Optional dependencies
---------------------

- ``netcdf4`` (`from Unidata <https://unidata.github.io/netcdf4-python/>`__; most likely needed for reading model/obs datasets)
- ``wrf-python`` (needed in order to use the WRF-Chem reader, currently the conda package seems to require Python < 3.12)
- ``typer`` (to use the :doc:`/cli`;
  add ``rich`` `for <https://typer.tiangolo.com/release-notes/#060>`__ fancy tracebacks and ``--help``)
- ``pooch`` (to enable automatic downloading of :doc:`tutorial datasets </examples/tutorial-data>`)

Incompatibilities
-----------------
- pandas=1 is incompatible with matplotlib 3.9+.
- wrf-python, at least in the official conda-forge package, is not available for Python 3.12+, until `this build issue <https://github.com/conda-forge/wrf-python-feedstock/pull/70>`__ is resolved.

.. _user-install-instructions:

General instructions
--------------------

If you are a user and are not planning to modify MELODIES MONET itself,
installing it is relatively simple. There are two methods available.

Option 1) Using Conda
^^^^^^^^^^^^^^^^^^^^^
We have recently created a conda-forge release of MELODIES MONET to make installation very simple 
with just 1 line of code below::

    $ conda create --name melodies-monet -y -c conda-forge python=3.9 melodies-monet wrf-python jupyterlab

.. note::
   Currently, the wrf-python conda package is not compatible with Apple Silicon (Apple machines using Intel should be fine). If you need to run the WRF-Chem reader and only have access to a machine using Apple Silicon, you can try compiling it from source code from the official repos.

Option 2) Using Conda and GitHub
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You are also welcome to install using our old method. The stable branch of MELODIES MONET (``main``) 
should always be the same as the conda-forge release of MELODIES MONET and be compatible with the
conda-forge releases of MONET/MONETIO.

First create and activate a conda environment::

    $ conda create --name melodies-monet python=3.9
    $ conda activate melodies-monet

Add dependencies from conda-forge::

    $ conda install -y -c conda-forge pyyaml pandas=1 'matplotlib-base<3.9' monet monetio netcdf4 wrf-python typer rich pooch jupyterlab
   
Now, install the stable branch of MELODIES MONET to the environment::

    $ pip install --no-deps https://github.com/NOAA-CSL/MELODIES-MONET/archive/main.zip


.. note::
   If you are interested in modifying what MELODIES MONET can do,
   take a look at the :doc:`/develop/developers_guide`.


.. [#yaml] Examples of `conda <https://conda.io>`__
   `environment.yml files <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>`__
   that include a record
   of all the dependencies (direct and indirect) are available via the GitHub:

   - `NCAR Casper environment.yml <https://github.com/NOAA-CSL/MELODIES-MONET/tree/develop/python_env_ymls/casper>`__
   - `NOAA Hera environment.yml <https://github.com/NOAA-CSL/MELODIES-MONET/tree/develop/python_env_ymls/hera>`__
