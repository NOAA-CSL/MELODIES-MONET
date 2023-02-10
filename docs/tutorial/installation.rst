Installation/Requirements
=========================

Required dependencies [#yaml]_
------------------------------

- Python 3.6+ (3.9 recommended)
- ``pyyaml`` (to read control files)
- ``monet``, which brings `many dependencies <https://monet-arl.readthedocs.io/en/stable/installing.html>`__
- ``monetio``, which brings `a few dependencies <https://monetio.readthedocs.io/en/stable/installing.html>`__

Optional dependencies
---------------------

- ``netcdf4`` (`from Unidata <https://unidata.github.io/netcdf4-python/>`__; most likely needed for reading model/obs datasets)
- ``wrf-python`` (needed in order to use the WRF-Chem reader)
- ``click`` (to use the :doc:`/cli`)
- ``pooch`` (to enable automatic downloading of :doc:`tutorial datasets </examples/tutorial-data>`)

.. _user-install-instructions:

General instructions
--------------------

If you are a user and are not planning to modify MELODIES MONET itself,
installing it is relatively simple.
The stable branch of MELODIES MONET (``main``) should always be compatible with the
conda-forge releases of MONET/MONETIO.
First create and activate a conda environment::

    $ conda create --name melodies-monet python=3.9
    $ conda activate melodies-monet

Add dependencies from conda-forge::

    $ conda install -y -c conda-forge pyyaml monet monetio netcdf4 wrf-python click pooch

Now, install the stable branch of MELODIES MONET to the environment::

    $ pip install --no-deps https://github.com/NOAA-CSL/MELODIES-MONET/archive/main.zip


.. note::
   If you are interested in modifying what MELODIES MONET can do,
   take a look at the :doc:`/develop/developers_guide`.


.. [#yaml] Examples of `conda <https://conda.io>`__
   `environment.yml files <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>`__
   that include a record
   of all the dependencies (direct and indirect) are available via the GitHub:

   - `NCAR Cheyenne environment.yml <https://github.com/NOAA-CSL/MELODIES-MONET/tree/develop/python_env_ymls/cheyenne>`__
   - `NOAA Hera environment.yml <https://github.com/NOAA-CSL/MELODIES-MONET/tree/develop/python_env_ymls/hera>`__
