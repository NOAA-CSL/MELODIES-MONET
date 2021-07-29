Installation/Requirements
=========================

Required dependencies
---------------------

- Python 3.6+
- `numpy <http://www.numpy.org/>`__ (1.11 or later)
- `pandas <http://pandas.pydata.org/>`__ (0.18.0 or later)
- `xarray <http://xarray.pydata.org/>`__ (0.10 or later)
- `dask <http://dask.pydata.org/>`__
- `netcdf4 <http://unidata.github.io/netcdf4-python/>`__
- `s3fs <https://github.com/dask/s3fs>`__

For parallel computing
~~~~~~~~~~~~~~~~~~~~~~

- `dask.array <http://dask.pydata.org>`__ (0.9.0 or later): required for

Conda YAML files
~~~~~~~~~~~~~~~~
Examples of conda configuration yaml files that include a record of all the dependencies are available via the GitHub:

- `MELODIES-MONET conda yml <https://github.com/NOAA-CSL/MELODIES-MONET/tree/develop_plots/monet_analysis/data/python_env_ymls>`__

General Instructions
--------------------

MELODIES-MONET itself is a pure Python package, but some of it's dependencies may not be.
First, set up a conda environment with all the dependencies, including MONET and MONETIO::

    $ conda create --name monet_py36 python=3.6
    $ conda activate monet_py36
    $ conda install netcdf4
    $ conda install -y -c conda-forge wrf-python
    $ conda install -y -c conda-forge jupyter
    $ conda install -y -c conda-forge monet
    $ conda install -y -c conda-forge monetio

Then, fork the MELODIES-MONET package to your own GitHub and clone the repository to the location you will be using the plotting package.


.. [1] For instructions for specific machine environments see Appendix:

      - `NCAR cheyenne <../tutorial/machine-specific-install.html#NCAR-HPC-cheyenne>`__

.. [2] For more information see the following references:

      - `Python 3 Statement <http://www.python3statement.org/>`__
      - `Tips on porting to Python 3 <https://docs.python.org/3/howto/pyporting.html>`__


