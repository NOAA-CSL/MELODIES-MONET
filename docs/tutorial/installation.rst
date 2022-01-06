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
Examples of conda configuration environment.yaml files that include a record 
of all the dependencies are available via the GitHub:

- `NCAR Cheyenne environment.yaml <https://github.com/NOAA-CSL/MELODIES-MONET/tree/develop/python_env_ymls/cheyenne>`__
- `NOAA Hera environment.yaml <https://github.com/NOAA-CSL/MELODIES-MONET/tree/develop/python_env_ymls/hera>`__

General Instructions
--------------------

Setting up MELODIES MONET on HPC machines can be complicated. To help users 
get started, instructions for specific HPC machine environments are in the 
Appendix. If you are installing MELODIES MONET on NCAR Cheyenne or NOAA Hera 
follow these machine specific instructions instead.

- `NCAR Cheyenne <../appendix/machine-specific-install.html#NCAR-HPC-cheyenne>`__
- `NOAA Hera <../appendix/machine-specific-install.html#NOAA-HPC-hera>`__

To install MELODIES MONET on your laptop or on HPC machines in general follow 
these instructions: 
 
(a) Set up a conda environment with all the dependencies, including MONET and 
MONETIO::

    $ conda create --name monet_py36 python=3.6
    $ conda activate monet_py36
    $ conda install netcdf4
    $ conda install -y -c conda-forge wrf-python
    $ conda install -y -c conda-forge jupyter
    $ conda install -y -c conda-forge monet
    $ conda install -y -c conda-forge monetio

(b) Clone and link the latest versions of MONET and MONETIO from github to 
your conda environment::

    $ git clone git@github.com:noaa-oar-arl/monet.git
    $ cd monet
    $ git checkout develop
    $ pip install -e .
    
    $ git clone git@github.com:noaa-oar-arl/monetio.git
    $ cd monetio
    $ git checkout development
    $ pip install -e .

\(c) Clone the MELODIES MONET package::

    $ git clone git@github.com:NOAA-CSL/MELODIES-MONET.git
    
**Note to developers:** In order to incorporate updates to MELODIES MONET, you 
will need to fork the repository to your own Github account, make changes, and 
submit a pull request. For details, see 
`How to incorporate updates to MELODIES MONET <../develop/developers_guide.html#How to incorporate updates to MELODIES MONET>`__


