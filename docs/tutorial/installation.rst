Installation/Requirements
=========================

Required dependencies
---------------------

- Python 3.9 (recommended)
- netcdf4
- wrf-python
- All required dependencies for `MONET <https://monet-arl.readthedocs.io/en/stable/installing.html>`__
- All required dependencies for `MONETIO <https://monetio.readthedocs.io/en/stable/installing.html>`__

Conda YAML files
~~~~~~~~~~~~~~~~
Examples of conda configuration environment.yml files that include a record
of all the dependencies are available via the GitHub:

- `NCAR Cheyenne environment.yml <https://github.com/NOAA-CSL/MELODIES-MONET/tree/develop/python_env_ymls/cheyenne>`__
- `NOAA Hera environment.yml <https://github.com/NOAA-CSL/MELODIES-MONET/tree/develop/python_env_ymls/hera>`__

General Instructions
--------------------

Setting up MELODIES MONET on HPC machines can be complicated. To help users 
get started, instructions for specific HPC machine environments are in the 
Appendix. If you are installing MELODIES MONET on NCAR Cheyenne or NOAA Hera 
follow these machine specific instructions instead.

- :ref:`NCAR Cheyenne <appendix/machine-specific-install:NCAR HPC Cheyenne/Casper>`
- :ref:`NOAA Hera <appendix/machine-specific-install:NOAA HPC Hera>`

To install MELODIES MONET on your laptop or on HPC machines in general follow 
these instructions: 
 
(a) Set up a conda environment with all the dependencies, including MONET and 
MONETIO::

    $ conda create --name melodies-monet python=3.9
    $ conda activate melodies-monet
    $ conda install -y -c conda-forge netcdf4 wrf-python jupyterlab monet monetio

(b) Clone and link the latest versions of MONET and MONETIO from GitHub to
your conda environment::

    $ git clone git@github.com:noaa-oar-arl/monet.git
    $ cd monet
    $ git checkout develop
    $ pip install --no-deps --editable .
    
    $ git clone git@github.com:noaa-oar-arl/monetio.git
    $ cd monetio
    $ git checkout develop
    $ pip install --no-deps --editable .

\(c) Clone the MELODIES MONET package::

    $ git clone git@github.com:NOAA-CSL/MELODIES-MONET.git
    
**Note to developers:** In order to incorporate updates to MELODIES MONET, you 
will need to fork the repository to your own GitHub account, make changes, and 
submit a pull request. For details, see 
:ref:`develop/developers_guide:How to incorporate updates to MELODIES MONET`.


