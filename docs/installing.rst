Installation
============

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

Instructions
------------

monetio itself is a pure Python package, but some of it's dependencies may not be.
The simplest way to install MONET is to install it from the channel bbakernoaa::

    $ conda install -c bbakernoaa monetio

This will install all of the dependencies needed by monetio along with monetio itself.

If you choose to install it manually you can install the dependencies we recommend using the the following command::

    $ conda config --add channels conda-forge
    $ conda install xarray dask netCDF4 numpy pandas s3fs

We recommend using the community maintained `conda-forge <https://conda-forge.github.io/>`_ channel
if you need difficult\-to\-build dependencies such as cartopy, pynio or PseudoNetCDF::

    $ conda install -c conda-forge xarray pandas matplotlib seaborn cartopy pseudonetcdf

To install MONET currently you must install with pip.  This can be done directly
from the github page::

    $ pip install git+https://github.com/noaa-oar-arl/monetio.git

or you can manually download it from GitHub and install it using the setup.py::

    $ git clone https://github.com/noaa-oar-arl/MONET.git
    $ cd MONET
    $ pip install setup.py

.. [1] MONET plans to drop support for python 2.7 sometime in 2019. This
   means that new releases of xarray published after this date will only be
   installable on python 3+ environments, but older versions of xarray will
   always be available to python 2.7 users. For more information see the
   following references:

      - `Python 3 Statement <http://www.python3statement.org/>`__
      - `Tips on porting to Python 3 <https://docs.python.org/3/howto/pyporting.html>`__
