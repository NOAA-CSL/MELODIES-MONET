Installation/Requirements
=========================

Required dependencies
---------------------

- Python 3.9 (recommended)
- ``pyyaml`` (to read control files)
- ``monet``, which brings `many dependencies <https://monet-arl.readthedocs.io/en/stable/installing.html>`__
- ``monetio``, which brings `a few dependencies <https://monetio.readthedocs.io/en/stable/installing.html>`__

Conda YAML files
~~~~~~~~~~~~~~~~

Examples of `conda <https://conda.io>`__
`environment.yml files <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>`__
that include a record
of all the dependencies are available via the GitHub:

- `NCAR Cheyenne environment.yml <https://github.com/NOAA-CSL/MELODIES-MONET/tree/develop/python_env_ymls/cheyenne>`__
- `NOAA Hera environment.yml <https://github.com/NOAA-CSL/MELODIES-MONET/tree/develop/python_env_ymls/hera>`__

Optional dependencies
---------------------

- ``netcdf4`` (`from Unidata <https://unidata.github.io/netcdf4-python/>`__; most likely needed for reading model/obs datasets)
- ``wrf-python`` (needed in order to use the WRF-Chem reader)
- ``click`` (to use the :doc:`/cli`)
- ``pooch`` (to enable automatic downloading of :doc:`tutorial datasets </examples/tutorial-data>`)

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


General Instructions (old)
--------------------------

Setting up MELODIES MONET on HPC machines can be complicated. To help users 
get started, instructions for specific HPC machine environments are in the 
Appendix. If you are installing MELODIES MONET on NCAR Cheyenne or NOAA Hera 
follow these machine specific instructions instead.

- :ref:`NCAR Cheyenne <appendix/machine-specific-install:NCAR HPC Cheyenne/Casper>`
- :ref:`NOAA Hera <appendix/machine-specific-install:NOAA HPC Hera>`

.. important::
   The instructions below are for cloning a repository using SSH.
   If you prefer, you can also clone the monet, monetio, and
   MELODIES-MONET repositories using HTTPS [#clone]_.

To install MELODIES MONET on your laptop or on HPC machines in general follow 
these instructions: 
 
(a) Set up a conda environment with all the dependencies, including MONET and 
MONETIO::

    $ conda create --name melodies-monet python=3.9
    $ conda activate melodies-monet
    $ conda install -y -c conda-forge netcdf4 wrf-python jupyterlab monet monetio

(b) Clone [#clone]_ and link the latest versions of MONET and MONETIO from GitHub to
your conda environment::

    $ git clone git@github.com:noaa-oar-arl/monet.git
    $ cd monet
    $ git checkout develop
    $ pip install --no-deps --editable .
    
    $ git clone git@github.com:noaa-oar-arl/monetio.git
    $ cd monetio
    $ git checkout develop
    $ pip install --no-deps --editable .

\(c) Clone [#clone]_ the MELODIES MONET package::

    $ git clone git@github.com:NOAA-CSL/MELODIES-MONET.git
    
**Note to developers:** In order to incorporate updates to MELODIES MONET, you 
will need to fork the repository to your own GitHub account, make changes, and 
submit a pull request. For details, see 
:ref:`develop/developers_guide:How to incorporate updates to MELODIES MONET`.


.. _clone-notes:
.. [#clone] Note that in order to do an SSH clone,
   e.g. ::

      $ git clone git@github.com:noaa-oar-arl/monet.git

   you must have already
   `added an SSH key <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`__
   to your GitHub account for your current machine.
   If you are new to GitHub, check out
   `this GitHub tutorial <https://jlord.us/git-it/>`__.
   We recommend the SSH method, but if you don't add an SSH key
   you can still clone the repositories via HTTPS, e.g. ::

       $ git clone https://github.com/noaa-oar-arl/monet.git
