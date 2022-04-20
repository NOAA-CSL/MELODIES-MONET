Machine-specific Install
========================

NCAR HPC Cheyenne/Casper
------------------------

Below is a specific recipe for how to set up all the necessary Python dependencies 
on the NCAR HPC, Cheyenne/Casper. Note: these are developer-specific instructions. 
There are three steps to complete before you can use and develop MELODIES MONET 
on Cheyenne/Casper: **Step 1.** Installing the conda package manager 
Miniconda, **Step 2.** Creating a conda environment with all required dependencies, 
and **Step 3.** Cloning MELODIES MONET code.

We will use the conda package manager system to create a contained Python environment 
for running and developing MELODIES MONET. This creates a “stand-alone” instance 
of interdependent packages that will not interfere with your access to the main 
installation of Python on the system. These instructions should work for a bash 
or tcsh shell.

**Step 0:** Log into the NCAR HPC

**Step 1 Miniconda:** The conda package manager is not installed on the Cheyenne 
or Casper clusters, but you can install Miniconda in your own user space 
(`CISL instructions <https://arc.ucar.edu/knowledge_base/83853599>`_):

(a) Download latest version::

    $ curl -o install_miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

(b) Verify download::

    $ sha256sum install_miniconda.sh

    Output (June 2021)::

    $ 1314b90489f154602fd794accfc90446111514a5a72fe1f71ab83e07de9504a7

(c) Install Miniconda::

    $ bash install_miniconda.sh

    If you already have miniconda and want to try and update::

    $ bash install_miniconda.sh -u

Accept the license and default install paths. Say yes to running ``conda init``. 
Log out of Cheyenne/Casper and log back in. You may see ``(base)`` as the start 
of your terminal prompt. If you have any issues with package conflicts you may 
need to uninstall miniconda and reinstall.
If you do not see ``(base)`` in your prompt, type ``bash`` and it should appear.

**Step 2 Dependent Python Packages:** Set up a conda environment with required 
dependencies.

(a) Set up and activate CONDA environment specific to MELODIES MONET. You can 
    call this environment whatever you like, we suggest 'melodies-monet'::

    $ conda create --name melodies-monet python=3.9
    $ conda activate melodies-monet

    You should see ``(melodies-monet)`` at the start of your terminal prompt.

(b) Install the following packages to the environment. Note they have sub-packages 
    that will be downloaded. The '-y' means you will not have to interactively
    choose 'y' to proceed with downloading packages::

    $ conda install -y -c conda-forge netcdf4 wrf-python jupyterlab monet monetio

    Some main packages that are downloaded with the monet install: *cartopy, 
    cython, dask, markupsafe, matplotlib-base, pandas, pydecorate, pyresample, 
    python-stratify, seaborn, xarray, xesmf*. Some main packages that are 
    downloaded with the jupyter install: *icu, ipykernel, ipython*. If other 
    systems have timeout errors, downloading these separately may improve the 
    setup performance.

(c) Once you set up the correct packages through conda, and while your conda 
    environment is still activated, get the most recent branches of MONET and 
    MONETIO using GitHub, and link them with conda. This is done because MONET
    and MONETIO are still in active development. Create a 'monet-base' folder
    (e.g. in your ``/glade/work/<user>`` directory, accessible on both Cheyenne
    and Casper).

    Set up and link MONET within monet-base [#clone]_ ::

    $ git clone git@github.com:noaa-oar-arl/monet.git
    $ cd monet
    $ git checkout develop
    $ pip install --force-reinstall --no-deps --editable .

    Set up and link MONETIO within monet-base::

    $ git clone git@github.com:noaa-oar-arl/monetio.git
    $ cd monetio
    $ git checkout develop
    $ pip install --force-reinstall --no-deps --editable .

**Step 3: Clone the MELODIES-MONET GitHub repository** [#clone]_ ::

    $ git clone git@github.com:NOAA-CSL/MELODIES-MONET.git

**You are ready to start developing MELODIES MONET!**

.. note::
   In the NCAR JupyterHub, in the Launcher, you should see an option
   resembling 'Python [conda env:melodies-monet]'.

NOAA HPC Hera
-------------

Below is a specific recipe for how to set up all the necessary Python 
dependencies on the NOAA Hera machine. There are three steps to complete 
before you can use and develop MELODIES MONET on hera: **Step 1.** Install 
the conda package manager Anaconda/Miniconda, **Step 2.** Install MELODIES MONET,
and **Step 3.** link cartopy files.

We will use the conda package manager system to create a contained Python 
environment for running and developing MELODIES MONET. 

#. **Install Anaconda/Miniconda:** Follow the instructions
   `on the RDHPCS wiki <https://rdhpcs-common-docs.rdhpcs.noaa.gov/wiki/index.php/Anaconda>`__
   to download and run Anaconda/Miniconda on Hera. Tips for success:

   * You will need a NOAA HPC account to access the RDHPCS wiki link above.

   * Both Anaconda/Miniconda will work well for MELODIES MONET. See
     `conda instructions <https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html#anaconda-or-miniconda>`__
     to determine, which is the best option for you.
     
   * Pick a directory for your download and run the following wget command with 
     modifications if needed: ::
     
     $ wget -nH -m -nd https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

   * Follow the instructions on `conda's website <https://conda.io/projects/conda/en/latest/user-guide/install/linux.html>`__,
     which is generally described below: ::
     
     $ bash Miniconda3-latest-Linux-x86_64.sh
     
     * Follow all prompts. Installing Anaconda/Miniconda on scratch is recommended 
       due to the limited space available on your home directory. Make sure you 
       change the default location.
     
     * Unless you want to initialize Anaconda/Miniconda yourself select "yes" 
       when asked "Do you wish the installer to initialize Miniconda3 by 
       running conda init?"

#. **Install MELODIES MONET:** If you plan to modify the MELODIES MONET or MONET/MONETIO
   codebases, follow the :ref:`dev install instructions <dev-install-instructions>`.
   Otherwise, follow the :ref:`user install instructions <user-install-instructions>`.

#. **Link the cartopy shapefiles:** Hera has download restrictions,
   so link the required cartopy shapefiles 
   for plotting by running the ``link_cartopy_files.sh`` script.

   If you have cloned the repo (e.g. following the dev install instructions)::
       
      $ cd MELODIES-MONET/python_env_ymls/hera
      $ ./link_cartopy_files.sh

   If you didn't clone the repo and don't want to::

      $ wget -O - https://raw.githubusercontent.com/NOAA-CSL/MELODIES-MONET/main/python_env_ymls/hera/link_cartopy_files.sh | bash


**You are ready to start using and developing MELODIES MONET!**


.. note::
   In the recent past [#hera_no_pypi]_, Hera did not allow downloading
   from PyPI. As a result, such ``pip install``\s commands failed since pip was not
   able to download setuptools from PyPI.
   As a (reluctant) workaround, ``python setup.py develop`` can be used instead
   for editable (development) installs of MELODIES MONET and MONET/MONETIO.


.. note::
   In the recent past, downloading a lot of dependent packages at once
   with conda on Hera led to stalling.
   To overcome this challange, try installing packages individually::
  
        $ conda create --name melodies-monet python=3.9
        $ conda activate melodies-monet
        $ conda install -c conda-forge jupyterlab
        $ conda install -c conda-forge netcdf4
        $ conda install -c conda-forge wrf-python
        $ conda install -c conda-forge cartopy
        $ conda install -c conda-forge esmf
        $ conda install -c conda-forge monet
        $ conda install -c conda-forge monetio    


.. [#clone] See :ref:`the cloning notes <clone-notes>` if you have
   trouble cloning the repositories this way.

.. [#hera_no_pypi] Recent as of 12-Apr-2022. See :issue:`79`.
