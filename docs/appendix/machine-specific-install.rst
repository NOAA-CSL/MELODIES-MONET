Machine-specific Install
========================

NCAR HPC Cheyenne/Casper
------------------------

Below are specific recipes for getting started with MELODIES MONET
on the NCAR HPC, Cheyenne/Casper.

**Official NCAR JupyterHub kernel**

One option is to use the "melodies-monet" kernel [#ncar_jhub_kernel]_ on the NCAR JupyterHub.

#. Navigate to https://jupyterhub.hpc.ucar.edu in your web browser.
#. Select "Production" and log in with your NCAR CIT credentials.
#. Start a server.
#. In the Notebook section of the Launcher, select "melodies-monet"
   to create a new notebook using the kernel.

Since you don't have to install anything yourself, this is the easiest way to get started.

**Personal conda environment**

Another option is to use conda to create your own MELODIES MONET installation
in a conda environment.
This creates a "stand-alone" instance 
of interdependent packages that will not interfere with your access to the main 
installation of Python on the system.
You can use the
`NCAR-maintained conda installation <https://arc.ucar.edu/knowledge_base/83853599>`__
to get access to ``conda`` by invoking::

    $ module load conda/latest

Or,
`install your own copy <https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html>`__
of Miniconda or Anaconda.

Once you have conda, you are ready to install MELODIES MONET.
If you plan to modify the MELODIES MONET or MONET/MONETIO
codebases, follow the :ref:`dev install instructions <dev-install-instructions>`.
Otherwise, follow the :ref:`user install instructions <user-install-instructions>`.

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


.. [#ncar_jhub_kernel] Maintained by NCAR CISL and members of the MELODIES MONET development team.

.. [#clone] See :ref:`the cloning notes <clone-notes>` if you have
   trouble cloning the repositories this way.

.. [#hera_no_pypi] Recent as of 12-Apr-2022. See :issue:`79`.
