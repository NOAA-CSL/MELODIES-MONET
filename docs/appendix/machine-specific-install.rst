Machine-specific Install
========================

NCAR HPC Cheyenne
-----------------

Below is a specific recipe for how to set up all the necessary Python dependencies 
on the NCAR HPC, Cheyenne/Casper. Note: these are developer-specific instructions. 
There are three steps to complete before you can use and develop MELODIES-MONET 
on Cheyenne/Casper: **Step 1.** Installing the conda package manager 
Miniconda, **Step 2.** Creating a conda environment with all required dependencies, 
and **Step 3.** Cloning MELODIES MONET code.

We will use the conda package manager system to create a contained Python environment 
for running and developing MELODIES-MONET. This creates a “stand-alone” instance 
of interdependent packages that will not interfere with your access to the main 
installation of Python on the system. These instructions should work for a bash 
or tcsh shell.

**Step 0:** Log into the NCAR HPC

**Step 1 Miniconda:** The conda package manager is not installed on the Cheyenne 
or Casper clusters, but you can install Miniconda in your own user space 
(`CISL instructions <https://www2.cisl.ucar.edu/resources/conda-environments>`_):

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

Accept the license and default install paths. Say yes to running conda init. 
Log out of Cheyenne/Casper and log back in. You may see '(base)' as the start 
of your terminal prompt. If you have any issues with package conflicts you may 
need to uninstall miniconda and reinstall.

**Step 2 Dependent Python Packages:** Set up a conda environment with required 
dependencies.

(a) Set up and activate CONDA environment specific to MELODIES-MONET. You can 
    call this environment whatever you like, we suggest ‘melodies-monet’::

    $ conda create --name melodies-monet python=3.9
    $ conda activate melodies-monet

    You should see ‘(melodies-monet)’ at the start of your terminal prompt.

(b) Download the following packages step-by-step. Note they have sub-packages 
    that will be downloaded. The ‘-y’ means you will not have to interactively 
    choose ‘y’ to proceed with downloading packages::

    $ conda install -y -c conda-forge netcdf4
    $ conda install -y -c conda-forge wrf-python
    $ conda install -y -c conda-forge notebook
    $ conda install -y -c conda-forge jupyterlab
    $ conda install -y -c conda-forge monet
    $ conda install -y -c conda-forge monetio

    Some main packages that are downloaded with the monet install: *cartopy, 
    cython, dask, markupsafe, matplotlib-base, pandas, pydecorate, pyresample, 
    python-stratify, seaborn, xarray, xesmf*. Some main packages that are 
    downloaded with the jupyter install: *icu, ipykernel, ipython*. If other 
    systems have timeout errors, downloading these separately may improve the 
    setup performance.

(c) Once you set up the correct packages through conda, and while your conda 
    environment is still activated, get the most recent branches of MONET and 
    MONETIO using github, and link them with conda. This is done because MONET 
    and MONETIO are still in active development. Create a ‘monet-base’ folder 
    (e.g. in your work location on cheyenne).

    Set up and link MONET within monet-base::

    $ git clone https://github.com/noaa-oar-arl/monet.git
    $ cd monet
    $ git checkout develop
    $ pip install -e .

    Set up and link MONET IO within monet-base::

    $ git clone https://github.com/noaa-oar-arl/monetio.git
    $ cd monetio
    $ git checkout development
    $ pip install -e .

**Step 3: Clone the MELODIES-MONET Github branch** ::

    $ git clone git@github.com:NOAA-CSL/MELODIES-MONET.git

    End step. At the end of working with MELODIES-MONET, deactivate the 
    conda environment::

    $ conda deactivate

**You are ready to start developing MELODIES-MONET!**

NOAA HPC Hera
-------------

Below is a specific recipe for how to set up all the necessary Python 
dependencies on the NOAA Hera machine. There are three steps to complete 
before you can use and develop MELODIES MONET on hera: **Step 1.** Install 
the conda package manager Anaconda/Miniconda, **Step 2.** Clone MELODIES MONET 
code, and **Step 3.** Create a conda environment with all required dependencies

We will use the conda package manager system to create a contained Python 
environment for running and developing MELODIES MONET. 

#. **Download Anaconda/Miniconda:** Follow the instructions
   `on the RDHPCS wiki <https://rdhpcs-common-docs.rdhpcs.noaa.gov/wiki/index.php/Anaconda>`__
   to download and run anaconda/miniconda on Hera. Tips for success:

   * You will need a NOAA HPC account to access the RDHPCS wiki link above.

   * Both anaconda/miniconda will work well for MELODIES MONET. See
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
     
     * Unless you want to initialize anaconda/miniconda yourself select "yes" 
       when asked "Do you wish the installer to initialize Miniconda3 by 
       running conda init?"

   * Follow the `github ssh key instructions <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`__
     to add an ssh key on Hera.

#. **Clone the MELODIES MONET Github package**::

    $ git clone git@github.com:NOAA-CSL/MELODIES-MONET.git

#. **Create a conda environment with the required dependencies on Hera:** 

   * Follow either option 1 below by using an example ``environment.yml`` file from 
     the MELODIES MONET package or follow option 2 below to set this up manually.
     
     **Option 1: Use an example environment.yaml file:**

       * Make a copy of the environment.yaml file for Hera stored in the 
         MELODIES MONET Github repository 
         (MELODIES_MONET/python_env_ymls/hera/environment.yml). If needed, 
         update the first line to change the default environment name. Also 
         update the last line to point to your own anaconda/miniconda directory 
         location and if needed update the default environment name.

       * Run the following, to create the environment. ::
    
          $ conda env create -f environment.yml

       * Verify the new environment exists ::
    
          $ conda env list 

       * Activate the new environment :: 
    
          $ conda activate py36_monet_def
     
     **Option 2: Manual method:** 
       
      * Downloading a lot of dependent packages at once on Hera leads to stalling. 
        To overcome this challange, either use Option 1 or install some of the 
        larger packages first and then install MONET and MONETIO like the following: ::
   
        $ conda create --name py39_monet_def python=3.9
        $ conda activate py39_monet_def
        $ conda install -c conda-forge notebook
        $ conda install -c conda-forge jupyterlab
        $ conda install -c conda-forge netcdf4
        $ conda install -c conda-forge wrf-python
        $ conda install -c conda-forge cartopy
        $ conda install -c conda-forge esmf
        $ conda install -c conda-forge monet
        $ conda install -c conda-forge monetio        
        
   * Note: There are instances where other packages will be needed. These are 
     just to download the basics, so if you get an error about missing a 
     package install it in your conda environment.
    
   * Once you have a working and activated conda environment, you will need to 
     link the latest versions of MONET and MONETIO from Github. ::
   
      $ git clone git@github.com:noaa-oar-arl/monet.git
      $ cd monet
      $ git checkout develop
      $ pip install -e .
    
      $ git clone git@github.com:noaa-oar-arl/monetio.git
      $ cd monetio
      $ git checkout development
      $ pip install -e .

   * Hera has download restrictions, so link the required cartopy shapefiles 
     for plotting by running the following script ::
       
      $ cd MELODIES-MONET/python_env_ymls/hera
      $ ./link_cartopy_files.sh

**You are ready to start using and developing MELODIES MONET!**
