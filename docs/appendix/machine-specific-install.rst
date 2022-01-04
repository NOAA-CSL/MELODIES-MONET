Machine-specific Install
=========================

NCAR HPC cheyenne
--------------------

Below is a specific recipe for how to set up all the necessary Python dependencies on the NCAR HPC, cheyenne/casper. Note: these are developer-specific instructions. There are three steps to complete before you can use and develop MELODIES-MONET on cheyenne/casper: **Step 1.** Installing the conda package manager Miniconda, **Step 2.** Creating a conda environment with all required dependencies, and **Step 3.** Cloning MELODIES-MONET code.

We will use the conda package manager system to create a contained Python environment for running and developing MELODIES-MONET. This creates a “stand-alone” instance of interdependent packages that will not interfere with your access to the main installation of Python on the system. These instructions should work for a bash or tcsh shell.

**Step 0:** Log into the NCAR HPC

**Step 1 Miniconda:** The conda package manager is not installed on the Cheyenne or Casper clusters, but you can install Miniconda in your own user space (`CISL instructions <https://www2.cisl.ucar.edu/resources/conda-environments>`_):

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

Accept the license and default install paths. Say yes to running conda init. Log out of cheyenne/casper and log back in. You may see '(base)' as the start of your terminal prompt. If you have any issues with package conflicts you may need to uninstall miniconda and reinstall.

**Step 2 Dependent Python Packages:** Set up a conda environment with required dependencies.

(a) Set up and activate CONDA environment specific to MELODIES-MONET. You can call this environment whatever you like, we suggest ‘monet_py36’::

    $ conda create --name monet_py36 python=3.6
    $ conda activate monet_py36

    You should see ‘(monet_py36)’ at the start of your terminal prompt.

(b) Download the following packages step-by-step. Note they have sub-packages that will be downloaded. The ‘-y’ means you will not have to interactively choose ‘y’ to proceed with downloading packages::

    $ conda install netcdf4
    $ conda install -y -c conda-forge monet
    $ conda install -y -c conda-forge monetio
    $ conda install -y -c conda-forge wrf-python
    $ conda install -y -c conda-forge jupyter

    Some main packages that are downloaded with the monet install: *cartopy, cython, dask, markupsafe, matplotlib-base, pandas, pydecorate, pyresample, python-stratify, seaborn, xarray, xesmf*. Some main packages that are downloaded with the jupyter install: *icu, ipykernel, ipython*. If other systems have timeout errors, downloading these separately may improve the setup performance.

(c) Once you set up the correct packages through conda, and while your conda environment is still activated, get the most recent branches of MONET and MONETIO using github, and link them with conda. This is done because MONET and MONETIO are still in active development. Create a ‘monet-base’ folder (e.g. in your work location on cheyenne).

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

**Step 3: Clone the MELODIES-MONET Github branch**

(a) Fork the Github repository to your own Github account:

    https://github.com/NOAA-CSL/MELODIES-MONET

    Note: you can pull updates from the main NOAA repository by using the “Fetch upstream” button on your fork. Alternatively::

    $ git pull upstream master
    $ git push origin master

(b) Navigate on cheyenne/casper to where you would like to keep the MELODIES-MONET code (e.g. in your work location) and clone your fork to cheyenne::

    $ git clone https://github.com/$GithubUsername/$ForkName.git

(c) Checkout the develop_plots branch - you need to do this with the remote branch as well as create a local tracking branch::

    $ git checkout origin/develop_plots
    $ git checkout develop_plots

    Then all development work will be in the monet_analysis folder::

    $ cd monet_analysis

    End step. At the end of working with MELODIES-MONET, deactivate the conda environment::

    $ conda deactivate


**You are ready to start developing MELODIES-MONET!**

NOAA HPC hera
-------------

Below is a specific recipe for how to set up all the necessary Python 
dependencies on the NOAA hera machine. Note: these are developer-specific 
instructions. There are three steps to complete before you can use and develop 
MELODIES-MONET on hera: **Step 1.** Install the conda package manager 
Anaconda/Miniconda, **Step 2.** Clone MELODIES-MONET code, and **Step 3.** 
Create a conda environment with all required dependencies

We will use the conda package manager system to create a contained Python 
environment for running and developing MELODIES-MONET. 

#. **Download Anaconda/Miniconda:** Follow the instructions `on the RDHPCS wiki`_ 
   to download and run anaconda/miniconda on hera. Tips for sucess:

   * You will need a NOAA HPC account to access the RDHPCS wiki link above.

   * Select **YES** at this step in the wiki: “4. After the installation is 
     complete, the installer will ask to initialize in your .bashrc/cshrc - 
     Select yes or no.” 

   * Both anaconda/miniconda will work well for MELODIES-MONET. See 
     `conda instructions'_ to determine, which is the best option for you.

   * Installing anaconda/miniconda on scratch is recommended due to the limited 
     space available on your home directory.

   * Follow the `github ssh key instructions`_ to add an ssh key on hera.

#. **Clone the MELODIES-MONET Github package**::

    $ git clone git@github.com:NOAA-CSL/MELODIES-MONET.git


#. **Create a conda environment with the required dependencies on hera:** 
   Due to download restrictions, hera cannot download a lot of dependent 
   python packages at once, so it is highly recommended to setup your 
   initial conda environment from a working environment.yml file as outline below:

   * Make a copy of the environment.yaml file for hera stored in the 
     MELODIES-MONET github repository 
     (MELODIES_MONET/python_env_ymls/hera/environment.yml). If needed, update 
     the first line to change the default environment name. Also update the 
     last line to point to your own anaconda/miniconda directory location and 
     if needed update the default environment name.

   * Run the following, to create the environment. Note this takes 15-30 
     minutes, so be patient. ::
    
    $ conda env create -f environment.yml

   * Verify the new environment exists::
    
    $ conda env list 

   * Activate the new environment:: 
    
    $ conda activate py36_monet_def

   * To use the latest versions of MONET and MONETIO from github. Clone and 
     link them to your conda environment::
   
    $ git clone git@github.com:noaa-oar-arl/monet.git
    $ cd monet
    $ git checkout develop
    $ pip install -e .
    
    $ git clone git@github.com:noaa-oar-arl/monetio.git
    $ cd monetio
    $ git checkout development
    $ pip install -e .

  * Link the required cartopy shapefiles for plotting::
    $ cd MELODIES-MONET/python_env_ymls/hera
    $ ./link_cartopy_files.sh

**You are ready to start developing MELODIES-MONET!**

.. _on the RDHPCS wiki: https://rdhpcs-common-docs.rdhpcs.noaa.gov/wiki/index.php/Anaconda#Installation
.. _conda instructions: https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html#anaconda-or-miniconda
.. _github ssh key instructions: https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account
