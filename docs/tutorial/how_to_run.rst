How To Run
==========

These are instructions for how to update the examples on GitHub and run 
MELODIES MONET both within a jupyter notebook or in a bash script submitted to 
an HPC machine. It is highly recommended to start first running MELODIES 
MONET in a jupyter notebook with a smaller dataset. Then you can expand on the 
analysis to produce hundreds of plots by submitting a bash script to your HPC 
machine. The basic commands in the jupyter notebook and bash script are exactly 
the same and both similarly call an input YAML file.

First, we describe how to `Prepare An Input YAML File`_.
Second, we define the `Basic Commands`_ to run MELODIES MONET.
Then, we describe how to start from either a 
the `Jupyter Notebook`_ or
the `Bash Script`_ example.

Prepare an Input YAML File
--------------------------
You will need to prepare an input YAML file to be read into MELODIES MONET. 
Example input YAML files to start from are provided in the 
``examples/yaml`` folder of the code on GitHub. There are a number 
of comments in these example input YAML files to get you started. The overall 
structure is the following:

* **analysis** -- All input related to the analysis class.

* **model** -- All input for each instance of the model class. The variables
  to be plotted are included in the model class in the "mapping" dictionary. 
  The model variable names are first (i.e. keys) and the observation variable 
  names are the second (i.e. values). Because the plots in MELODIES MONET 
  will plot multiple models with one observation, the observation variables 
  listed in the mapping dictionary must be consistent across all models. 
  For example, if you want to plot the results of multiple model datasets 
  against the AirNow observations for "OZONE" and "PM2.5", you must 
  provide the model variable names for "OZONE" and "PM2.5" in the mapping 
  dictionary for all models. Say if you only provide the model variable 
  names for "OZONE" for one of the models, MELODIES MONET will error. 
* **obs** -- All input for each instance of the observation class.

* **plots** -- All input for each plotting group. A plotting group consists 
  of one plotting type. The plotting types are described in 
  :doc:`/background/supported_plots`. All model /
  observational pairs and domains specified for the plotting group will be 
  included. You may include as many plotting groups as you like.

* **stats** -- All input needed to calculate the statistics. The supported
  statistics available in MELODIES MONET are described in 
  :doc:`/background/supported_stats`. All model /
  observational pairs and domains specified will be included. You may 
  include as many statistics as you like. Note however that the calculation 
  of the statistics is relatively slow right now. Optimizing this code is 
  under development.

A detailed description of all the options in the input YAML file is provided 
in the Appendix under :doc:`/appendix/yaml`.

Basic Commands
--------------

First, you will import the MELODIES MONET driver. You will need to update the path
below to point to the location of your MELODIES-MONET GitHub repository. ::

    import sys; sys.path.append("../../")
    from melodies_monet import driver

Then you will create an instance of the python driver analysis class. The 
analysis class consists of 4 main parts; information read in by the input YAML 
file, instances of the model class, instances of the observation class, 
instances of the paired class. This instance of the analysis class will allow 
us to store all the information we need to create the plots and calculate the 
statistics. ::

    an = driver.analysis()

Then you will provide the path and name to your input YAML file and read the 
information from the input YAML file into the analysis class. ::

    an.control = 'control_cmaq.yaml'
    an.read_control()

Then you will read in all models listed in the "model" section of the input 
YAML file. This will create an instance of the model class for each model, 
which includes information like the label, type, mapping table, filenames, 
and file paths. Note: multiple files can be opened at the same time using hot 
keys ::
    
     an.open_models()

Then you will read in all observations listed in the "obs" section of the input 
YAML file. This will create an instance of the observation class for each 
observation dataset, which includes information like the label, type, filenames, 
and file paths. MONET can automatically download observational data. For now, 
MELODIES MONET will open preprocessed data because some HPC platforms have 
download restrictions. We will work on automating this process further in the 
future. See :doc:`downloading_obs` to learn how to
preprocess the observational data for MELODIES MONET ::

     an.open_obs()
     
Pair all of the models with all of the observations listed in the input YAML 
file creating instances of the paired class. Model results are directly paired 
in both time and space to the observational dataset for all variables listed 
in the mapping dictionary specified in the model class. ::

     an.pair_data()
     
For all plotting groups, specified in the input YAML file, create all plots 
looping over all specified variables in the mapping table, all model/observational
pairs, and all domains. ::

     an.plotting()

Calculate all statistics specified in the input YAML file looping over all 
specified variables in the mapping table, all model/observational pairs, and 
all domains. ::

     an.stats()
     
Jupyter Notebook
----------------

Jupyter notebook examples explaining how to run MELODIES MONET are in the 
``examples/jupyter_notebooks`` folder of the code on GitHub. It is
highly recommended for new and expert users to first use MELODIES MONET on a 
subset of the analysis in jupyter notebook. The jupyter notebook examples 
demonstrate how to print different pieces of the analysis class instance in 
order to help trouble shoot problems. In order to print figures in the jupyter 
notebook in the analysis section of the input YAML file, set debug = True.

Bash Script
-----------

Jupyter notebooks are great for quick analysis and ensuring you have set up the 
configuration properly, but if you want to perform the analysis for a longer 
time period or create hundreds of plots submitting a bash script as a job on 
your HPC computer is preferred. Bash script examples for running MELODIES MONET 
are in the ``examples/submit_jobs`` folder of the code on GitHub.

   * If you are using a model like WRF-Chem, CMAQ, or RRFS-CMAQ that is run in 
     forecasting mode and you want to combine model results across multiple 
     days or even over an entire month, you may need to link model data into 
     a directory first to ensure that you have sequential model results to 
     incorporate into MELODIES MONET. Examples of bash scripts for doing this 
     are provided (link_files_*.sh). 

   * Then you will need to copy and update the run_melodies_monet.py script.
     
        - Update to include the path and file name for your input YAML file. 
        - This script defaults to running both the plotting and stats routines. 
          If you only want to perform one or the other, comment one of them out. 

   * Then copy and edit the submit_hera.sh script. This is an example of how to 
     submit the job to the NOAA Hera machine. Edit this script to be appropriate 
     for your HPC machine. Note: you may need to use a larger memory node to run 
     MELODIES MONET.
     
        - Update the location of your conda environment. 
        - Also update the location and name of your run_melodies_monet.py script.

   * Submit the submit_hera.sh script (e.g., on NOAA hera: sbatch submit_hera.sh)	 

