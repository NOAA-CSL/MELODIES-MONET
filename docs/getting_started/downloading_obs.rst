Downloading Observations
========================

MONET can automatically download observational data. Some HPC platforms including 
the NOAA Hera machine have download restrictions that prevent us from using the 
automatic download feature. So for now, MELODIES MONET has separate scripts to 
preprocess the observational datasets. These preprocess scripts are also useful
so that users do not have to re-download observational data over and over again 
for the same analysis period. We will work on automating this process further 
in the future. 

Examples of preprocessed observational data for MELODIES MONET are here:
`MELODIES MONET Example Datasets <https://csl.noaa.gov/groups/csl4/modeldata/melodies-monet/>`_. 

In order to preprocess the observational data for additional time periods 
follow the instructions in the jupyter notebooks in the 
``examples/process_obs`` folder of the code on GitHub. Examples for 
the following observational datasets are provided.

   * Aeronet
   * AirNow
   * Improve

Adapting these scripts for other observational datasets should be straight 
forward.

Note: for users using MELODIES MONET on the NOAA Hera machine (or other machines 
with download restrictions), you will need to run these jupyter notebooks on a 
machine without download restrictions and manually copy the netCDF files produced 
onto the NOAA Hera machine.




