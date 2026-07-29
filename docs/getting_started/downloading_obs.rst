Downloading Observations
========================

As described below, some observations can be directly used in the MELODIES MONET tool as is 
and some need a preprocessing step to convert them into a consistent data format.

Surface
-------

Surface datasets commonly used for air quality and atmospheric composition applications are all in different 
formats and occasionally some HPC platforms including the NOAA Hera machine have download restrictions 
that prevent us from using the automatic download features available in MONET. So for now, 
MELODIES MONET has separate scripts to preprocess the surface observational datasets and save the output to an 
intermediate NetCDF file. These preprocess scripts are also useful so that users do not have to re-download 
observational data over and over again for the same analysis period. We will work on automating this process further 
in the future.

The MELODIES MONET tool has a Command Line Interface (CLI) that can be used to download and create 
MELODIES MONET-ready datasets for: AirNow, AERONET, AQS, ISH, ISH-Lite, and OpenAQ. For users experiencing extended 
download times, refer to the note below. 

The Command Line Interface allows users to very easily download datasets with a single command line argument. 
Generally, users only need to select which subcommand to use (i.e., which observational data set you want to download) 
and then specify the start date and end date like that below to download all US EPA AQS observations in August 2023::

    $ melodies-monet get-aqs -s 2023-08-01 -e 2023-09-01

The other datasets can be downloaded in the same way::

    $ melodies-monet get-aeronet -s 2023-08-01 -e 2023-09-01
    $ melodies-monet get-airnow -s 2023-08-01 -e 2023-09-01
    $ melodies-monet get-ish -s 2023-08-01 -e 2023-09-01
    $ melodies-monet get-ish-lite -s 2023-08-01 -e 2023-09-01
    $ melodies-monet get-openaq -s 2023-08-01 -e 2023-09-01

The Command Line Interface will default to compressing the dataset, which can significantly save space. However, this
compression step also takes time and some users have run into problems. Users can easily turn this compression off 
by adding ``--no-compress``::

    $ melodies-monet get-aqs -s 2023-08-01 -e 2023-09-01 --no-compress

There are many other optional features available that are fully described in the Appendix :doc:`/cli`.

.. note::
   For users using MELODIES MONET on the NOAA Hera machine (or other machines 
   with download restrictions), you will need to use the MELODIES MONET Command Line Interface on a 
   machine without download restrictions and manually copy the netCDF files produced 
   onto the NOAA Hera machine.

.. note::
   For users using ISH and ISH-LITE, download times can be significant. Running this command in the CLI is be reccomended 
   if you are experiencing extended wait times while downloading any dataset. Connection times between MELODIE-MONET and the 
   data server can time out. With this loop, the download can fail, but will automatically try again::

   $ until melodies-monet get-ish-lite -s 2017-07-01 -e 2017-07-03 --no-compress; do     echo "Download failed, retrying in 10 seconds...";     sleep 10; done

   More generally::

   $ until melodies-monet get-<dataset> -s <start_date> -e <end_date> --no-compress; do     echo "Download failed, retrying in 10 seconds...";     sleep 10; done

Aircraft, Sonde, Mobile, and Ground Campaign Data
-------------------------------------------------

Aircraft, sonde, mobile, and ground campaign data can be used directly in the tool as long 
as the data format is NetCDF, `ICARTT <https://www-air.larc.nasa.gov/missions/etc/IcarttDataFormat.htm>`_, or CSV. Users download their own observational data. 
No pre-processing is required for these datasets.

Satellite
---------

For satellite data, users download their own observational data. No pre-processing is required 
for these datasets.