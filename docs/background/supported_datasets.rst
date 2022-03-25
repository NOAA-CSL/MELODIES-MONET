Supported Datasets
==================

Supported Models and Observations are below. Please see
:ref:`Adding New Datasets <develop/datasets:Adding New Datasets>`
for advice on how to add new model and observational datasets to MELODIES MONET.

Supported Models
----------------

   * `CMAQ <https://www.epa.gov/cmaq/>`_
   * `WRF-Chem <https://ruc.noaa.gov/wrf/wrf-chem/>`_
   * `CESM/CAM-chem <https://www2.acom.ucar.edu/gcm/cam-chem>`_
   * RRFS-CMAQ (under development)
   * GSD-Chem (under development)

In general, processing requires input to be in netCDF format. For the above 
models, scripts to configure the model data into a standard format for 
MELODIES MONET are available. If input datasets are in netCDF format and  
define latitude, longitude, altitude, and a datetime object, MELODIES MONET may be able 
to directly read the input files.

See the `Expand models <https://github.com/NOAA-CSL/MELODIES-MONET/projects/8>`_ 
project on GitHub to learn about current and future development.

Supported Observations
----------------------

Surface
^^^^^^^

   * `AirNow <https://www.airnow.gov/>`_ 
   * `AERONET <https://aeronet.gsfc.nasa.gov/>`_ (under development)
   * `IMPROVE <http://vista.cira.colostate.edu/Improve/>`_ (under development)
   * `AQS <https://www.epa.gov/aqs/>`_ (in MONET, coming soon to MELODIES MONET)
   * `CRN <https://www.ncdc.noaa.gov/crn/>`_ (in MONET, coming soon to MELODIES MONET)
   * `TOLNet <https://www-air.larc.nasa.gov/missions/TOLNet/>`_ 
     (in MONET, coming soon to MELODIES MONET)
   * `CEMS <https://www.epa.gov/emc/emc-continuous-emission-monitoring-systems/>`_ 
     (in MONET, coming soon to MELODIES MONET)
   * `ISD <https://www.ncdc.noaa.gov/isd/>`_ (in MONET, coming soon to MELODIES MONET)
   
See the `Expand Surface Observations <https://github.com/NOAA-CSL/MELODIES-MONET/projects/7>`_ 
project on GitHub to learn about current and future development.

Aircraft (under development)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * `FIREX-AQ <https://csl.noaa.gov/projects/firex-aq/>`_ (under development)
   * `ATom <https://espo.nasa.gov/atom/content/ATom>`_ (under development)
   
See the `Incorporate Aircraft Evaluation <https://github.com/NOAA-CSL/MELODIES-MONET/projects/4>`_ 
project on GitHub to learn about current and future development.

Satellite (under development)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See the `Incorporate Satellite Evaluation <https://github.com/NOAA-CSL/MELODIES-MONET/projects/6>`_ 
project on GitHub to learn about current and future development.
