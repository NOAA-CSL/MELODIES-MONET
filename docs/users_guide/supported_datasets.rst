Supported Datasets
==================

Supported Models and Observations are below. Please see
:ref:`Adding New Datasets <develop/datasets:Adding New Datasets>`
for advice on how to add new model and observational datasets to MELODIES MONET.

Supported Models
----------------

.. list-table:: Currently Connected Capabilities for Model Readers
   :widths: 20 20 30 30
   :header-rows: 1

   * - Model
     - Surface
     - Aircraft
     - Satellite
   * - `MERRA2 <https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/>`_
     - Yes
     - Needs testing
     - MODIS
   * - `WRF-Chem <https://www2.acom.ucar.edu/wrf-chem>`_
     - Yes
     - Yes
     - TROPOMI, TEMPO
   * - `CESM/CAM-chem FV <https://www2.acom.ucar.edu/gcm/cam-chem>`_
     - Yes
     - Needs testing
     - Needs testing
   * - `CESM/CAM-chem SE <https://www2.acom.ucar.edu/gcm/cam-chem>`_
     - Yes
     - | Needs testing & to 
       | add unstructured 
       | grid capabilities
     - | Needs testing & to 
       | add unstructured 
       | grid capabilities
   * - `CMAQ <https://www.epa.gov/cmaq/>`_
     - Yes
     - Needs testing
     - Needs testing
   * - `UFS-AQM (RRFS) <https://github.com/ufs-community/ufs-srweather-app/wiki/Air-Quality-Modeling>`_
     - Yes
     - Yes
     - Needs testing
   * - `CAMx <https://www.camx.com/>`_
     - Yes
     - Needs testing
     - TROPOMI, TEMPO
   * - `RAQMS <http://raqms-ops.ssec.wisc.edu/>`_
     - Yes
     - Needs testing
     - MOPITT, OMPS

In general, processing requires input to be in netCDF format. For the above 
models, scripts to configure the model data into a standard format for 
MELODIES MONET are available. If input datasets are in netCDF format and  
define latitude, longitude, altitude, and a datetime object, MELODIES MONET may be able 
to directly read the input files.

See the `Expand models in MELODIES-MONET <https://github.com/orgs/NOAA-CSL/projects/6>`_ 
project on GitHub to learn about current and future development.

Supported Observations
----------------------

Surface
^^^^^^^

   * `AirNow <https://www.airnow.gov/>`_ 
   * `AERONET <https://aeronet.gsfc.nasa.gov/>`_
   * `IMPROVE <http://vista.cira.colostate.edu/Improve/>`_ (under development)
   * `AQS <https://www.epa.gov/aqs/>`_ (in MONET, coming soon to MELODIES MONET)
   * `CRN <https://www.ncdc.noaa.gov/crn/>`_ (in MONET, coming soon to MELODIES MONET)
   * `TOLNet <https://www-air.larc.nasa.gov/missions/TOLNet/>`_ 
     (in MONET, coming soon to MELODIES MONET)
   * `CEMS <https://www.epa.gov/emc/emc-continuous-emission-monitoring-systems/>`_ 
     (in MONET, coming soon to MELODIES MONET)
   * `ISD <https://www.ncei.noaa.gov/products/land-based-station/integrated-surface-database>`_
     (in MONET, coming soon to MELODIES MONET)

See the `Expand Surface Observations in MELODIES-MONET <https://github.com/orgs/NOAA-CSL/projects/6>`_ 
project on GitHub to learn about current and future development.

.. note::

   The :doc:`/cli` can be used to download and create MELODIES MONET-ready datasets for:
   AirNow, AERONET, AQS, and ISD.

Aircraft (under development)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * `FIREX-AQ <https://csl.noaa.gov/projects/firex-aq/>`_ (under development)
   * `ATom <https://espo.nasa.gov/atom/content/ATom>`_ (under development)
   
See the `Incorporate Aircraft Evaluation in MELODIES-MONET <https://github.com/orgs/NOAA-CSL/projects/6>`_ 
project on GitHub to learn about current and future development.

Satellite (under development)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See the `Incorporate Satellite Evaluation in MELODIES-MONET <https://github.com/orgs/NOAA-CSL/projects/6>`_ 
project on GitHub to learn about current and future development.
