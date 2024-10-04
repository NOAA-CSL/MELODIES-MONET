Adding New Datasets
===================

Observations
------------

Examples for how to read in data from Aeronet, AirNow, and Improve are in the
``examples/process_obs`` folder in the MELODIES MONET repository
`on GitHub <https://github.com/NOAA-CSL/MELODIES-MONET>`__.
Use these examples as reference in order to add new surface observational datasets.

Instructions for reading in aircraft and satellite observations are under development. 

Models
------
Examples for reading model datasets can be
found in the ``monetio/models`` folder in the MONETIO repository
`on GitHub <https://github.com/noaa-oar-arl/monetio>`__.
These include e.g., _cesm_fv_mm.py, _cmaq_mm.py, and _wrfchem_mm.py.
While a part of the MONETIO repository,
the private MELODIES MONET readers are designated with prefix ``_mm``.

Support for additional models is also under developed.

Standard variables are required to be computed in each model reader for each capability including surface, aircraft, and satellite as specified in the table below.

.. list-table:: Required Variables for Model Readers
   :widths: 10 30 30 30
   :header-rows: 1

   * - Capability
     - | Variable Name 
       | in Code
     - Description
     - Additional Requirements
   * - Surface
     - | ``time``
       | ``latitude``
       | ``longitude``
     - | Time in ``datetime64[ns]`` format
       | Latitude in degrees
       | Longitude in degrees
     - | Provide only surface model data 
       | or if provide vertical model data, 
       | first level must be the level 
       | nearest to the surface.
       | All gases are in ppb and 
       | all aerosols are in µg/m3.
   * - Aircraft
     - | ``time``
       | ``latitude``
       | ``longitude``
       | ``pres_pa_mid``
       | ``temperature_k``
     - | Time in ``datetime64[ns]`` format
       | Latitude in degrees
       | Longitude in degrees
       | Mid-level pressure in pascals (Pa)
       | Mid-level temperature in kelvin (K)
     - | Provide vertical model data. 
       | All gases are in ppb and 
       | all aerosols are in µg/m3.
   * - Satellites
     - | ``time``
       | ``latitude``
       | ``longitude``
       | ``pres_pa_mid``
       | ``temperature_k``
       | ``dz_m``
       | ``surfpres_pa``
     - | Time in ``datetime64[ns]`` format
       | Latitude in degrees
       | Longitude in degrees
       | Mid-level pressure in pascals (Pa)
       | Mid-level temperature in kelvin (K)
       | Layer thickness in meters (m)
       | Surface pressure in pascals (Pa)
     - | Provide vertical model data.