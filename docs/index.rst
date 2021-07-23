Model and ObservatioN Evaluation Toolkit (MONET) Input Output (IO)
=================================================

**MONETIO** is an open source project and Python package that aims to create a
common platform for atmospheric composition data for weather and
air quality models.

MONET was developed to evaluate the Community Multiscale Air Quality Model (CMAQ)
for the NOAA National Air Quality Forecast Capability (NAQFC) modeling system.
From MONET version 2.1.4, MONETIO was broken off from MONET to be its own dedicated repository.
MONETIO is built to work in unison with MONET. For more information on MONET please refer to
https://monet-arl.readthedocs.io.

Our goals is to provide easy tools to retrieve and read atmospheric composition data in
order to speed scientific research.  Currently, MONETIO is able to process
several models and observations related to air composition and meteorology.

If you use MONETIO please reference the package.

Reference
^^^^^^^^^

Baker, Barry; Pan, Li. 2017. “Overview of the Model and Observation
Evaluation Toolkit (MONET) Version 1.0 for Evaluating Atmospheric
Transport Models.” Atmosphere 8, no. 11: 210

What's New
^^^^^^^^^^

MONETIO v0.1 has been released.  MONETIO has been split off from MONET to be able to
create a more focused repository for the I/O of models and observational data.

Features include:

  * fixes to observational datasets including Airnow, AQS, Aeronet, and ISH
  * New ICARTT reader developed to read into Xarray.Dataset
  * OpenAQ reader to read directly from the OpenAQ Amazon S3 server

.. toctree::
   :maxdepth: 4
   :caption: Getting Started

   why-monet
   installing
   observations
   models
   tutorial
   monetio_wcoss

Get in touch
------------

- Ask questions, suggest features or view source code `on GitHub`_.

.. _on GitHub: https://github.com/noaa-oar-arl/MONET


Supported datasets
------------------

**Supported Models**

* `HYSPLIT <https://www.ready.noaa.gov/HYSPLIT.php/>`_
* `CMAQ <https://www.epa.gov/cmaq/>`_
* `CAMx <http://www.camx.com/about/default.aspx/>`_
* FV3-CHEM (comming soon)
* WRF-CHEM (comming soon)

**Supported Observations**

* `AirNow <https://www.airnow.gov/>`_
* `AQS <https://www.epa.gov/aqs/>`_
* `AERONET <https://aeronet.gsfc.nasa.gov/>`_
* `CRN <https://www.ncdc.noaa.gov/crn/>`_
* `TOLNet <https://www-air.larc.nasa.gov/missions/TOLNet/>`_
* `CEMS <https://www.epa.gov/emc/emc-continuous-emission-monitoring-systems/>`_
* `IMPROVE <http://vista.cira.colostate.edu/Improve/>`_
* `ISH <https://www.ncdc.noaa.gov/isd/>`_





**Help & Reference**

.. toctree::
   :maxdepth: 10
   :caption: Help * Reference

   api
