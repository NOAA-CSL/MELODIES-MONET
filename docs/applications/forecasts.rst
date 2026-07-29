Evaluating Air Quality Forecasts
================================

RAP-chem
--------

Jordan Schnell at NOAA GSL is using MELODIES MONET to evaluate the performance 
of a variety of research and operational air quality forecast models 
including a newly developed research forecast model called RAP-chem every day 
in near real time against the AirNow surface observations of ozone, PM\ :sub:`2.5`\, 
CO, and NO\ :sub:`2` and AERONET AOD measurements.

Check out the full analysis for each forecast `here <https://rapidrefresh.noaa.gov/monet_rrfs_verif/>`__.

This includes a new feature developed by NOAA GSL summer student, Mackenzie Arnold,
to `interactively view plots for individual surface sites online <https://rapidrefresh.noaa.gov/monet_rrfs_verif/>`__.

The code to produce this analysis using MELODIES MONET is in the
``examples/forecast_evaluation`` folder on GitHub.

Example plots for ozone and PM\ :sub:`2.5` for the forecast on January 18th, 2022 
are below.

.. figure:: /_static/figures/OZONE_EPA_f00_rapchem.png

.. figure:: /_static/figures/PM25_EPA_f00_rapchem.png

UFS-AQM (AQM-Eval)
------------------

The UFS-AQM scientific team, in collaboration with NOAA-EPIC, developed a standardized MELODIES MONET model evaluation wrapper optimized for rocoto-based workflows running in HPC environments. The suite comprises packages targeting key measures of atmospheric composition model performance: chemistry, meteorology, particulate matter, and volatile organic compounds (VOCs). Automated workflow generation is currently implemented for the `UFS-Short-range Weather App (UFS-SRW) <https://ufs-srweather-app.readthedocs.io/en/develop/UsersGuide/BuildingRunningTesting/AQM.html#melodies-monet-mm-evaluation>`__ with plans to generalize workflow creation and configuration to connect to a wider variety of modeling systems.

Code and wiki-based documentation are hosted on GitHub in NOAA-EPIC’s `AQM-Eval <https://github.com/NOAA-EPIC/AQM-Eval/wiki/aqm%E2%80%90mm%E2%80%90eval>`__ repository.

Below is a sequence diagram providing an overview of major operations in the MM wrapper.

.. figure:: /_static/figures/aqm-mm-eval-sequence.png