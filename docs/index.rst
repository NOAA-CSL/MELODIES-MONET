MELODIES MONET
==============

**MELODIES MONET** is a joint project between NCAR and NOAA to develop a 
modular framework that integrates existing and future diverse atmospheric 
chemistry observational datasets with chemistry model results for the 
evaluation of air quality and atmospheric composition. MELODIES MONET combines 
the Model EvaLuation using Observations, DIagnostics and Experiments Software 
(MELODIES) project at NCAR with the Model and ObservatioN Evaluation Toolkit 
(MONET) project at NOAA to develop a python diagnostic package that is open 
source, generic, portable, and model-agnostic. Overall, the project provides a 
framework for evaluating a wide range of models in a more consistent manner. 
The tool is accessible to everyone including university and national 
laboratory researchers, as well as graduate students and postdocs.

The goal is to evaluate research, operational, and regulatory models against 
a variety of observations including surface, aircraft, and satellite data all 
within a common framework. MELODIES MONET uses the functionality already 
developed by MONETIO to read in multiple observational and model datasets and 
MONET to do pairing/analysis/plotting. For more information on MONET and 
MONETIO please refer to:
 
- https://monet-arl.readthedocs.io
- https://monetio.readthedocs.io


.. warning::
   MELODIES MONET is currently under development. The code is public to 
   encourage collaboration amongst the community. Do not publish results using 
   MELODIES MONET without consulting the development team.
   
.. note::
   Please cite the following to acknowledge use of MELODIES MONET

   - Baker, Barry; Pan, Li. 2017. “Overview of the Model and Observation Evaluation Toolkit (MONET) Version 1.0 for Evaluating Atmospheric Transport Models.” Atmosphere 8, no. 11: 210

   - MELODIES - TBD

.. toctree::
   :hidden:
   :maxdepth: 4
   :caption: User's Guide

   users_guide/introduction
   users_guide/description
   users_guide/supported_datasets
   users_guide/supported_diagnostics
   users_guide/supported_plots
   users_guide/supported_stats
   users_guide/time_chunking
   users_guide/gridded_datasets

.. toctree::
   :hidden:
   :maxdepth: 4
   :caption: Getting Started

   getting_started/installation
   getting_started/software_architecture
   getting_started/downloading_obs
   getting_started/how_to_run

.. toctree::
   :maxdepth: 1
   :caption: Contribute
   
   develop/contribute
   develop/other_resources
   develop/development_team
   develop/developers_guide
   develop/datasets

.. toctree::
   :hidden:
   :maxdepth: 4
   :caption: Examples

   examples/intro_examples
   examples/tutorial-data
   examples/airnow_wrfchem
   examples/airnow_wrfchem_reg
   examples/airnow_camchem_se
   examples/camchem
   examples/idealized
   examples/save_paired_data
   examples/read_paired_data
   examples/aircraft_pairing
   examples/AEROMMA_UFS-AQM_Plots
   
.. toctree::
   :hidden:
   :maxdepth: 4
   :caption: Current Applications

   applications/publications
   applications/forecasts

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Help and Reference

   api
   cli
   appendix/machine-specific-install
   appendix/yaml
   appendix/troubleshooting


