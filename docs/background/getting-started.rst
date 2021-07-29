Getting Started
===============

The user interface
------------------
YAML (YAML Ain't Markup Language) is the format of the MELODIES configuration file. The yaml configuration file is how the general user interacts with MELODIES-MONET. A user edits the yaml file and then uses the driver python script to run the yaml script. In the configuration file, the user specifies which evaluations will be performed, i.e. which measurements to use and which comparisons to make.

Links to example yaml files can be found here:

- Example 1


How MELODIES-MONET works
------------------------

A single driver file reads the yaml configuration file and handles sending all the appropriate calls to the processing code. The driver functions are initiated from a python notebook. The analysis system:

 - transforms model and observations into standard formats
 - includes a standard set of observational data
 - automatically generates diagnostics and plots
 - is easily altered to produce publication-quality plots

Processing requires input to be in netCDF format. MONETIO can be used to convert input files into netCDF prior to using MELODIES-MONET. Data needs to have latitude, longitude, and a datetime object defined.

.. figure:: /_static/MM_classes_connections.png
  :alt: diagram showing the Python classes defined and used in MELODIES-MONET.
  
  Schematic of the classes defined in the Python code and used by MELODIES-MONET.
