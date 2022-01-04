Description
===========

Aims of MELODIES-MONET
----------------------
MELODIES-MONET aims to help researchers systematically and reproducibly compare 
global and regional atmospheric chemistry model output with measurements and/or 
other models. It aims to overcome several challenges, including providing access 
to the complex atmospheric chemistry datasets (numerous compounds with complex, 
non-standardized names, various time and spatial sampling, different instruments 
and platforms), as well as routines to extract model results at appropriate time 
and spatial resolutions to quantitatively compare to the observations. While it 
is not an exhaustive analysis system, MELODIES-MONET provides examples of 
analysis and tools to make accurate, easy comparisons. The software system 
can be added to with additional models and measurements. Our goals are:

- well documented for easy use by the community
- written to be user-friendly with easy customization
- written in open-source programming languages, mainly Python 
- made available through GitHub: https://github.com/NOAA-CSL/MELODIES-MONET

Connection with NOAA's MONET
----------------------------

At the same time that a need for MELODIES was identified, NOAA was developing 
a tool to be used for the automated comparison of model forecasts and 
measurements at ground stations in the USA. Thus, the Model and ObservatioN 
Evaluation Toolkit (MONET) was created (Baker and Pan, 2017). The usefulness 
for the broader community was quickly realized and the clear overlap between 
the two projects allowed collaboration between NCAR and NOAA and the creation 
of the MELODIES-MONET project.

MONET specifically encompasses two packages:

- MONET I/O package that focuses on importing and exporting model and observational 
  data into the necessary formats for further analysis (https://monetio.readthedocs.io),  and 
- MONET base package that sets up universal functions such as re-gridding and nearest 
  neighbor calculations (https://monet-arl.readthedocs.io/).

The two MONET packages above are used by MELODIES-MONET to complete analysis 
and plotting. Most users will only need to work with the MELODIES-MONET 
code-set to complete their desired analysis.

.. figure:: /_static/MM_diagram_connection.png
  :alt: diagram showing connections between MONET, MONETIO and MELODIES-MONET
  
  The structure and function of MONET and MONET I/O, and MELODIES-MONET interaction. 
  Adapted and updated from Figure 1 in Baker and Pan (2017).

