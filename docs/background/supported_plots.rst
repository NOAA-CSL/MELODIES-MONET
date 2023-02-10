Supported Plots
===============

Model to Model Comparisons
--------------------------
Under development. 

See the `Spatial Verification <https://github.com/NOAA-CSL/MELODIES-MONET/projects/3>`_ 
project on GitHub to learn about current and future development.

Model to Observation Comparisons
---------------------------------

Surface Evaluation
^^^^^^^^^^^^^^^^^^
.. figure:: /_static/figures/plot_grp1.timeseries.PM2.5.2019-08-01_12.2019-08-11_12.all.CONUS.png

   **Timeseries** - Plot comparing one or more model results with one
   observation (y-axis) versus time (x-axis) over the analysis window with
   options to specify the domain, time (local or UTC), and averaging window.

.. figure:: /_static/figures/plot_grp2.taylor.OZONE.2019-08-01_12.2019-08-11_12.all.CONUS.png
   :scale: 25 %

   **Taylor** - Taylor plot comparing one or more model results with one
   observation over the analysis window with options to specify the domain.      
     
.. figure:: /_static/figures/plot_grp3.spatial_bias.OZONE.2019-08-01_12.2019-08-11_12.all.CONUS.airnow_rrfs_diurnal_fires.png

   **Spatial Bias** - Difference plot of model - observations averaged over
   the analysis window with options to specify the domain. Defaults to average,
   but users can also optionally plot percentiles.

.. figure:: /_static/figures/plot_grp6.spatial_bias_exceedance.OZONE_reg.2019-08-01_13.2019-08-31_12.all.CONUS.airnow_cmaq_oper_exceedance.png

   **Spatial Bias Exceedance** - Difference plot of model - observations for the number of
   exceedances greater than the regulatory standard within the analysis window with options to specify
   the domain. This only works for regulatory calculations (regulatory = True) for variables "OZONE" and "PM2.5".
   An exceedance occurs when MDA8 ozone is greater than 70 ppbv or 24 hour averaged PM\ :sub:`2.5` \ is
   greater than 35 ug m\ :sup:`-3`\.
     
.. figure:: /_static/figures/plot_grp4.spatial_overlay.OZONE.2019-08-01_12.2019-08-11_12.all.CONUS.airnow_rrfs_diurnal_fires.png

   **Spatial Overlay** - Model results in contours with observational
   results overlaid in markers averaged over the analysis window with
   options to specify the domain.  
  
.. figure:: /_static/figures/plot_grp5.boxplot.OZONE.2019-08-01_12.2019-08-11_12.all.CONUS.png
   :scale: 25 %

   **BOX-plot** - BOX-plot comparing one or more model results with one
   observation over the analysis window with options to specify the domain.

See the `Expand Surface Observations <https://github.com/NOAA-CSL/MELODIES-MONET/projects/7>`_ 
project on GitHub to learn about current and future development.

Aircraft Evaluation 
^^^^^^^^^^^^^^^^^^^
Under development. 

See the `Incorporate Aircraft Evaluation <https://github.com/NOAA-CSL/MELODIES-MONET/projects/4>`_ 
project on GitHub to learn about current and future development.

Satellite Evaluation 
^^^^^^^^^^^^^^^^^^^^
Under development.

See the `Incorporate Satellite Evaluation <https://github.com/NOAA-CSL/MELODIES-MONET/projects/6>`_ 
project on GitHub to learn about current and future development.
