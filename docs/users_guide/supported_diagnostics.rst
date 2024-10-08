Supported Diagnostics
=====================

Supported data analysis options in MELODIES MONET are explained below.

Calculating Regulatory Metrics
------------------------------

MDA8 (8-hour daily maximum) ozone and 24 hour average PM\ :sub:`2.5` \ can be
calculated within MELODIES MONET from hourly observational and model data. All plots
except for the ``spatial_overlay`` and ``scorecard`` plots and all stats will work with the regulatory
metrics.

The MDA8 ozone metric is calculated based on the following:
  1) Local time is used
  2) For each day, rolling 8-hour averages are calculated for each period with at least 6 hours of data
  3) The MDA8 value is the highest of the 8-hour averages in a given day
  4) The MDA8 value is only used for a given day if at least 18 of the 24 possible 8-hour averages are available

As described in the EPA report on
`Health Risk and Exposure Assessment for Ozone <https://www3.epa.gov/ttn/naaqs/standards/ozone/data/20140131healthrea4a.pdf>`__
(page 8).

The 24 hour average PM\ :sub:`2.5` \ metric is calculated based on the following:
  1) Local time is used
  2) The average PM\ :sub:`2.5` \ value is calculated over a given day
  3) The PM\ :sub:`2.5` \ value is only used for a given day if at least 18 of the 24 possible hours are available

In order to calculate the regulatory metric, add ``regulatory: True`` into the input
YAML file under the "obs" section for each variable that you want to apply the calculation.
Currently, this option only works for "OZONE" and "PM2.5" variables. There are separate
plotting characteristics for the regulatory options ("ylabel_reg_plot", "vmin_reg_plot",
"vmax_reg_plot", and "vdiff_reg_plot") that can also be specified in the input YAML file.
An example input YAML file that calculates MDA8 ozone and 24 hr PM\ :sub:`2.5` \ is in
``examples/yaml/control_rrfs_cmaq_airnow_reg.yaml``. These input YAML file options are
also further described in the Appendix under :doc:`/appendix/yaml`.


