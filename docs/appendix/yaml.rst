Description of All YAML Options
===============================

General Rules
-------------

* Any key that is specific for a plot type will begin with one of the 
  following descriptors:
  - ts for Timeseries
  - ty for Taylor
* When a key is optional it will be followed by #Opt 
* All plots use data over the entire analysis window from the "start_time" 
  to the "end_time" specified in the "analysis" section. 
  - timeseries - average over window provided by "ts_avg_window"
  - taylor - calculated over entire analysis window
  - spatial_bias - average over entire analysis window
  - spatial_overlay - average over entire analysis window
  - boxplot - calculated over entire analysis window
* If "set_axis" = True in "data_proc" section of each "plot_grp", the y-axis 
  for that plot_grp will be set based on the values specified in the "obs" 
  section for each "variable". If "set_axis" = False, then the automatic
  scaling in Matplotlib will be used. 'vmin_plot' and 'vmax_plot' are needed
  for 'timeseries', 'spatial_overlay', and 'boxplot'. 'vdiff_plot' is needed
  for 'spatial_bias' plots and 'ty_scale' is needed for 'taylor' plots. 
  'nlevels' or the number of levels used in the contour plot can also 
  optionally be provided for 'spatial_overlay' plot. If set_axis = True and 
  the proper limits are not provided in the 'obs' section, a warning will 
  print, and the plot will be created using the automatic scaling in
  Matplotlib.

Analysis
--------
All input related to the analysis class.

**start_time:** The start time in UTC of the analysis window.
(e.g., '2019-08-02-12:00:00')

**end_time:** The end time in UTC of the analysis window.
(e.g., '2019-08-03-12:00:00')

**output_dir**: This is optional. This is the directory where the plots are saved. 
If not specified plots stored in code directory. 

**debug:** This is an option to print out plots and more options for trouble 
shooting. If you want plots to print in jupyter notebooks select this to True.
Set this to False, when you are submitting MELODIES MONET as a job to an HPC
machine to avoid display errors. 

Models
------
All input for each instance of the model class. First level should be the model 
label. Then under each model label provide the following:

**files:** The file directory location and name(s). Hotkeys are allowed.

**files_vert:** This is for CMAQ only. If you want to calculate vertical info, 
please provide location of ``*.metcro3d.ncf`` files here.

**files_surf:** This is for CMAQ only. If you want to calculate vertical info, 
please provide location of ``*.metcro2d.ncf`` files here.

**mod_type:** The model type. Options are: 'cmaq', 'wrfchem', 'rrfs', 'gsdchem'. 
If you specify another name, MELODIES MONET will try to read in the data using
xarray.open_mfdataset and xarray.open_dataset().

**mod_kwargs**: This is an optional dictionary to include information to 
provide to the model dataset reader scripts in MONETIO (temporarily in the 
MELODIES-MONET/melodies_monet/new_monetio folder on GitHub). For example, you
can provide mechanism information (e.g., mech: 'cb6r3_ae6_aq') or for some models, 
in order to reduce processing time you can only pull in the surface data 
(e.g., surf_only: True).

**radius_of_influence:** The "radius of influence" used for pairing in MONET. 
Typically this is set at the horizontal resolution of your model * 1.5. Setting 
this to a smaller value will speed up the pairing process. 

**mapping:** This is the mapping dictionary for all variables to be plotted. 
For each observational dataset, add a mapping dictionary where the model 
variable name is first (i.e., key) and the observation variable name is second 
(i.e., value). Because the plots in MELODIES MONET will plot multiple models 
with one observation, the observation variables listed in the mapping dictionary 
must be consistent across all models. For example, if you want to plot the 
results of multiple model datasets against the AirNow observations for "OZONE" 
and "PM2.5", you must provide the model variable names for "OZONE" and "PM2.5" 
in the mapping dictionary for all models. Say if you only provide the model 
variable names for "OZONE" for one of the models, MELODIES MONET will error. Be 
careful that if variable names like NO are a command in python to add 'NO' to 
indicate that it should be interpreted as a string.

For example, ::

  mapping:
    airnow:
      CO: 'CO'
      NO2: 'NO2'
      'NO': 'NO' 
      PM25_TOT: 'PM2.5'
      O3: 'OZONE'
    
**projection:** Not used currently in the code. This is under development and 
likely to be moved to the "plots" section

**plot_kwargs:** This is optional. If you do not provide this, MELODIES MONET 
will use a default list of colors. Add a dictionary of plotting characteristics
to be read in by Matplotlib. 

For example, ::

  plot_kwargs: #Opt
    color: 'magenta'
    marker: 'o'
    linestyle: '--'
  
Copy that above and update the model label for all the models you would like 
to include in the analysis.

Observations
------------
All input for each instance of the observation class. First level should be the 
observation label. Then under each observation label provide the following:

**use_airnow:** If the observations are AirNow set to True, else set to False. 
Generalizing this to include other surface observations is under development.

**filename:**  The file directory location and name. These observations need 
to be preprocessed prior to incorporating them into MELODIES MONET. See 
:doc:`../tutorial/downloading_obs` for more details.

**obs_type:** The observation type. Options are: "pt_sfc" or point surface. Adding 
options for Aircraft and Satellite observations are under development.

**variables:** This is all optional. For each observational variable you can 
include the following information to handle unit conversions, min/max values, 
NaNs, and add optional plotting information. The obs_min, obs_max, and 
nan_values are set to NaN first and then the unit conversion is applied.

   * **unit_scale:** The value for unit conversion.
   * **unit_scale_method:** The method for unit conversion. Options are: 
     Multiply = '*' , Add = '+', subtract = '-', divide = '/'. 
   * **obs_min:** Set all values less than this value to NaN
   * **obs_max:** Set all values greater than this value to NaN
   * **nan_value:** -1.0 # Set this value to NaN
   * **ylabel_plot:** String to use as ylabel in plot. Useful for adding units
     or instrument information.
   * **ty_scale:** Scaling to be used in Taylor plots. 
   * **vmin_plot:** Minimum for y-axis during plotting. To apply to a plot, 
     change set_axis = True in plot_group.
   * **vmax_plot:** Maximum for y-axis during plotting. To apply to a plot, 
     change set_axis = True in plot_group.
   * **vdiff_plot:** The range (+/-) to use in bias plots. To apply to a 
     plot, change set_axis = True in plot_group.
   * **nlevels_plot:** The number of levels used in colorbar for contourf plot. To 
     apply to a plot, change set_axis = True in plot_group.

For example, ::

  PM2.5:
    unit_scale: 1
    unit_scale_method: '*'
    obs_min: 0 
    obs_max: 100
    nan_value: -1.0
    ylabel_plot: 'PM2.5 (ug/m3)'
    ty_scale: 2.0 
    vmin_plot: 0.0 
    vmax_plot: 22.0 
    vdiff_plot: 15.0 
    nlevels_plot: 23

Copy that above and update the observation label for all the observations you 
would like to include in the analysis. Note that all models are paired with all 
observations. At this point MELODIES MONET does not pair observations with each 
other.

Plots
-----
All input for each plotting group. A plotting group consists of one plotting 
type. The plotting types are described in 
:doc:`/background/supported_plots`. All model /
observational pairs and domains specified for the plotting group will be 
included. You may include as many plotting groups as you like.

For each plotting group, update the label and include the following information.
Note: the labels need to be unique, but otherwise are not used.

**type:** The model type. Options are: 'timeseries', 'taylor', 'spatial_bias',
'spatial_overlay', 'boxplot'

**fig_kwargs:** This is optional to provide a dictionary with figure 
characteristics to be read in by Matplotlib. 

For example, ::

  fig_kwargs:
    figsize: [14,6]

**default_plot_kwargs:** This is optional to provide a dictionary with plotting 
characteristics to be read in by Matplotlib. Note that the "plot_kwargs" in the 
"model" section will overwrite these. This is a good method to set the line width 
and marker size for the plot.

For example, ::

  default_plot_kwargs:
    linewidth: 2.0
    markersize: 2.

**text_kwargs:** This is optional to provide a dictionary with text 
characteristics to be read in by Matplotlib.

For example, ::

  text_kwargs:
    fontsize: 18.

**domain_type:** List of domain types to be plotted. These correspond with
the columns in the observation file. (e.g., airnow: epa_region, state_name, 
siteid, etc.).

**domain_name:** List of domain names to be plotted. If domain_type = all, all 
data will be used and the domain_name is used only in the plot title. If 
domain_type is not equal to all, MELODIES MONET will query all of the data 
where domain_type is equal to domain_name.

**data:** This a list of model / observation pairs to be plotted where the 
observation label is first and the model label is second 
(e.g., ['airnow_cmaq_expt', 'airnow_rrfs_13km', 'airnow_wrfchem_v4.2'])

**data_proc:** This section stores all of the data processing information.

   * **rem_obs_nan:** If True, remove all points where model or obs variable is 
     NaN. If False, remove only points where model variable is NaN.
   * **set_axis:** If = True, use the axis constraints described in the 
     observation class (e.g., ty_scale, vmin_plot, vmax_plot, vdiff_plot, 
     nlevels_plot). If = False, use automatic scaling in matplotlib.
   * **ts_select_time:** This is for timeseries plots only. This is the time 
     used for averaging and plotting. Options are 'time' for UTC or 'time_local' 
     for local time
   * **ts_avg_window:** This is for timeseries plots only. This is the averaging 
     window applied to the data. Options are None for no averaging or a pandas 
     resample rule (e.g., 'H' is hourly, 'D' is daily).
   
Stats
-----
All input needed to calculate the statistics. The supported statistics available 
in MELODIES MONET are described in 
:doc:`/background/supported_stats`. All model /
observational pairs and domains specified will be included. You may include as 
many statistics as you like. Note however that the calculation of the statistics 
is relatively slow right now. Optimizing this code is under development.

The statistics require positive numbers, so if you want to calculate temperature 
use Kelvin. Wind direction has special calculations for AirNow if the observation 
name is 'WD'. 

**stat_list:** List of acronyms of statistics to calculate as defined in 
:doc:`/background/supported_stats`. (e.g., ['MB', 'MdnB',
'NMB', 'NMdnB','R2', 'RMSE']). A dictionary of definitions is also included in 
MELODIES-MONET/melodies_monet/stats/proc_stats.py. 

**round_output:** This is optional. This is the integer provided to Pandas 
round function defining the number of decimal places to which to round each 
value. Defaults to 3 (i.e., rounds to 3rd decimal place).

**output_table:** This is optional. The statistics will always output a table in 
.csv format. If True, a matplotlib table figure is also output.

**output_table_kwargs:** This is optional. This is a dictionary defining all
of the characteristics of the matplotlib table figure. This is completely 
customizable because optimal sizes will depend on the number of pairs and 
statistics included.

For example, ::

  output_table_kwargs:
    figsize: [7, 3]
    fontsize: 12.
    xscale: 1.4
    yscale: 1.4
    edges: 'horizontal'


**domain_type:** List of domain types to be plotted. These correspond with
the columns in the observation file. (e.g., airnow: epa_region, state_name, 
siteid, etc.).

**domain_name:** List of domain names to be plotted. If domain_type = all, all 
data will be used and the domain_name is used only in the plot title. If 
domain_type is not equal to all, MELODIES MONET will query all of the data 
where domain_type is equal to domain_name.

**data:** This a list of model / observation pairs to be plotted where the 
observation label is first and the model label is second 
(e.g., ['airnow_cmaq_expt', 'airnow_rrfs_13km', 'airnow_wrfchem_v4.2'])
