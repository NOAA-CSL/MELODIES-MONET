#!/usr/bin/env python
# coding: utf-8

# # MELODIES-MONET-analysis dev 
# First lets just import the driver to see how it happens and so we can play around with it a little 

#This is needed to tell matplotlib to use a non-interactive backend and avoid display errors.
import matplotlib
matplotlib.use('Agg')
from melodies_monet import driver
import os
import dask
import sys

do_stats=sys.argv[1]
#print(sys.argv[1])
#print(do_stats)
#do_stats=0

# ### Driver class
# Now lets create an instance of the python driver analysis class. It consists of 4 main parts; model instances, observation instances, a paired instance of both.  This will allow us to move things around the plotting function for spatial and overlays and more complex plots.

an = driver.analysis()

# ### Control File
# 
# set the yaml control file and begin by reading the file

an.control = 'control.yaml'        # 'control.yaml.rapchemtest'
an.read_control() # control='control.yaml')

dask.config.set(**{'array.slicing.split_large_chunks': True})

# ### Loading the model data 
# 
# driver will automatically loop through the "models" found in the model section of the yaml file and create an instance of the driver.model class for each that includes the label, mapping information, and xarray object as well as the filenames.  Note it can open multiple files easily by including hot keys 

an.open_models()

#All the info in the analysis class can also be called.
print(an.start_time)
print(an.end_time)

an.open_obs()

#All the info in the observation class can also be called.
an.obs['airnow'].obj

#This just pairs the data
print("pairing")
an.pair_data()

print('plotting')
#And this generates all the plots.
an.plotting()

# And the stats
if(do_stats=='1'):
  print('doing stats')
  an.stats()
