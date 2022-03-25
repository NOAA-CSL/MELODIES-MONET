# This code uses MELODIES-MONET to read in a .yaml file 
# and produces plots. For an interactive script see 
# jupyter notebooks in main directory.

#This is needed to tell matplotlib to use a non-interactive backend and avoid display errors.
import matplotlib
matplotlib.use('Agg')
import sys; sys.path.append("../../")
from melodies_monet import driver
import os
import dask
an = driver.analysis()
# -- Update the yaml file below
an.control = '../yaml/control_cmaq-rrfs_surface-all-short_test_jupyter.yaml'
an.read_control()
# -- Lines below make a copy of the namelist in the plot directory for reference later
cmd = 'cp ' + an.control + ' ' + an.control_dict['analysis']['output_dir']
os.system(cmd)
dask.config.set(**{'array.slicing.split_large_chunks': True})
an.open_models()
an.open_obs()
an.pair_data()
# -- Optionally comment out lines below if only want to create plots or stats
an.plotting()
an.stats()
