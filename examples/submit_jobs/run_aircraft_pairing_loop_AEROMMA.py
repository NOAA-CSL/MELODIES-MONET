# This code uses MELODIES-MONET to read in a .yaml file 
# and produces sets of paired aircraft observations. 
# For an interactive script see jupyter notebooks in main directory.

from melodies_monet import driver
from melodies_monet.util.tools import loop_pairing
import os
import dask
an = driver.analysis()
# -- Update the yaml file below
control_fn = 'control_aircraft_looping_AEROMMA_UFSAQM-submit.yaml'
file_pairs_yaml='supplementary_aircraft_looping_file_pairs_AEROMMA_UFSAQM-submit.yaml'
an.control = control_fn
an.read_control()

# -- Lines below make a copy of the namelist in the plot directory for reference later
cmd = 'cp ' + an.control + ' ' + an.control_dict['analysis']['output_dir']
os.system(cmd)

cmd = 'cp ' + file_pairs_yaml + ' ' + an.control_dict['analysis']['output_dir']
os.system(cmd)

dask.config.set(**{'array.slicing.split_large_chunks': True})
loop_pairing(control=control_fn,file_pairs_yaml=file_pairs_yaml)

