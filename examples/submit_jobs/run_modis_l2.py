import os
import sys
sys.path.append('../../')
import driver
from util.pair_obs import pair_obs

an = driver.analysis()
an.control = '../yaml/control_modis_l2.yaml'
an.read_control()
an.open_obs()

# to be added to the analysis class as an.pair_obs()
pair_obs(an)
