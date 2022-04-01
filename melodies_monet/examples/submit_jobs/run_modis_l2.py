import os
import sys
sys.path.append('../../')
import driver

an = driver.analysis()
an.control = '../yaml/control_modis_l2.yaml'
an.read_control()
an.open_obs()

for obs_dataset in an.obs.keys():
    for granule in an.obs[obs_dataset].obj.keys():
        print(an.obs[obs_dataset].obj[granule])
