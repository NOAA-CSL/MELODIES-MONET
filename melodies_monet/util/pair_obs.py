"""
file:  pair_obs.py

method:  pair_obs
    to be incorporated into
    class analysis
"""
import numpy as np
import xarray as xr


def pair_obs(an):
    for obs_name in an.obs.keys():
        for granule in an.obs[obs_name].obj.keys():
            print(an.obs[obs_name].obj[granule])
