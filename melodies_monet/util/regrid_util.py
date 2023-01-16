# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#

"""
file: regrid_util.py
"""
import os
import xarray as xr
import xesmf as xe


def setup_obs_regridder(config):
    """
    Setup regridder for observations

    Parameters
        config (dict): configuration dictionary

    Returns
        regridder (xe.Regridder): regridder instance
        ds_target (xr.Dataset): target grid dataset
    """

    base_file = os.path.expandvars(config['obs']['regrid']['base_grid'])
    ds_base = xr.open_dataset(base_file)
    target_file = os.path.expandvars(config['obs']['regrid']['target_grid'])
    ds_target = xr.open_dataset(target_file)

    regridder = xe.Regridder(ds_base, ds_target, config['obs']['regrid']['method'])

    return regridder, ds_target
