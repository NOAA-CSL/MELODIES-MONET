# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#

"""
file: regrid_util.py
"""
import os
import xarray as xr
import xesmf as xe


def setup_regridder(config, config_type='obs'):
    """
    Setup regridder for observations or model

    Parameters
        config (dict): configuration dictionary

    Returns
        regridder (xe.Regridder): regridder instance
        ds_target (xr.Dataset): target grid dataset
    """

    base_file = os.path.expandvars(config[config_type]['regrid']['base_grid'])
    config_typeds_base = xr.open_dataset(base_file)
    target_file = os.path.expandvars(config[config_type]['regrid']['target_grid'])
    ds_target = xr.open_dataset(target_file)

    regridder = xe.Regridder(ds_base, ds_target, config[config_type]['regrid']['method'])

    return regridder, ds_target

