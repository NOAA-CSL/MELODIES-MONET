# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#

"""
file: regrid_util.py
"""
import os
import xarray as xr
import xesmf as xe


def setup_regridder(config, config_group='obs'):
    """
    Setup regridder for observations or model

    Parameters
        config (dict): configuration dictionary

    Returns
        regridder (dict of xe.Regridder): dictionary of regridder instances
    """
    target_file = os.path.expandvars(config['analysis']['target_grid'])
    ds_target = xr.open_dataset(target_file)

    regridder_dict = dict()

    for name in config[config_group]:
        base_file = os.path.expandvars(config[config_group][name]['regrid']['base_grid'])
        ds_base = xr.open_dataset(base_file)
        method = config[config_group][name]['regrid']['method']
        regridder = xe.Regridder(ds_base, ds_target, method)
        regridder_dict[name] = regridder

    return regridder_dict

