# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#

"""
file: regrid_util.py
"""
import os
import xarray as xr


def setup_regridder(config, config_group='obs', target_grid=None):
    """
    Setup regridder for observations or model

    Parameters
        config (dict): configuration dictionary

    Returns
        regridder (dict of xe.Regridder): dictionary of regridder instances
    """
    try:
        import xesmf as xe
    except ImportError as e:
        print('regrid_util: xesmf module not found')
        raise

    print('setup_regridder.target_grid')
    print(target_grid)

    if target_grid is not None:
        ds_target = target_grid
    else:
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


def filename_regrid(filename, regridder):
    """
    Construct modified filename for regridded dataset

    Parameters
        filename (str): filename of dataset
        regridder (xe.Regridder): regridder instance

    Returns
        filename_regrid (str): filename of regridded dataset
    """
    filename_regrid = filename.replace('.nc', '_regrid.nc')

    return filename_regrid

