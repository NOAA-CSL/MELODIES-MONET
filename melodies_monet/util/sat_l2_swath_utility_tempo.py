# Copyright (C) 2024 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#

# read all swath data for the time range
# developed for TEMPO Level2 NO2
#

import xesmf as xe
import numpy as np
import xarray as xr
from datetime import datetime
import monet
import collections

import logging
numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)


def create_swathdefinition(obsobj):
    """Creates a pyresample SwathDefinition with the TEMPO swath grid

    Partameters
    -----------
    obsobj : xr.Dataset 
         Object containing lat andsdf:q
         lon of satellite observations

    Returns
    -------
    pyresample.geometry.SwathDefinition
        A pyresample SwathDefinition with the swath data
    """

    swathdef = monet.util.interp_util.lonlat_to_swathdefinition(obsobj["lon"].values, obsobj["lat"].values)
    return swathdef

def tempo_interp_mod2swath(obsobj, modobj, method="bilinear"):
    """Interpolate model to satellite swath/swaths

    Parameters
    ----------
    obsobj : xr.Dataset | collections.OrderedDict
        satellite with swath data. If type is xr.Dataset, a single 
        swath is assumed. If type is collections.OrderedDict, a key containing
        each swath in a scan is assumed.
    modobj : xr.Dataset
        model data (with no2 col calculated)

    Returns
    -------
    modswath : xr.Dataset | collections.OrderedDict
        Regridded model data at swath or swaths. If type is xr.Dataset, a single 
        swath is returned. If type is collections.OrderedDict, it returns an
        OrderedDict in which each time represents the reference time of the swath.
    """

    if obsobj.isinstance(collections.OrderedDict):
        modswath = collections.OrderedDict()
        for key in obsobj.keys():
            modswath[key] = _interp_mod2swath(obsobj[key], modobj, method=method)
    elif obsobj.isinstance(xr.Dataset):
        modswath = _interp_mod2swath(obsobj, modobj, method)
    return modswath


def _interp_mod2swath(swathobj, modobj, method="bilinear"):
    """Interpolate model to swath,

    Parameters
    ----------
    swathobj: xr.Dataset
        It contains the referencetime, lat and lon data of the swath
    modobj: xr.Dataset
        It contains the data to be regridded (time, lat and lon)

    Returns
    -------
    modswath : xr.Dataset
        Dataset containing the regridded model
    """

    mod_at_swathtime = modobj.interp(time=swathobj.time.mean())
    regridder = xe.Regridder(mod_at_swathtime, swathobj, method, unmapped_to_nan=True)
    modswath = regridder(modobj)
    return modswath

def _calc_dp(obsobj):
    """Calculate delta pressure in satellite layers

    Parameters
    ----------
    obsobj : xr.Dataset
        satellite observations containing pressure (in Pa)

    Returns
    -------
    dp: xr.DataArray
        Pressure difference in layer
    """

    dp_vals = (obsobj["pressure"].isel(swt_level_stagg=slice(1, None)).values 
               - obsobj["pressure"].isel(swt_level_stagg=slice(1, None)).values)
    dp = xr.DataArray(
        data=dp_vals,
        dims = ("x", "y", "swt_level"),
        coords = {"lon": (("x", "y"), obsobj["lon"].values),
                  "lat": (("x", "y"), obsobj["lat"].values)},
        attrs = {"units": "Pa",
                 "description": "Delta pressure in layer",
                 "long_name": "delta_p"}
    )
    return dp

