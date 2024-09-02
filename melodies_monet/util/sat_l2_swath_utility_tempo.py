# Copyright (C) 2024 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#

# read all swath data for the time range
# developed for TEMPO Level2 NO2
#

import collections
import logging
import warnings
from datetime import datetime

import monet
import numba
import numpy as np
import xarray as xr
import xesmf as xe

numba_logger = logging.getLogger("numba")
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

    swathdef = monet.util.interp_util.lonlat_to_swathdefinition(
        obsobj["lon"].values, obsobj["lat"].values
    )
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
    swathobj : xr.Dataset
        It contains the referencetime, lat and lon data of the swath
    modobj : xr.Dataset
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
    dp : xr.DataArray
        Pressure difference in layer
    """

    dp_vals = (
        obsobj["pressure"].isel(swt_level_stagg=slice(1, None)).values
        - obsobj["pressure"].isel(swt_level_stagg=slice(1, None)).values
    )
    dp = xr.DataArray(
        data=dp_vals,
        dims=("x", "y", "swt_level"),
        coords={
            "lon": (("x", "y"), obsobj["lon"].values),
            "lat": (("x", "y"), obsobj["lat"].values),
        },
        attrs={
            "units": "Pa",
            "description": "Delta pressure in layer",
            "long_name": "delta_p",
        },
    )
    return dp


@numba.jit(nopython=True)
def _fast_interp_vert(orig, target, data):
    """Performs the numpy interpolation. It is separated from other functions
    for the sake of using the numba jit.

    Parameters:
    -----------
    orig : np.ndarray
        Original grid from which to interpolate. The expected dimensions are (time, z, x, y),
        in that order. The horizontal and time dimensions are expected to be previously
        interpolated. The original pressure levels should be in decreasing order.
    target : np.ndarray
        Target data with vertical grid information. The expected dimensions are (time, z, x, y),
        in that order. The target pressure layers should be in decreasing order.
    data : np.ndarray
        Data to be interpolated. It should have the same grid (including vertical) and dimensions
        as orig.

    Returns
    -------
    interp : np.ndarray
        Interpolated data
    """
    assert orig.shape == data.shape, "Grid shape does not match data"
    nt, nz, nx, ny = target.shape
    interp = np.zeros((nt, nz, nx, ny))
    for t in nt:
        for x in nx:
            for y in ny:
                interp[t, :, x, y] = np.flip(
                    np.interp(
                        np.flip(target[t, :, x, y]),
                        np.flip(orig[t, :, x, y]),
                        np.flip(data[t, :, x, y]),
                    )
                )
    return interp


def interp_vertical_mod2swath(obsobj, modobj, vars=["no2_col"]):
    """Interpolates model vertical layers to TEMPO vertical layers

    Paramenters
    -----------
    modobj : xr.Dataset
        Model data (as provided by MONETIO)
    obsobj : xr.Dataset
        TEMPO data (as provided by MONETIO). Must include pressure.
    vars : list[str]
        Variables to interpolate.

    Returns
    -------
    modsatlayers : xr.Dataset
        Model data (interpolated to TEMPO vertical layers
    """
    assert modobj["longitude"].fillna(0).values == obsobj["lon"].fillna(0).values
    assert modobj["latitude"].fillna(0).values == obsobj["lat"].fillna(0).values

    modsatlayers = xr.Dataset()
    p_mid_tempo = (
        obsobj["pressure"].isel(swt_level_stagg=slice(1, None)).values
        - obsobj["pressure"].isel(swt_level_stagg=slice(None, -1)).values
    )
    p_orig = modobj["pres_pa_mid"]
    dimensions = ("time", "z", "lon", "lat")
    coords = {
        "time": (("time",), modobj["time"].values),
        "lon": (("x", "y"), modobj["longitude"].values),
        "lat": (("x", "y"), modobj["latitude"].values),
    }
    for var in vars:
        interpolated = _interp_mod2swath(p_orig, p_mid_tempo, modobj[var].values)
        modsatlayers[var] = xr.DataArray(
            data=interpolated, dims=dimensions, coords=coords, attrs=modobj[var].attrs
        )
    modsatlayers["p_mid_tempo"] = xr.DataArray(
        data=p_mid_tempo,
        dims=dimensions,
        coords=coords,
        attrs=modobj["p_mid_tempo"].attrs,
    )
    _interp_description = "Mid layer pressure interpolated to TEMPO mid swt_layer pressures"
    modsatlayers["p_mid_tempo"].attrs["description"] = _interp_description
    return modsatlayers


def _calc_partialcolumn(modobj, var='NO2'):
    """Calculates the partial column of a species from its concentration.

    Parameters
    ----------
    modobj : xr.Dataset
        Model data

    Returns
    -------
    partial_col : xr.DataArray
        DataArray containing the partial column of the species.
    """
    R = 8.314 # m3 * Pa / K / mol
    layer_thickness = _calc_layer_thickness(modobj)
    partial_col = (modobj["pres_pa_mid"] * modobj[var] * layer_thickness) / (R * modobj["tk"])

    return partial_col

def _calc_layer_thickness(modobj):
    """Calculates layer thickness

    Parameters
    ----------
    modbj : xr.Dataset
        Model data as produced by the MONETIO reader.

    Returns
    -------
    layer_thickness : xr.DataArray
        Layer thickness in m.
    """
    layer_thickness = xr.zeros_like(modobj["height_agl"])
    layer_thickness[:, 0, :, :] = modobj["height_agl"].isel(z=0).values
    layer_thickness[:, 1:, :, :] = (
        modobj["height_agl"].isel(z=slice(2, None)).values
        - modobj["height_agl"].isel(z=slice(1, -1)).values
    )
    return layer_thickness
