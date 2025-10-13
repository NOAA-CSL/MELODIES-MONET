# SPDX-License-Identifier: Apache-2.0
#

# read all swath data for the time range
# developed for TEMPO Level2 NO2
#

"""Python utility for TROPOMI use."""

import warnings
import logging

import numba
import numpy as np
import pandas as pd
import xarray as xr
import xesmf as xe

from .sat_l2_swath_utility_tempo import (  # calc_grid_corners,
    calc_altitude_from_thickness,
    calc_dz_m_from_altitude,
)
from .satellite_utilities import mod_to_overpasstime
from .tools import calc_partialcolumn

# import warnings


numba_logger = logging.getLogger("numba")
numba_logger.setLevel(logging.WARNING)


default_ak_variable_names = {
    "tropomi_l2_no2": {
        "averaging_kernel": "averaging_kernel",
        "tropospheric_averaging_kernel_calc": True,
        "airmass_factor_total": "air_mass_factor_total",
        "airmass_factor_troposphere": "air_mass_factor_troposphere",
    },
    "tropomi_l2_hcho": {
        "averaging_kernel": "averaging_kernel",
        "tropospheric_averaging_kernel_calc": False,
        "airmass_factor_total": "formaldehyde_clear_air_mass_factor",
        "airmass_factor_troposphere": "formaldehyde_tropospheric_air_mass_factor",
    },
}

default_mod_variable_names = {
    "tropomi_l2_no2": "NO2",
    "tropomi_l2_hcho": "HCHO",
    "tropomi_l2_co": "CO",
}


def interp_horizontal_mod2sat(obsobj, modobj, method="bilinear", isglobal=False, **kwargs):
    """Interpolates model horizontally to satellite

    Parameters
    ----------
    obsobj : xr.Dataset
        xr.Dataset with a granule of the obsobj.
    modobj : xr.Dataset
        Dataset with satellite data as formatted by monetio, already
        interpolated to correct time.
    method : str
        Method of regridding, any method supported by xesmf should work.
    isglobal : bool
        Whether the model is global. If True, xe.Regridder will be set
        to periodic=True
    **kwargs
        Extra arguments to pass to xe.Regridder

    Returns
    -------
    dict[np.datetime64, xr.Dataset]
        Dictionary with model data in obsobj grid.
    """

    regridder = xe.Regridder(
        modobj,
        obsobj,
        ignore_degenerate=True,
        unmapped_to_nan=True,
        method=method,
        periodic=isglobal,
        **kwargs,
    )
    return regridder(modobj)


def mod_to_overpasstime(modobj, opass_tms):
    """
    Interpolate model to satellite overpass time.

    Parameters
    ----------

    modobj : xarray.Dataset
        model data
    opass_tms : pandas.DatetimeIndex
        satellite overpass local time
    partial_col : str
        variable to calculate partial columns for
    Output
    ------
    outmod : xarray.Dataset
        revised model data at local overpass time
    """

    (nst,) = opass_tms.shape
    # nmt, = modobj.time.shape
    # ny,nx = modobj.longitude.shape

    # Determine local time offset
    local_utc_offset = (modobj["longitude"] / 15).round(1).astype("timedelta64[h]")
    # initialize local time as variable
    modobj["localtime"] = modobj["time"] + local_utc_offset

    # initialize new model object with satellite datetimes
    outmod = []

    for ti in np.arange(nst):
        # Apply filter to select model data within +/- 1 output time step of the overpass time
        tempmod = modobj.where(
            np.abs(modobj["localtime"] - opass_tms[ti].to_datetime64())
            < (modobj.time[1] - modobj.time[0])
        )

        # determine factors for linear interpolation in time
        tfac = 1 - (
            np.abs(tempmod["localtime"] - opass_tms[ti].to_datetime64())
            / (modobj.time[1] - modobj.time[0])
        )
        tempmod = tempmod.drop_vars("localtime")
        # Carry out time interpolation
        ## Note regarding current behavior: will only carry out time interpolation if at least 2 model timesteps
        outmod.append((tfac * tempmod).sum(dim="time", min_count=2, keep_attrs=True))
    # print(outmod)
    outmod = xr.concat(outmod, dim="time")
    outmod["time_local"] = (["time"], opass_tms)
    outmod.coords["date"] = opass_tms.date
    return outmod


@numba.jit(nopython=True)
def _interpolate_time(orig, target, data):
    pass


def interpolate_time(modelobj, overpass_time=None):
    """Interpolates data in time from orig to target times.

    Parameters
    ----------
    modelobj : xr.Dataset
        Model data to interpolate. It must contain a time dimension
        and a longitude coordinate or variable.
    overpass_time : int | None
        If provided, the overpass time used for the interpolation
        in hours. If None, 13.5 is used.

    Returns
    -------
    xr.Dataset
        Interpolated data.
    """
    overpass = 13.5 if overpass_time is None else overpass_time
    overpass_ns = int(overpass * 3600 * 1e9)
    utc_offset_nanoseconds = (modelobj["longitude"] / 15 * 3600 * 1e9).astype(int)

    days = np.unique(modelobj["time"].dt.floor("D"))
    interpolated_data = []
    for day in days:
        day_min = day - np.timedelta64(1, "D")
        day_max = day + np.timedelta64(1, "D")
        modelobj_day = modelobj.sel(time=slice(day_min, day_max))
        target_time = day + np.timedelta64(overpass_ns, "ns")
        localtime = modelobj_day["time"] + utc_offset_nanoseconds
        previous_index, next_index = _calculate_previous_and_next_indices(localtime, target_time)
        previous_weight, next_weight = _calculate_time_weights(
            localtime, target_time, previous_index, next_index
        )
        previous = modelobj_day.isel(time=previous_index)
        previous["time"] = (("time"), [target_time])
        import pdb; pdb.set_trace()
        interp = (previous_weight * modelobj_day.isel(time=previous_index)) + (
            next_weight * modelobj_day.isel(time=next_index)
        )
        interp_with_time = interp.expand_dims('time', axis=0).assign_coords(time=[target_time])
        import pdb; pdb.set_trace()
        interp.time = (('time'), [target_time])
        interpolated_data.append(interp)
    concat_data = xr.concat(interpolated_data, dim="time")
    return concat_data


def _calculate_previous_and_next_indices(localtime, target_time):
    """Calculates the indices of the time immmediately before and after the target time.
    Parameters
    ----------
    localtime : xr.DataArray
        Local time of the model data.
    target_time : np.datetime64
        Target time to interpolate to.
    Returns
    -------
    tuple[xr.DataArray, xr.DataArray]
        Indices of the previous and next times.
    """
    timediff = (localtime - target_time).astype("float64")
    previous_timediff = timediff.where(timediff <= 0, np.nan)
    next_timediff = timediff.where(timediff > 0, np.nan)
    previous_index = np.abs(previous_timediff).argmin(dim="time", skipna=True)
    next_index = np.abs(next_timediff).argmin(dim="time", skipna=True)
    return previous_index, next_index


def _calculate_time_weights(localtime, target_time, previous_index, next_index):
    """Calculates the weights for the time interpolation.

    Parameters
    ----------
    localtime : xr.DataArray
        Local time of the model data.
    target_time : np.datetime64
        Target time to interpolate to.
    previous_index : xr.DataArray
        Index of the previous time.
    next_index : xr.DataArray
        Index of the next time.

    Returns
    -------
    tuple[xr.DataArray, xr.DataArray]
        Weights for the previous and next times.
    """
    previous_time = localtime.isel(time=previous_index)
    next_time = localtime.isel(time=next_index)
    total_diff = (next_time - previous_time).astype("float64")
    previous_diff = (target_time - previous_time).astype("float64")
    next_diff = (next_time - target_time).astype("float64")
    previous_weight = 1 - previous_diff / total_diff
    next_weight = 1 - next_diff / total_diff
    return previous_weight, next_weight


@numba.jit(nopython=True)
def _interp_vert(orig, target, data):
    """Performs the numpy interpolation. It is separated from other functions
    for the sake of using the numba jit.

    Parameters:
    -----------
    orig : np.ndarray
        Original grid from which to interpolate. The expected dimensions are (z, x, y),
        in that order. The horizontal and time dimensions are expected to be previously
        interpolated. The original pressure levels should be in decreasing order.
    target : np.ndarray
        Target data with vertical grid information. The expected dimensions are (z, x, y),
        in that order. The target pressure layers should be in decreasing order.
    data : np.ndarray
        Data to be interpolated. It should have the same grid (including vertical) and dimensions
        as orig.

    Returns
    -------
    np.ndarray
        Interpolated data
    """
    assert orig.shape == data.shape, "Grid shape does not match data"
    nz, ny, nx = target.shape
    interp = np.zeros((nz, ny, nx))
    for y in range(ny):
        for x in range(nx):
            interp[:, y, x] = np.flip(
                np.interp(
                    np.flip(target[:, y, x]),
                    np.flip(orig[:, y, x]),
                    np.flip(data[:, y, x]),
                )
            )
    return interp


def interp_vertical_mod2swath(obsobj, modobj, variables="NO2_col"):
    """Interpolates model vertical layers to TEMPO vertical layers

    Parameters
    ----------
    modobj : xr.Dataset
        Model data (as provided by MONETIO)
    obsobj : xr.Dataset
        TEMPO data (as provided by the reader). Must include pressure.
    variables : str | list[str]
        Variables to interpolate.

    Returns
    -------
    xr.Dataset
        Model data (interpolated to TEMPO vertical layers
    """
    assert np.all(modobj["longitude"].fillna(0).values == obsobj["longitude"].fillna(0).values)
    assert np.all(modobj["latitude"].fillna(0).values == obsobj["latitude"].fillna(0).values)

    modsatlayers = xr.Dataset()
    p_mid_tropomi = (
        obsobj["pressure"].isel(swt_level_stagg=slice(None, -1)).values
        + obsobj["pressure"].isel(swt_level_stagg=slice(1, None)).values
    ) / 2
    p_orig = modobj["pres_pa_mid"].values
    dimensions = ("z", "x", "y")
    coords = {
        "longitude": (("x", "y"), modobj["longitude"].values),
        "latitude": (("x", "y"), modobj["latitude"].values),
    }
    for var in list(variables):
        interpolated = _interp_vert(p_orig, p_mid_tropomi, modobj[var].values)
        modsatlayers[var] = xr.DataArray(
            data=interpolated, dims=dimensions, coords=coords, attrs=modobj[var].attrs
        )
    modsatlayers["pres_pa_mid"] = xr.DataArray(
        data=p_mid_tropomi,
        dims=dimensions,
        coords=coords,
        attrs=modobj["pres_pa_mid"].attrs,
    )
    _interp_description = "Mid layer pressure interpolated to tropomi mid vertical layer pressures"
    modsatlayers["pres_pa_mid"].attrs["description"] = _interp_description
    return modsatlayers


def apply_averaging_kernel(modobj, obsobj, sat_type, varname=None, averaging_kernel_params=None):
    """Applies the averaging kernel and calculates the column

    Parameters
    ----------
    modobj : xr.Dataset
        DataArray containing the model information. It has to be
        previously regridded to satellite space.
    obsobj : xr.Dataset
        Dataset containing all the observational data, including the
        variables related to the averaging kernel.
    sat_type : str
        string of satellite type. Currently, tropomi_l2_no2,
        tropomi_l2_hcho and tropomi_l2_co are supported
    averaging_kernel_params : dict[str, str]
        dictionary containing the names the keys "averaging_kernel" and
        "tropospheric_averaging_kernel_calc" plus, optionally,
        "airmass_factor_total" and "airmass_factor_troposphere".

    Returns
    -------
    xr.DataArray
        DataArray containing the model columns after applying the averaging kernel.
    """
    if averaging_kernel_params is not None:
        ak_params = {**default_ak_variable_names[sat_type], **averaging_kernel_params}
    else:
        ak_params = default_ak_variable_names[sat_type]
    if ak_params["tropospheric_averaging_kernel_calc"]:
        ak = (
            obsobj[ak_params["airmass_factor_total"]]
            / obsobj[ak_params["airmass_factor_troposphere"]]
            * obsobj[ak_params["averaging_kernel"]]
        )
    else:
        ak = obsobj[ak_params["averaging_kernel"]]
    if "tm5_tropopause_pressure" in obsobj:
        ak = ak.where(obsobj["pres_pa_mid"] >= obsobj["tm5_tropopapuse_pressure"], other=0)
    if varname is None:
        varname = {
            "tropomi_l2_no2": "NO2",
            "tropomi_l2_hcho": "HCHO",
            "tropomi_l2_co": "CO",
        }[sat_type]
    partial_cols = calc_partialcolumn(modobj, varname)

    column_data_model = xr.dot(ak, partial_cols, dim="z")
    column_data_model.attrs["description"] = "column after applying averaging kernel"
    return column_data_model


def crop_obsobj(obsobj, bounds):
    """Discards the observations outside the bounds.

    Parameters
    ----------
    obsobj : xr.Dataset
        Dataset containing all the observational data.
    bounds : list[float]
        List containing the bounds in the order [min_lon, max_lon, min_lat, max_lat].

    Returns
    -------
    xr.Dataset
        Dataset with the observations outside the bounds discarded.
    """
    return obsobj.where(
        (obsobj["longitude"] >= bounds[0])
        & (obsobj["longitude"] <= bounds[1])
        & (obsobj["latitude"] >= bounds[2])
        & (obsobj["latitude"] <= bounds[3]),
        drop=True,
    )


def within_model_domain(obsobj, bounds):
    """Checks if any of the observations are within the model domain.

    Parameters
    ----------
    obsobj : xr.Dataset
        Dataset containing all the observational data.
    bounds : list[float]
        List containing the bounds in the order [min_lon, max_lon, min_lat, max_lat].

    Returns
    -------
    bool
        True if all observations are within the model domain, False otherwise.
    """
    return (
        (
            (obsobj["longitude"] >= bounds[0])
            & (obsobj["longitude"] <= bounds[1])
            & (obsobj["latitude"] >= bounds[2])
            & (obsobj["latitude"] <= bounds[3])
        )
        .any()
        .item()
    )


def calc_local_geodate(time, longitude):
    """Calculates the geographical date based on longitude

    Parameters
    ----------
    time : xr.DataArray
        DataArray containing time for each pixel
    longitude : xr.DataArray
        Longitude for each pixel. It has to be
        [-180; 180]

    Returns
    -------
    xr.DataArray
        DataArray containing the local time based on longitude
        for each granule.
    """
    assert ((longitude <= 180) | longitude.isnull()).all()
    if len(time.shape) == 2:
        # Some of the TROPOMI datasets have delta_time depending only
        # on scanline (e.g. NO2)
        return time + np.timedelta64(longitude.isel(x=0).values * 240, "s")
    if len(time.shape) == 3:
        # Some of the TROPOMI datasets have delta_time depending
        # on scanline and ground_pixel (e.g. HCHO)
        return time + np.timedelta64(longitude.values * 240, "s")
    raise Exception("time variable (e.g., time_granule) has wrong dimension number.")


def _regrid_and_apply_ak(modobj, obsobj):
    """Regrids and applies AK to one swath.

    Parameters
    ----------
    modobj : xr.Dataset
        model dataset, as read in by MELODIES-MONET. Model should be
        already at overpass time
    obsobj : dict[np.datetime64, xr.Dataset]
        Dictionary containing observations

    Returns
    -------
    dict[str, xr.Dataset]
        Dictionary containing the same keys as the obs, and the model
        data in satellite space after applying the ak as values
    """

    output_pair = {}
    overpass_time = pd.date_range(
        pd.to_datetime(modobj["time"][0].values).replace(hour=13, minute=30),
        pd.to_datetime(modobj["time"][-1].values).replace(hour=13, minute=30),
        freq="D",
    )
    modobj_at_overpass_time = mod_to_overpasstime(modobj, overpass_time)
    modobj_at_overpass_time["altitude"] = calc_altitude_from_thickness(
        modobj_at_overpass_time["dz_m"]
    )
    bounds = [
        modobj["longitude"].min().values,
        modobj["longitude"].max().values,
        modobj["latitude"].min().values,
        modobj["latitude"].max().values,
    ]
    # TODO: Accelerate this by using the bounds to select the obsobj
    obsobj_dates = np.unique(obsobj["time_granule"].dt.date.values)
    modobj_dates = modobj_at_overpass_time["time_local"].dt.date.values
    import pdb

    pdb.set_trace()  # noqa: E999
    for d in obsobj_dates:
        if d not in modobj_dates:
            warnings.warn(f"Model does not have data for {d}, skipping.")
            continue
        obsobj_cropped = crop_obsobj(obsobj[d], bounds)
        modobj_regrid = interp_horizontal_mod2sat(
            obsobj_cropped, modobj_at_overpass_time.sel(date=d)
        )
        # Apply vertical interpolation
        modobj_regrid = interp_vertical_mod2swath(obsobj_cropped, modobj_regrid)
        # Apply averaging kernel
        modobj_regrid[default_mod_variable_names["tropomi_l2_no2"]] = apply_averaging_kernel(
            modobj_regrid, obsobj_cropped, "tropomi_l2_no2"
        )
        output_pair[d] = modobj_regrid


if __name__ == "_main__":
    # For debugging purposes
    import read_sat_l2_tropomi as reader

    satellite_data = reader.read_tropomi_l2_no2(
        "/glade/derecho/scratch/plichtig/TROPOMI/S5P_RPRO_L2__NO2____20220720T000819_20220720T014949_24695_03_020400_20230203T042051.nc",
        variable_names={
            "latitude": "latitude",
            "longitude": "longitude",
            "time_granule": "time",
            "pressure": "pressure",
            "NO2_column_number_density": "NO2_column_number_density",
            "quality_flag": "qa_value",
            "cloud_fraction": "cloud_fraction",
            "surface_albedo": "surface_albedo",
            "solar_zenith_angle": "solar_zenith_angle",
            "viewing_zenith_angle": "viewing_zenith_angle",
            "relative_azimuth_angle": "relative_azimuth_angle",
            # Averaging kernel related variables
            "averaging_kernel": "averaging_kernel",
            "airmass_factor_total": "air_mass_factor_total",
            "airmass_factor_troposphere": "air_mass_factor_troposphere",
        },
        time_variable_name="time_granule",
        time_variable_type="datetime64[ns]",
        time_variable_units=None,
        time_variable_calendar=None,
        lat_variable_name="latitude",
        lon_variable_name="longitude",
        data_variable_names=["NO2_column_number_density"],
        qa_variable_name="qa_value",
        cloud_variable_name="cloud_fraction",
        pressure_variable_name="pressure",
        extra_variables=[
            "surface_albedo",
            "solar_zenith_angle",
            "viewing_zenith_angle",
            "relative_azimuth_angle",
            # Averaging kernel related variables
            "averaging_kernel",
            "air_mass_factor_total",
            "air_mass_factor_troposphere",
        ],
    )
    breakpoint()
