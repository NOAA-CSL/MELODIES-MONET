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
import xarray as xr
import xesmf as xe

from .sat_l2_swath_utility_tempo import (  # calc_grid_corners,
    calc_altitude_from_thickness,
    calc_dz_m_from_altitude,
)
from .tools import calc_partialcolumn, N_A

# import warnings


numba_logger = logging.getLogger("numba")
numba_logger.setLevel(logging.WARNING)


M2TOCM2 = 1e4
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


def tropomi_mol_m2_to_molec_cm2(column_data):
    """Converts column data from mol/m2 to molec/cm2

    Parameters
    ----------
    column_data : xr.DataArray
        DataArray containing the column data in mol/m2

    Returns
    -------
    xr.DataArray
        DataArray containing the column data in molec/cm2
    """
    original_units = column_data.attrs.get("units", "no unit attribute")
    if original_units.lower() in ["molec/cm2", "molec cm-2", "molec/cm^2", "molec cm^-2"]:
        return column_data
    if original_units.lower() not in ["mol/m2", "mol m-2", "mol/m^2", "mol m^-2"]:
        raise ValueError(f"Input column data units are not mol/m2, found {original_units}.")
    with xr.set_options(keep_attrs=True):
        if "multiplication_factor_to_convert_to_molecules_percm2" in column_data.attrs:
            column_data_molec_cm2 = (
                column_data
                * column_data.attrs["multiplication_factor_to_convert_to_molecules_percm2"]
            )
            column_data_molec_cm2.attrs.pop("multiplication_factor_to_convert_to_molecules_percm2")
        else:
            column_data_molec_cm2 = column_data * N_A / M2TOCM2
    column_data_molec_cm2.attrs["units"] = "molec/cm2"
    return column_data_molec_cm2


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
        modelobj_day = modelobj.sel(
            time=slice(day - np.timedelta64(1, "D"), day + np.timedelta64(1, "D"))
        )
        target_time = day + np.timedelta64(overpass_ns, "ns")
        localtime = modelobj_day["time"] + utc_offset_nanoseconds
        if localtime.min() > target_time or localtime.max() < target_time:
            warnings.warn(
                f"Target time {target_time} is outside model time range "
                f"({np.datetime_as_string(localtime.min().values)} - "
                f"{np.datetime_as_string(localtime.max().values)}), skipping."
            )
            continue
        interp = _interpolate_time(target_time, localtime, modelobj_day)
        interp = interp.expand_dims("time", axis=0).assign_coords(time=[target_time])
        interp["time_utc"] = target_time - utc_offset_nanoseconds
        interpolated_data.append(interp)
    concat_data = xr.concat(interpolated_data, dim="time")
    return concat_data


def _interpolate_time(target_time, localtime, data):
    """Applies time interpolation to the data.

    Parameters
    ----------
    target_time : np.datetime64
        Target time to interpolate to.
    localtime : xr.DataArray
        Local time of the model data.
    data : xr.DataArray
        Data to interpolate.

    Returns
    -------
    xr.DataArray
        Interpolated data.
    """
    previous_index, next_index = _calculate_previous_and_next_indices(localtime, target_time)
    previous_weight, next_weight = _calculate_time_weights(
        localtime, target_time, previous_index, next_index
    )
    previous_data = data.isel(time=previous_index).drop_vars("time")
    next_data = data.isel(time=next_index).drop_vars("time")
    interp = (previous_weight * previous_data) + (next_weight * next_data)
    return interp


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
    previous_diff = (target_time - previous_time).astype("float64").drop_vars("time")
    next_diff = (next_time - target_time).astype("float64").drop_vars("time")
    previous_weight = 1 - previous_diff / total_diff
    next_weight = 1 - next_diff / total_diff
    return previous_weight, next_weight


# @numba.jit(nopython=True)
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
    nt, nz, ny, nx = target.shape
    interp = np.zeros((nt, nz, ny, nx))
    for t in range(nt):
        for y in range(ny):
            for x in range(nx):
                interp[t, :, y, x] = np.flip(
                    np.interp(
                        np.flip(target[t, :, y, x]),
                        np.flip(orig[t, :, y, x]),
                        np.flip(data[t, :, y, x]),
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
    p_mid_tropomi = obsobj["pres_pa_mid"].values
    p_orig = modobj["pres_pa_mid"].values
    altitude = calc_altitude_from_thickness(modobj["dz_m"])
    dimensions = ("time", "z", "y", "x")
    coords = {
        "longitude": (("y", "x"), modobj["longitude"].values),
        "latitude": (("y", "x"), modobj["latitude"].values),
    }
    if isinstance(variables, str):
        variables = [variables]
    all_needed_vars = variables + ["temperature_k"]
    for var in all_needed_vars:
        interpolated = _interp_vert(p_orig, p_mid_tropomi, modobj[var].values)
        modsatlayers[var] = xr.DataArray(
            data=interpolated, dims=dimensions, coords=coords, attrs=modobj[var].attrs
        )
    altitude = xr.DataArray(
        data=_interp_vert(p_orig, p_mid_tropomi, altitude.values), dims=dimensions, coords=coords
    )
    modsatlayers["pres_pa_mid"] = xr.DataArray(
        data=p_mid_tropomi,
        dims=dimensions,
        coords=coords,
        attrs=modobj["pres_pa_mid"].attrs,
    )
    modsatlayers["dz_m"] = xr.DataArray(
        data=calc_dz_m_from_altitude(altitude),
        dims=dimensions,
        coords=coords,
        attrs=modobj["dz_m"].attrs,
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
        dictionary containing the keys "averaging_kernel" and
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
    if varname is None:
        varname = default_mod_variable_names[sat_type]
        warnings.warn(f"Variable name not provided, assuming {varname}.")

    mod_p_cols = calc_partialcolumn(modobj, varname, unit="mol/m2")
    if sat_type == "tropomi_l2_no2":
        return apply_averaging_kernel_no2(mod_p_cols, obsobj, averaging_kernel_params=ak_params)
    if sat_type == "tropomi_l2_hcho":
        return apply_averaging_kernel_hcho(mod_p_cols, obsobj, averaging_kernel_params=ak_params)

    return mod_p_cols


def apply_averaging_kernel_no2(mod_p_cols, obsobj, averaging_kernel_params=None):
    """Applies the averaging kernel for TROPOMI NO2 and calculates the column

    Parameters
    ----------
    mod_p_cols : xr.DataArray
        DataArray containing the model NO2 partial columns. It has to be
        previously regridded to satellite space.
    obsobj : xr.Dataset
        Dataset containing all the observational data, including the
        variables related to the averaging kernel.
    averaging_kernel_params : dict[str, str]
        dictionary containing the keys "averaging_kernel" and
        "tropospheric_averaging_kernel_calc" plus, optionally,
        "airmass_factor_total" and "airmass_factor_troposphere".

    Returns
    -------
    xr.DataArray
        DataArray containing the model columns after applying the averaging kernel.
    """
    if averaging_kernel_params["tropospheric_averaging_kernel_calc"]:
        ak = (
            obsobj[averaging_kernel_params["airmass_factor_total"]]
            / obsobj[averaging_kernel_params["airmass_factor_troposphere"]]
            * obsobj[averaging_kernel_params["averaging_kernel"]]
        )
    else:
        ak = obsobj[averaging_kernel_params["averaging_kernel"]]
    if "tm5_tropopause_pressure" in obsobj:
        ak = ak.where(obsobj["pres_pa_mid"] >= obsobj["tm5_tropopause_pressure"], other=0)

    column_data_model = xr.dot(ak, mod_p_cols, dim="z") * N_A / M2TOCM2
    column_data_model.attrs = {
        "description": "Tropospheric column after applying averaging kernel",
        "units": "molec/cm2",
    }
    return column_data_model


def apply_averaging_kernel_hcho(mod_p_cols, obsobj, varname=None, averaging_kernel_params=None):
    """Applies the averaging kernel for TROPOMI HCHO and calculates the column.
    It is in the ATBD documentation, instead of the user guide.

    Parameters
    ----------
    mod_p_cols : xr.DataArray
        DataArray containing the model HCHO partial columns. It has to be
        previously regridded to satellite space.
    obsobj : xr.Dataset
        Dataset containing all the observational data, including the
        variables related to the averaging kernel.
    averaging_kernel_params : dict[str, str]
        dictionary containing the keys "averaging_kernel" and
        "tropospheric_averaging_kernel_calc" plus, optionally,
        "airmass_factor_total" and "airmass_factor_troposphere".

    Returns
    -------
    xr.DataArray
        DataArray containing the model columns after applying the averaging kernel.
    """
    ak = obsobj[averaging_kernel_params["averaging_kernel"]]
    if "tm5_tropopause_pressure" in obsobj:
        ak = ak.where(obsobj["pres_pa_mid"] >= obsobj["tm5_tropopause_pressure"], other=0)

    column_data_model = xr.dot(ak, mod_p_cols, dim="z") * N_A / M2TOCM2
    column_data_model.attrs = {
        "description": f"Tropospheric column of model {varname} after applying averaging kernel",
        "units": "molec/cm2",
    }
    return column_data_model


def apply_averaging_kernel_hcho_or_co(
    mod_p_cols, obsobj, varname=None, averaging_kernel_params=None
):
    """Applies the averaging kernel for TROPOMI HCHO and calculates the column.
    It is in the ATBD documentation, instead of the user guide.

    Parameters
    ----------
    mod_p_cols : xr.DataArray
        DataArray containing the model HCHO partial columns. It has to be
        previously regridded to satellite space.
    obsobj : xr.Dataset
        Dataset containing all the observational data, including the
        variables related to the averaging kernel.
    averaging_kernel_params : dict[str, str]
        dictionary containing the keys "averaging_kernel" and
        "tropospheric_averaging_kernel_calc" plus, optionally,
        "airmass_factor_total" and "airmass_factor_troposphere".

    Returns
    -------
    xr.DataArray
        DataArray containing the model columns after applying the averaging kernel.
    """
    ak = obsobj[averaging_kernel_params["averaging_kernel"]]
    if "tm5_tropopause_pressure" in obsobj:
        ak = ak.where(obsobj["pres_pa_mid"] >= obsobj["tm5_tropopause_pressure"], other=0)

    column_data_model = xr.dot(ak, mod_p_cols, dim="z") * N_A / M2TOCM2
    column_data_model.attrs = {
        "description": f"Tropospheric column of model {varname} after applying averaging kernel",
        "units": "molec/cm2",
    }
    return column_data_model


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
    raise ValueError("time variable (e.g., time_granule) has wrong dimension number.")


def crop_obsobj(obsobj, modobj):
    """Subselects the observations depending on the model domain.

    Parameters
    ----------
    obsobj : xr.Dataset
        Dataset containing all the observational data.
    modobj : xr.Dataset
        Model dataset, as read in by MELODIES-MONET.

    Returns
    -------
    xr.Dataset | None
        Dataset with the observations outside the bounds discarded.
        If no observations are within the model domain, None is returned.
    """
    lonlat_mask = (
        (obsobj.longitude >= modobj.longitude.min())
        & (obsobj.longitude <= modobj.longitude.max())
        & (obsobj.latitude >= modobj.latitude.min())
        & (obsobj.latitude <= modobj.latitude.max())
    )

    valid_x_indices = lonlat_mask.any(dim="y").values.nonzero()[0]
    if valid_x_indices.size == 0:
        return None
    x_min, x_max = valid_x_indices.min(), valid_x_indices.max()

    valid_y_indices = lonlat_mask.any(dim="x").values.nonzero()[0]
    if valid_y_indices.size == 0:
        return None
    y_min, y_max = valid_y_indices.min(), valid_y_indices.max()

    cropped_obsobj = obsobj.isel(x=slice(x_min, x_max), y=slice(y_min, y_max))
    return cropped_obsobj


def select_swaths_overlapping_model(obsobj, modobj):
    """Selects the swaths overlapping with the model domain.

    Parameters
    ----------
    obsobj : dict[str, xr.Dataset]
        Dictionary containing observations
    modobj : xr.Dataset
        Model dataset, as read in by MELODIES-MONET.

    Returns
    -------
    dict[np.datetime64, xr.Dataset]
        Dictionary containing the same keys as the obs, and the
        observations within the model domain as values.
    """
    output_pair = {}
    bounds = [
        modobj["longitude"].min().values,
        modobj["longitude"].max().values,
        modobj["latitude"].min().values,
        modobj["latitude"].max().values,
    ]
    for k in obsobj.keys():
        if within_model_domain(obsobj[k], bounds):
            output_pair[k] = crop_obsobj(obsobj[k], bounds)
        else:
            warnings.warn(f"Swath {k} is outside model domain, skipping.")
    return output_pair


def _regrid_and_apply_ak(
    modobj,
    obsobj,
    mod_var="NO2",
    sat_var="nitrogendioxide_tropospheric_column",
    sat_type="tropomi_l2_no2",
):
    """Regrids and applies AK to one swath.

    Parameters
    ----------
    modobj : xr.Dataset
        model dataset, as read in by MELODIES-MONET. Model should be
        already at overpass time
    obsobj : dict[np.datetime64, xr.Dataset]
        Dictionary containing observations
    mod_var : str
        Variable name in the model dataset
    sat_var : str
        Variable name in the satellite dataset

    Returns
    -------
    dict[str, xr.Dataset]
        Dictionary containing the same keys as the obs, and the model
        data in satellite space after applying the ak as values
    """

    print("obsobj", obsobj)
    output_pair = {}
    modobj_at_overpass_time = interpolate_time(modobj, overpass_time=13.5)
    modobj_at_overpass_time["altitude"] = calc_altitude_from_thickness(
        modobj_at_overpass_time["dz_m"]
    )

    obsobj_dates = np.unique(obsobj["time_granule"].dt.floor("D"))
    modobj_dates_granules = modobj_at_overpass_time["time_utc"].dt.floor("D")
    for d in obsobj_dates:
        if d not in modobj_dates_granules:
            warnings.warn(f"Model does not have data for {d}, skipping.")
            continue
        obsobj_cropped = crop_obsobj(obsobj, modobj)
        if obsobj_cropped is None:
            warnings.warn(f"Swath on {d} is outside model domain, skipping.")
            continue
        # NOTE: We still need to check if this works accross the dateline
        modobj_at_date = modobj_at_overpass_time.where(
            modobj_dates_granules == d, drop=True
        ).drop_vars("time_utc")
        modobj_regrid = interp_horizontal_mod2sat(obsobj_cropped, modobj_at_date)
        modobj_regrid = interp_vertical_mod2swath(obsobj_cropped, modobj_regrid, mod_var)
        # Apply averaging kernel
        modobj_regrid[mod_var] = apply_averaging_kernel(
            modobj_regrid, obsobj_cropped, sat_type, varname=mod_var
        )
        starttime_swath = np.datetime_as_string(obsobj["time_granule"].min().values)
        output_dataset = xr.Dataset()
        output_dataset[mod_var] = modobj_regrid[mod_var]
        output_dataset[sat_var] = tropomi_mol_m2_to_molec_cm2(obsobj_cropped[sat_var])
        output_pair[starttime_swath] = output_dataset
    return output_pair


def regrid_and_apply_ak(
    obsobj_dict,
    modobj,
    mod_var="NO2",
    sat_var="nitrogendioxide_tropospheric_column",
    sat_type="tropomi_l2_no2",
):
    """Regrids and applies AK to multiple swaths.

    Parameters
    ----------
    obsobj_dict : dict[np.datetime64, xr.Dataset]
        Dictionary containing observations
    modobj : xr.Dataset
        model dataset, as read in by MELODIES-MONET. Model should be
        already at overpass time

    Returns
    -------
    dict[str, xr.Dataset]
        Dictionary containing the same keys as the obs, and the model
        data in satellite space after applying the ak as values
    """

    output_pair = {}
    for k in obsobj_dict.keys():
        regridded_swath = _regrid_and_apply_ak(
            modobj, obsobj_dict[k], mod_var=mod_var, sat_var=sat_var, sat_type=sat_type
        )
        output_pair.update(regridded_swath)
    if len(output_pair) == 0:
        raise ValueError("Output pair is empty. Are you sure that there is matching data?")
    return output_pair


def back_to_structured_grid(paired_object, target_grid):
    """Reformats the output dictionary to a structured grid.

    Parameters
    ----------
    paired_object : dict[str, xr.Dataset]
        Dictionary containing the same keys as the obs, and the model
        data in satellite space after applying the ak as values
    target_grid : xr.Dataset
        Dataset containing the target grid information.

    Returns
    -------
    xr.Dataset
        Dataset containing the paired data in the structured grid.
    """

    output_all = []

    for k in paired_object.keys():
        regridder = xe.Regridder(
            paired_object[k],
            target_grid,
            ignore_degenerate=True,
            unmapped_to_nan=True,
            method="bilinear",
        )
        regridded_pair = regridder(paired_object[k])
        output_all.append(regridded_pair)

    output_pair = xr.concat(output_all, dim="time")
    if len(output_pair.time) == 1:
        return output_pair
    return output_pair.groupby("time").mean()
