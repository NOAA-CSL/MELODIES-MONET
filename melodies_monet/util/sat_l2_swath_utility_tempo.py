# Copyright (C) 2024 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#

# read all swath data for the time range
# developed for TEMPO Level2 NO2
#

import collections
import glob
import logging
import warnings

import numba
import numpy as np
import xarray as xr
import xesmf as xe

numba_logger = logging.getLogger("numba")
numba_logger.setLevel(logging.WARNING)


def tempo_interp_mod2swath(obsobj, modobj, method="bilinear"):
    """Interpolate model to satellite swath/swaths

    Parameters
    ----------
    obsobj : xr.Dataset
        satellite with swath data.
    modobj : xr.Dataset
        model data (with no2 col calculated)

    Returns
    -------
    xr.Dataset
        Regridded model data at swath or swaths. If type is xr.Dataset, a single
        swath is returned. If type is collections.OrderedDict, it returns an
        OrderedDict in which each time represents the reference time of the swath.
    """

    mod_at_swathtime = modobj.interp(time=obsobj.time.mean())
    regridder = xe.Regridder(mod_at_swathtime, obsobj, method, unmapped_to_nan=True)
    modswath = regridder(mod_at_swathtime)
    return modswath


def _calc_dp(obsobj):
    """Calculate delta pressure in satellite layers

    Parameters
    ----------
    obsobj : xr.Dataset
        satellite observations containing pressure (in Pa)

    Returns
    -------
    xr.DataArray
        Pressure difference in layer
    """

    # REMINDER: pressure is higher at lower vertical levels: dp is positive
    # only if defined as lower - higher.
    dp_vals = (
        obsobj["pressure"].isel(swt_level_stagg=slice(None, -1)).values
        - obsobj["pressure"].isel(swt_level_stagg=slice(1, None)).values
    )
    dp = xr.DataArray(
        data=dp_vals,
        dims=("swt_level", "x", "y"),
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
    nz, nx, ny = target.shape
    interp = np.zeros((nz, nx, ny))
    # for t in nt:
    for x in range(nx):
        for y in range(ny):
            interp[:, x, y] = np.flip(
                np.interp(
                    np.flip(target[:, x, y]),
                    np.flip(orig[:, x, y]),
                    np.flip(data[:, x, y]),
                )
            )
    return interp


def interp_vertical_mod2swath(obsobj, modobj, vars=["NO2_col"]):
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
    xr.Dataset
        Model data (interpolated to TEMPO vertical layers
    """
    assert np.all(modobj["lon"].fillna(0).values == obsobj["lon"].fillna(0).values)
    assert np.all(modobj["lat"].fillna(0).values == obsobj["lat"].fillna(0).values)

    modsatlayers = xr.Dataset()
    p_mid_tempo = (
        obsobj["pressure"].isel(swt_level_stagg=slice(None, -1)).values
        + obsobj["pressure"].isel(swt_level_stagg=slice(1, None)).values
    ) / 2
    p_orig = modobj["pres_pa_mid"].values
    # dimensions = ("time", "z", "lon", "lat")
    dimensions = ("z", "x", "y")
    coords = {
        # "time": (("time",), modobj["time"].values),
        "lon": (("x", "y"), modobj["lon"].values),
        "lat": (("x", "y"), modobj["lat"].values),
    }
    for var in vars:
        interpolated = _interp_vert(p_orig, p_mid_tempo, modobj[var].values)
        modsatlayers[var] = xr.DataArray(
            data=interpolated, dims=dimensions, coords=coords, attrs=modobj[var].attrs
        )
    modsatlayers["p_mid_tempo"] = xr.DataArray(
        data=p_mid_tempo,
        dims=dimensions,
        coords=coords,
        attrs=modobj["pres_pa_mid"].attrs,
    )
    _interp_description = "Mid layer pressure interpolated to TEMPO mid swt_layer pressures"
    modsatlayers["p_mid_tempo"].attrs["description"] = _interp_description
    return modsatlayers


def calc_partialcolumn(modobj, var="NO2"):
    """Calculates the partial column of a species from its concentration.

    Parameters
    ----------
    modobj : xr.Dataset
        Model data

    Returns
    -------
    xr.DataArray
        DataArray containing the partial column of the species.
    """
    R = 8.314  # m3 * Pa / K / mol
    NA = 6.022e23
    PPBTOMOLMOL = 1e-9
    M2TOCM2 = 1e4
    fac_units = PPBTOMOLMOL * NA / M2TOCM2
    layer_thickness = _calc_layer_thickness(modobj)
    partial_col = (
        modobj["NO2"]
        * modobj["pres_pa_mid"]
        * layer_thickness
        * fac_units
        / (R * modobj["temperature_k"])
    )
    return partial_col


def _calc_layer_thickness(modobj):
    """Calculates layer thickness

    Parameters
    ----------
    modbj : xr.Dataset
        Model data as produced by the MONETIO reader.

    Returns
    -------
    xr.DataArray
        Layer thickness in m.
    """
    height_agl = modobj["layer_height_agl"]
    layer_thickness = xr.zeros_like(height_agl)
    layer_thickness[0, :, :] = height_agl.isel(z=0).values
    layer_thickness[1:, :, :] = (
        height_agl.isel(z=slice(1, None)).values - height_agl.isel(z=slice(0, -1)).values
    )
    return layer_thickness


def apply_weights_mod2tempo_no2_hydrostatic(obsobj, modobj):
    """Apply the scattering weights and air mass factors accordint to
    Cooper et. al, 2020, doi: https://doi.org/10.5194/acp-20-7231-2020,
    assuming the hydrostatic equation. It does not require temperature
    nor geometric layer_thickness.

    Parameters
    ----------
    obsobj : xr.Dataset
        TEMPO data, including pressure and scattering weights
    modobj : xr.Dataset
        Model data, already interpolated to TEMPO grid

    Returns
    -------
    xr.DataArray
        A xr.DataArray containing the NO2 model data after applying
        the air mass factors and scattering weights
    """
    unit_c = 6.022e23 * 9.8 / 1e4
    dp = _calc_dp(obsobj).rename({"swt_level": "z"})
    PPBTOMOLMOL = 1e-9
    # import pdb; pdb.set_trace()
    tropopause_pressure = obsobj["tropopause_pressure"]
    scattering_weights = obsobj["scattering_weights"].transpose("swt_level", "x", "y")
    scattering_weights = scattering_weights.rename({"swt_level": "z"})
    scattering_weights = scattering_weights.where(modobj["p_mid_tempo"] >= tropopause_pressure)
    modno2 = modobj["NO2"].where(modobj["p_mid_tempo"] >= tropopause_pressure)
    amf_troposphere = obsobj["amf_troposphere"]
    modno2col_trfmd = (dp * scattering_weights * modno2).sum(dim="z") * unit_c * PPBTOMOLMOL
    modno2col_trfmd = modno2col_trfmd.where(modno2.isel(z=0).notnull())
    modno2col_trfmd = modno2col_trfmd / amf_troposphere
    return modno2col_trfmd


def apply_weights_mod2tempo_no2(obsobj, modobj):
    """Apply the scattering weights and air mass factors according to
    Cooper et. al, 2020, doi: https://doi.org/10.5194/acp-20-7231-2020

    Parameters
    ----------
    obsobj : xr.Dataset
        TEMPO data, including pressure and scattering weights
    modobj : xr.Dataset
        Model data, already interpolated to TEMPO grid

    Returns
    -------
    xr.DataArray
        A xr.DataArray containing the NO2 model data after applying
        the air mass factors and scattering weights
    """
    partial_col = modobj["NO2_col"]

    # import pdb; pdb.set_trace()
    tropopause_pressure = obsobj["tropopause_pressure"] * 100
    scattering_weights = obsobj["scattering_weights"].transpose("swt_level", "x", "y")
    scattering_weights = scattering_weights.rename({"swt_level": "z"})
    scattering_weights = scattering_weights.where(modobj["p_mid_tempo"] >= tropopause_pressure)
    amf_troposphere = obsobj["amf_troposphere"]
    modno2col_trfmd = (scattering_weights * partial_col).sum(dim="z") / amf_troposphere
    modno2col_trfmd = modno2col_trfmd.where(modobj["NO2_col"].isel(z=0).notnull())
    modno2col_trfmd.attrs = {
        "units": "molecules/cm2",
        "description": "model tropospheric column after applying TEMPO scattering weights and AMF",
        "history": "Created by MELODIES-MONET, apply_weights_mod2tempo_no2, TEMPO util",
    }
    return modno2col_trfmd


def _regrid_and_apply_weights(obsobj, modobj):
    """Does the complete process of regridding and
    applying scattering weights. Assumes that obsobj is a Dataset

    Parameters
    ----------
    obsobj : xr.Dataset
        TEMPO observations
    modobj : xr.Dataset
        Model data

    Returns
    -------
    xr.DataArray
        Model data regridded to the TEMPO grid,
        with the averaging kernel.
    """
    modobj_hs = tempo_interp_mod2swath(obsobj, modobj)
    if "layer_height_agl" in list(modobj.variables):
        modobj_hs["NO2_col"] = calc_partialcolumn(modobj_hs)
        modobj_swath = interp_vertical_mod2swath(obsobj, modobj_hs, ["NO2_col"])
        da_out = apply_weights_mod2tempo_no2(obsobj, modobj_swath)
    else:
        warnings.warn(
            "There is no layer_height_agl variable, and the partial column"
            + "cannot be directly calculated. Assuming hydrostatic equation."
        )
        modobj_swath = interp_vertical_mod2swath(obsobj, modobj_hs, ["NO2"])
        da_out = apply_weights_mod2tempo_no2_hydrostatic(obsobj, modobj_swath)
    return da_out


def regrid_and_apply_weights(obsobj, modobj, pair=True):
    """Does the complete process of regridding
    and applying scattering weights.

    Parameters
    ----------
    obsobj : xr.Dataset | collections.OrderedDict
    modobj : xr.Dataset

    Returns
    -------
    xr.Dataset | collections.OrderedDict
        Model with regridded data. If obsobj is of type collections.OrderedDict,
        an OrderedDict is returned.
    """

    if isinstance(obsobj, xr.Dataset):
        regridded = _regrid_and_apply_weights(obsobj, modobj)
        output = regridded.to_dataset(name="NO2_col_wsct")
        output.attrs["reference_time_string"] = obsobj.attrs["reference_time_string"]
        if pair:
            output = xr.merge([output, obsobj["vertical_column_troposphere"]])
        if "lat" in output.variables:
            output = output.rename({"lat": "latitude", "lon": "longitude"})
        return output
    elif isinstance(obsobj, collections.OrderedDict):
        output_multiple = collections.OrderedDict()
        for ref_time in obsobj.keys():
            output_multiple[ref_time] = _regrid_and_apply_weights(
                obsobj[ref_time],
                modobj,
            ).to_dataset(name="NO2_col_wsct")
            output_multiple[ref_time].attrs["reference_time_string"] = ref_time
            output_multiple[ref_time].attrs["scan_num"] = obsobj[ref_time].attrs["scan_num"]
            output_multiple[ref_time].attrs["granule_number"] = obsobj[ref_time].attrs[
                "granule_number"
            ]
            if pair:
                output_multiple[ref_time] = xr.merge(
                    [output_multiple[ref_time], obsobj[ref_time]["vertical_column_troposphere"]]
                )
            if "lat" in output_multiple[ref_time].variables:
                output_multiple[ref_time] = output_multiple[ref_time].rename(
                    {"lat": "latitude", "lon": "longitude"}
                )
        return output_multiple
    else:
        raise Exception("Obsobj must be xr.Dataset or collections.OrderedDict")


def back_to_modgrid(
    paireddict,
    modobj,
    keys_to_merge="all",
    add_time=True,
    to_netcdf=False,
    out_name="Regridded_object_XYZ.nc",
):
    """Grids object in sat-space to modgrid. Designed to grid back to modgrid after applying
    the scattering weights and air mass factors. It is designed for a single scan.

    Parameters
    ----------
    paireddict : collections.OrderedDict[str, xr.Dataset]
        An OrderedDict with time_reference strings as keys.
    modobj : xr.Dataset
        A modobj including the modgrid.
    keys_to_merge : str | list[str]
        If 'all', all keys are assumed to be part of the same scan and merged.
        Else, only the keys provided are merged.
    add_time : bool
        If True, add reference time as a coordinate for the scan.
        Can be useful to concatenate later if multiple scans are required.
    to_netcdf : bool
        If True, save a netcdf with the paired data
    out_name : str
        The base name to save the files if to_netcdf is True. XX will be replaced
        with the scan number and the reference time. If to_netcdf is False, this will
        be ignored.

    Returns
    -------
    xr.Dataset
        Dataset with obj2grid regridded to modobj.
    """
    if keys_to_merge == "all":
        ordered_keys = sorted(list(paireddict.keys()))
    else:
        ordered_keys = sorted(list(keys_to_merge))
    concatenated = paireddict[ordered_keys[0]]
    scan_num = concatenated.attrs["scan_num"]
    granules = [concatenated.attrs["granule_number"]]
    ref_times = [concatenated.attrs["reference_time_string"][:-1]]  # Remove unneded Z
    if len(paireddict) > 1:
        for k in ordered_keys[1:]:
            ds_to_add = paireddict[k]
            if ds_to_add.attrs["scan_num"] != scan_num:
                raise Exception(
                    "back_to_modgrid is prepared to work with data of a single scan. "
                    + f"However, {ordered_keys[0]} is from scan {scan_num} and "
                    + f"{k} if from scan {ds_to_add.attrs['scan_num']}."
                )
            concatenated = xr.concat([concatenated, paireddict[k]], dim="x")
            granules.append(paireddict[k].attrs["granule_number"])
            ref_times.append(paireddict[k].attrs["reference_time_string"][:-1])
    regridder = xe.Regridder(concatenated, modobj, method="bilinear", unmapped_to_nan=True)
    out_regridded = regridder(concatenated)
    # out_regridded = out_regridded.rename({"longitude": "lon", "latitude": "lat"})
    for v in out_regridded.variables:
        if v in concatenated.variables:
            out_regridded[v].attrs = concatenated[v].attrs
        else:
            warnings.warn(f"Variable {v} not found in mod2grid nor obs2grid. Continuing.")
    out_regridded.attrs["reference_time_string"] = ref_times
    out_regridded.attrs["granules"] = np.array(granules)
    scan_num = concatenated.attrs["scan_num"]
    out_regridded.attrs["scan_num"] = scan_num
    if add_time:
        time = [np.array(ref_times[0], dtype="datetime64[ns]")]
        da_time = xr.DataArray(
            name="time",
            data=time,
            dims=["time"],
            attrs={"description": "Reference start time of first selected granule in scan."},
            coords={"time": (("time",), time)},
        )
        out_regridded = out_regridded.expand_dims(time=da_time)
        out_regridded["end_time"] = (
            ("time",),
            [np.array(ref_times[-1], dtype="datetime64[ns]")],
        )
        out_regridded["end_time"].attrs = {
            "description": "time at which the last swath of the scan starts"
        }
    if to_netcdf:
        if "XYZ" in out_name:
            scan_num = out_regridded.attrs["scan_num"]
            out_regridded.to_netcdf(
                out_name.replace(
                    "XYZ",
                    f"S{scan_num:03d}_{out_regridded['time'].values.astype(str)[0][0:19]}",
                )
            )
        else:
            out_regridded.to_netcdf(out_name)

    return out_regridded


def back_to_modgrid_multiscan(
    paireddict, modobj, to_netcdf=False, out_name="Regridded_object_XYZ.nc"
):
    """Grids object in sat-space to modgrid. Designed to grid back to modgrid after applying
    the scattering weights and air mass factors. It is designed for multiple scans, and uses
    back_to_modgrid under the hood. Generally, back_to_modgrid should only be used if
    you can ensure that you are reading only one scan at a time.

    Parameters
    ----------
    paireddict : collections.OrderedDict[str, xr.Dataset]
        An OrderedDict with time_reference strings as keys.
    modobj : xr.Dataset
        A modobj including the modgrid.
    subset_var : str | list[str]
        A list containing a subset of variables to regrid.
        If it is 'all', all variables are regridded
    to_netcdf : bool
        If True, save a netcdf with the paired data
    out_name : str
        The base name to save the files if to_netcdf is True. XX will be replaced
        with the first and list times. If to_netcdf is False, this will
        be ignored.

    Returns
    -------
    xr.Dataset
        Dataset with obj2grid regridded to modobj.
    """
    out_regridded = xr.Dataset()
    ordered_keys = sorted(list(paireddict.keys()))
    scan_num = paireddict[ordered_keys[0]].attrs["scan_num"]
    keys_in_scan = [ordered_keys[0]]
    if len(ordered_keys > 1):
        for k in ordered_keys[1:]:
            if paireddict[k].attrs["scan_num"] == scan_num:
                keys_in_scan.append(k)
            else:
                regridded_scan = back_to_modgrid(paireddict, modobj, keys_in_scan)
                out_regridded = xr.merge([out_regridded, regridded_scan])
                scan_num = paireddict[ordered_keys[0]].attrs["scan_num"]
                keys_in_scan = [ordered_keys[0]]
        regridded_scan = back_to_modgrid(paireddict, modobj, keys_in_scan)
        out_regridded = xr.merge([out_regridded, regridded_scan])

    if to_netcdf:
        if "XYZ" in out_name:
            first_time = out_regridded["time"][0].values.astype(str)[0:19]
            last_time = out_regridded["time"][-1].values.astype(str)[0:19]
            scan_num = out_regridded.attrs["scan_num"]
            out_regridded.to_netcdf(out_name.replace("XYZ", f"{first_time}_{last_time}"))
        else:
            out_regridded.to_netcdf(out_name)

    return out_regridded


def _subset_ds(ds, subset_vars="all"):
    """Subset of variables to select from an xr.Dataset. If subset_vars is all, it just returns
    the Dataset whithout doing anything

    Parameters
    ----------
    ds : xr.Dataset
        Dataset to subset
    subset_vars : str | list[str]
        String or list of strings with the names of the variables to subset.

    Returns
    -------
    xr.Dataset
        Dataset with the subset of variables of interest
    """
    if subset_vars == "all":
        return ds
    else:
        if isinstance(subset_vars, str):
            subset_vars = [subset_vars]
        return ds[subset_vars]


def read_paired_gridded_tempo_model(path):
    """Reads in paired gridded tempo and model data

    Parameters
    ----------
    path : str or globobject

    Returns
    -------
    combined dataset with paired tempo and gridded data.
    """

    pathlist = sorted(glob.glob(path))
    return xr.open_mfdataset(pathlist)


def save_swath(moddict, path="Paired_swath_XYZ.nc"):
    """Saves each swath individually

    Parameters
    ----------
    moddict : collections.OrderedDict[str, xr.Dataset]
        Ordered dict containing all of the paired model and swath.
    path : str | list[str]
        Path to save the swath. If XYZ is present, it will be replaced
        by the date, number of scan and number of granule. If it is a list,
        it will use the list elements to save the swaths.

    Returns
    -------
    None
    """
    for i, k in enumerate(moddict.keys()):
        scan_num = moddict[k].attrs["scan_num"]
        gran_num = moddict[k].attrs["granule_number"]
        if isinstance(path, str):
            if "XYZ" in path:
                pathout = path.replace("XYZ", f"{k[:-1]}_S{scan_num:03d}G{gran_num:03d}")
            else:
                pathout = "1" + path
        elif isinstance(path, list):
            pathout = list[i]
        moddict[k].to_netcdf(pathout)
