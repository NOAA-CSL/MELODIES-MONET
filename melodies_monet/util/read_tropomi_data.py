""" Read TROPOMI data into MELODIES-MONET
"""

import glob
import warnings

import netCDF4 as nc4
import numpy as np
import xarray as xr

MILISECONDS_TO_SECONDS = 0.001


def _open_one_dataset(fname, variable_dict):
    """Opens only one dataset
    Based on MONETIO

    Parameters
    ----------
    fname : str
        Input file name
    variable_dict : dict[str, dict]
        Dictionary containing dictionaries with info for each variable
    """
    assert isinstance(variable_dict, dict)
    print(f"reading {fname}")

    ds = xr.Dataset()
    dso = nc4.Dataset(fname, "r")
    ds.attrs = dso.__dict__
    time = open_var_no_format("time", dso)  # base unit in seconds since...
    dtime = open_var_no_format("delta_time", dso)  # in ms
    time_cast = xr.DataArray(data=time[:], dims=time.dimensions, attrs=time.__dict__)
    ds["time"] = xr.conventions.decode_cf_variable("time", time_cast.variable)
    lon = open_var_no_format("longitude", dso)
    lat = open_var_no_format("latitude", dso)
    ds["pres_pa_mid"], ds["pres_pa_int"] = _calc_pressure_levels(dso)
    ds = ds.assign_coords(
        {
            "time": (("time",), ds["time"].values),
            "latitude": (("y", "x"), lat[:].squeeze()),
            "longitude": (("y", "x"), lon[:].squeeze()),
        }
    )
    _set_latlon(ds, lat[:], lon[:])
    ds["time_granule"] = _add_time_granule(time, dtime)

    ds = ds.assign_coords()
    for variable in variable_dict:
        if variable not in ["pres_pa_mid", "tm5_tropopause_pressure"]:
            ds[variable] = _add_variable(variable, dso)
        elif variable == "tm5_tropopause_pressure":
            ds["tm5_tropopause_pressure"] = _calc_tm5_tropopause_pressure(ds, dso)
        if "quality_flag" in variable_dict[variable]:
            ds[variable].attrs["quality_flag"] = variable_dict[variable]["quality_flag"]
            if "qa_thresh_min" in variable_dict[variable]:
                ds[variable].attrs["qa_thresh_min"] = variable_dict[variable]["qa_thresh_min"]
            if "qa_thresh_max" in variable_dict[variable]:
                ds[variable].attrs["qa_thresh_max"] = variable_dict[variable]["qa_thresh_max"]
            ds[variable] = apply_quality_flag(ds[variable], dso)
    dimensions = []
    for x in ["time", "z", "y", "x"]:
        if x in ds.dims:
            dimensions.append(x)
    dso.close()
    ds = ensure_increasing_altitude(ds)
    return ds.transpose(*dimensions, ...)


def _set_latlon(ds, lat, lon):
    """Sets latitude y longitude on TROPOMI ds inplace

    Parameters
    ----------
    ds : xr.Dataset
        Dataset to add latitude and longitude
    lat : np.array
        latitude (range: -90, 90)
    lon : np.array
        longitude (range: -180, 180)

    Returns
    -------
    None
    """
    ds["latitude"] = xr.DataArray(
        data=lat.squeeze(),
        dims=("y", "x"),
        attrs={"units": "degrees_north"},
        coords={"latitude": (("y", "x"), lat[:].squeeze())},
    )
    ds["longitude"] = xr.DataArray(
        data=lon.squeeze(),
        dims=("y", "x"),
        attrs={"units": "degrees_east"},
        coords={"longitude": (("y", "x"), lon[:].squeeze())},
    )


def ensure_increasing_altitude(ds):
    """Ensures that the altitude is increasing (i.e, the pressure should
    decrease as z increases)

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with the satellite data. If pressure is not included,
        nothing will be done.

    Returns
    -------
    xr.Dataset
        Dataset with corrected pressure
    """
    if ("pres_pa_mid" not in ds) and ("pres_pa_int" not in ds):
        warnings.warn("Missing pressure information. Ignoring vertical directionality check")
        return ds
    vertical_dim = {"pres_pa_mid": "z", "pres_pa_int": "z_stagg"}
    for pres_var, vert_dim in vertical_dim.items():
        if (ds[pres_var].isel(time=0).isel(**{vert_dim: slice(0, 10)}).diff(dim="z") > 0).any():
            ds = ds.sel(**{vert_dim: slice(None, None, -1)})
    return ds


def _add_time_granule(time, dtime):
    """Sets time of granule/pixel from TROPOMI

    Parameters
    ----------
    time : nc4.Variable
        reference time in seconds since reference
    dtime : nc4.variable
        Dataset with original data

    Returns
    -------
    xr.DataArray
        DataArray containing the time of each granule
    """
    _time_granule = time[:] + dtime[:] * MILISECONDS_TO_SECONDS
    if len(_time_granule.shape) == 2:
        time_granule = xr.DataArray(data=_time_granule, dims=("time", "y"), attrs=time.__dict__)
    elif len(_time_granule.shape) == 3:
        time_granule = xr.DataArray(
            data=_time_granule, dims=("time", "y", "x"), attrs=time.__dict__
        )
    else:
        raise ValueError("Could not assign the time of each granule. Check data dimensions.")

    time_granule = xr.conventions.decode_cf_variable("time_granule", time_granule.variable)
    return time_granule


def _walktree_search(variable, netcdf_dataset, path=""):
    """Recursive search for each variable

    Parameters
    ----------
    variable : str
        Variable name
    netcdf_dataset : nc4.Dataset
        netCDF4 dataset with data to search
    path : str, optional
        Path within the netCDF4 dataset, by default ""

    Returns
    -------
    str
        str path that was searched for.

    Raises
    ------
    ValueError
        If variable is not found
    """
    if variable in netcdf_dataset.variables:
        return f"{path}/{variable}"
    groups = netcdf_dataset.groups
    for group in groups:
        try:
            return _walktree_search(variable, netcdf_dataset[group], f"{path}/{group}")
        except ValueError:
            continue
    raise ValueError(f"{variable} not found in {netcdf_dataset}")


def _add_variable(variable, netcdf_dataset):
    """Creates xr.DataArray formatted for MELODIES-MONET.

    Parameters
    ----------
    variable : str
        Variable name
    netcdf_dataset : nc4.Dataset
        nc4.Dataset with data to search

    Returns
    -------
    xr.DataArray
        DataArray with the variable that was searched for, formatted
        for MELODIES-MONET
    """
    var = open_var_no_format(variable, netcdf_dataset)
    _replacements = {"layer": "z", "scanline": "y", "ground_pixel": "x"}
    _dimensions = list(var.dimensions)
    dimensions = [_replacements[x] if x in _replacements else x for x in _dimensions]
    dtype = var[:].dtype
    if np.issubdtype(dtype, np.integer):
        var_values = var[:].filled(np.iinfo(dtype).min)
        da = xr.DataArray(data=var_values, dims=dimensions, attrs=var.__dict__).astype(dtype)
    else:
        da = xr.DataArray(data=var[:], dims=dimensions, attrs=var.__dict__).astype(dtype)
    return da


def open_var_no_format(variable, netcdf_dataset):
    """Opens only one variable from a netCDF4 dataset

    Parameters
    ----------
    variable : str
        Variable name
    netcdf_dataset : nc4.Dataset
        nc4.Dataset with data to search

    Returns
    -------
    xr.DataArray
        DataArray with the variable that was searched for"""
    return netcdf_dataset[_walktree_search(variable, netcdf_dataset)]


def _calc_pressure_levels(netcdf_tropomi, product="check"):
    """Calculates pressure levels.

    Parameters
    ----------
    netcdf_tropomi : nc4.Dataset
        Dataset from TROPOMI files
    product : str  ("no2" "hcho")
        str indicating which product it is. If 'check', all are tried

    Returns
    -------
    (xr.DataArray, xr.DataArray)
        Two DataArrays containing the midlevel pressure and the
        pressure at the layer interface respectively.
    """
    if ("id" in netcdf_tropomi.ncattrs()) and ("_CO_" in netcdf_tropomi.id):
        pressure_level_bottom = _add_variable("pressure_levels", netcdf_tropomi)
        return _calc_pressure_tropomi_co(pressure_level_bottom)
    tm5_constant_a = _add_variable("tm5_constant_a", netcdf_tropomi)
    tm5_constant_b = _add_variable("tm5_constant_b", netcdf_tropomi)
    surface_pressure = _add_variable("surface_pressure", netcdf_tropomi)

    dims = tm5_constant_a.dims
    if "vertices" in dims or product == "no2":
        return _calc_pressure_tropomi_no2(tm5_constant_a, tm5_constant_b, surface_pressure)
    if "time" in dims or product == "hcho":
        return _calc_pressure_tropomi_hcho(tm5_constant_a, tm5_constant_b, surface_pressure)
    raise ValueError(f"Dims in tm5_constant_a {dims=} do not match expectations.")


def _calc_pressure_tropomi_no2(tm5_constant_a, tm5_constant_b, surface_pressure):
    """Calculates the pressure for the TROPOMI L2 NO2 product

    Parameters
    ----------
    tm5_constant_a : xr.DataArray
        constant A to calculate pressure.
    tm5_constant_b : xr.DataArray
        constant b to calculate pressure.
    surface_pressure : xr.DataArray
        surface pressure of tropomi.

    Returns
    -------
    xr.DataArray, xr.DataArray
        Midlayer pressure and interface pressure
    """
    num_times, num_y, num_x = surface_pressure.shape
    num_layers = len(tm5_constant_a[:, 0])
    midlayer_pressure = xr.DataArray(
        data=np.zeros((num_times, num_layers, num_y, num_x), dtype=np.float64),
        dims=("time", "z", "y", "x"),
    )
    for i in range(num_layers):
        midlayer_pressure[:, i, :, :] = (
            tm5_constant_a[i, 0].values
            + tm5_constant_b[i, 0].values * surface_pressure[:].values
            + tm5_constant_a[i, 1].values
            + tm5_constant_b[i, 1].values * surface_pressure[:].values
        ) / 2
    midlayer_pressure.attrs = {"units": "Pa", "long_name": "midlayer_pressure_in_pa"}
    interface_pressure = xr.DataArray(
        data=np.zeros((num_times, num_layers + 1, num_y, num_x), dtype=np.float64),
        dims=("time", "z_stagg", "y", "x"),
    )
    interface_pressure[:, 0, :, :] = surface_pressure[:]
    for i in range(0, num_layers):
        interface_pressure[:, i + 1, :, :] = (
            tm5_constant_a[i, 1] + tm5_constant_b[i, 1] * surface_pressure[:]
        )
    return midlayer_pressure, interface_pressure


def _calc_pressure_tropomi_hcho(tm5_constant_a, tm5_constant_b, surface_pressure):
    """Calculates the pressure for the TROPOMI L2 NO2 product

    Parameters
    ----------
    tm5_constant_a : xr.DataArray
        constant A to calculate pressure.
    tm5_constant_b : xr.DataArray
        constant b to calculate pressure.
    surface_pressure : xr.DataArray
        surface pressure of tropomi.

    Returns
    -------
    xr.DataArray, xr.DataArray
        Midlayer pressure and interface pressure
    """
    num_times, num_y, num_x = surface_pressure.shape
    num_layers = len(tm5_constant_a.isel(time=0))
    interface_pressure = xr.DataArray(
        data=np.zeros((num_times, num_layers + 1, num_y, num_x), dtype=np.float64),
        dims=("time", "z_stagg", "y", "x"),
    )
    interface_pressure[:, 0, :, :] = surface_pressure[:]
    for i in range(0, num_layers):
        interface_pressure[:, i + 1, :, :] = (
            tm5_constant_a[0, i].values + tm5_constant_b[0, i].values * surface_pressure[:]
        )
    midlayer_pressure = xr.DataArray(
        data=np.zeros((num_times, num_layers, num_y, num_x), dtype=np.float64),
        dims=("time", "z", "y", "x"),
    )
    for i in range(num_layers):
        midlayer_pressure[:, i, :, :] = (
            interface_pressure[:, i, :, :].values + interface_pressure[:, i + 1, :, :].values
        ) / 2
    midlayer_pressure.attrs = {"units": "Pa", "long_name": "midlayer_pressure_in_pa"}
    return midlayer_pressure, interface_pressure


def _calc_pressure_tropomi_co(pressure_level_bottom):
    """Calculates interface and midlayer pressure for CO.

    Parameters
    ----------
    pressure_level_bottom : xr.DataArray
        DataArray containing all the pressure at mid layer

    Returns
    -------
    xr.DataArray, xr.DataArray
        DataArrays containing the pressure at the interface and at midlevel
    """
    pressure_level_bottom_transpose = pressure_level_bottom.transpose("time", "z", "y", "x")
    num_times, num_layers, num_y, num_x = pressure_level_bottom_transpose.shape
    interface_pressure = xr.DataArray(
        data=np.zeros((num_times, num_layers + 1, num_y, num_x), dtype=np.float64),
        dims=("time", "z_stagg", "y", "x"),
        attrs={"long_name": "pressure_interface", "units": "Pa"},
    )
    interface_pressure[:, 1:, :, :] = pressure_level_bottom_transpose.values
    midlayer_pressure = xr.DataArray(
        data=np.zeros((num_times, num_layers, num_y, num_x), dtype=np.float64),
        dims=("time", "z", "y", "x"),
        attrs={"long_name": "pressure_midlayer", "units": "Pa"},
    )
    midlayer_pressure[:, :, :, :] = (
        interface_pressure[:, :-1, :, :].values + interface_pressure[:, 1:, :, :].values
    ) / 2
    return midlayer_pressure, interface_pressure


def _calc_tm5_tropopause_pressure(processed_data, netcdf_tropomi):
    """Calculates the TM5 tropopause pressure

    Parameters
    ----------
    processsed_data : xr.Dataset
        Dataset containing processed data. It has to include
        'pres_pa_mid' as a variable
    netcdf_tropomi : nc4.Dataset
        Dataset containing the netCDF4 tropomi file

    Returns
    -------
    xr.DataArray
        DataArray with the tropopause pressure.
    """
    tm5_tropopause_pressure_idx = _add_variable("tm5_tropopause_layer_index", netcdf_tropomi)
    tm5_tropopause_pressure_idx = tm5_tropopause_pressure_idx.where(
        (tm5_tropopause_pressure_idx > 0) & (tm5_tropopause_pressure_idx < 10000), other=-1
    )
    tropopause_pressure = processed_data["pres_pa_mid"].isel(z=tm5_tropopause_pressure_idx)
    return tropopause_pressure


def apply_quality_flag(variable, netcdf_tropomi):
    """Applies quality_flags inplace

    Parameters
    ----------
    variable : str
        Variable containing the attribute qa_thersh_min.
    netcdf_tropomi : nc4.Dataset
        Dataset containing the netCDF4 tropomi file

    Returns
    -------
    xr.DataArray
        DataArray with applied quality flag
    """
    assert "quality_flag" in variable.attrs, f"quality_flag not in {variable.name}"
    assert ("qa_thresh_min" in variable.attrs) or (
        "qa_thresh_max" in variable.attrs
    ), f"Neither qa_thresh_min nor qa_thresh_max in {variable.name}"
    qa = _add_variable(variable.attrs["quality_flag"], netcdf_tropomi)
    if "qa_thresh_min" in variable.attrs:
        variable = variable.where(qa >= variable.attrs["qa_thresh_min"])
    if "qa_thresh_max" in variable.attrs:
        variable = variable.where(qa <= variable.attrs["qa_thresh_max"])
    return variable


def open_datasets(all_files, variable_dict):
    """Creates a dict containing all the datasets

    Parameters
    ----------
    all_files : str, list[str]
        String or list of strings containing all the files that should
        be opened. Wildcards are supported.
    variable_dict : dict[str, dict]
        Dictionary of dictionaries for each variable.

    Returns
    -------
    dict[str, xr.Dataset]
        Dictionary with time reference as keys and xr.Dataset containing
        satellite information as values.
    """

    if isinstance(all_files, str):
        datasets = sorted(glob.glob(all_files))
    elif isinstance(all_files, list):
        datasets = []
        for ds in all_files:
            datasets.extend(glob.glob(ds))
        datasets = sorted(datasets)
    ds_collection = {}
    for data in datasets:
        d = _open_one_dataset(data, variable_dict)
        # Select time coverage start as key, removing the trailing Z
        ds_collection[d.attrs["time_coverage_start"].replace("Z", "")] = d
    return ds_collection
