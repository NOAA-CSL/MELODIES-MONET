#!/usr/bin/env python
"""Boulder Air"""

# This is written to concatenate boulder air data
# and add it for usage

import argparse
import glob
import warnings

import pandas as pd
import xarray as xr


def create_file_with_coords(
    path_coords, path_data_root, variables="all", resample_freq=None, method="mean"
):
    """Read coordinates from BoulderAir file.
    File should be a csv have columns site_abbreviation, lat, lon

    Parameters:
    -----------
    path_coords : str
        string containing the path to the coordinates file.
    path_data_root : str
        string containing the path to the directory where the files
        containing the measurements are located. The name of those files
        should contain the abbreviation XYZ or wildcards where the abbreviation
        of the station is located.
    variables : str | list[str]
        list of variables to concatenate
    resample_freq : str
        pandas-like resample frequency. pandas is used under the hood.
    method : str
        If resample_freq exists, it will be used to do the interpolation
        to instantaneous values, min, max, mean or median. Possible values:
        'inst', 'min', 'max', 'mean', 'median'.

    Returns
    -------
    pd.DataFrame
        DataFrame containing all the relevant information
    """

    coordinates_data = pd.read_csv(path_coords, comment="#")
    abbrevs = coordinates_data["site_abbreviation"]

    compile_data = xr.Dataset()
    print(abbrevs)
    for x, abbrev in enumerate(abbrevs):
        if "XYZ" in path_data_root:
            path_data = path_data_root.replace("XYZ", abbrev)
        all_data_paths = sorted(glob.glob(path_data))

        # Empty lists are False, following PEP8
        if not all_data_paths:
            warnings.warn(f"{path_data} does not exist. Skipping.")
            continue
        data = _concat_data(all_data_paths, resample_freq, method)
        time_vals = pd.to_datetime(data.index.values)
        ds = xr.Dataset()
        ds["time"] = xr.DataArray(
            data=time_vals,
            dims=["time"],
            coords={"time": (["time"], time_vals)},
        )
        site_coords = coordinates_data.loc[abbrevs == abbrev]
        ds["siteid"] = xr.DataArray(data=[abbrev], dims=["x"], coords={"x": (["x"], [x])})
        ds["longitude"] = xr.DataArray(
            data=site_coords["lon"].values,
            dims=["x"],
            coords={"x": (["x"], [x])},
            attrs={"units": "degrees_east"},
        )
        ds["latitude"] = xr.DataArray(
            data=site_coords["lat"].values,
            dims=["x"],
            coords={"x": (["x"], [x])},
            attrs={"units": "degrees_north"},
        )
        ds["elevation"] = xr.DataArray(
            data=site_coords["m_asl"], dims=["x"], coords={"x": (["x"], [x])}
        )
        if variables == "all":
            variables = list(data.keys())
        elif isinstance(variables, str):
            variables = [variables]

        for variable in variables:
            if (variable == "nox") and ("nox" not in data):
                if ("no2" in data) and ("no" in data):
                    data["nox"] = data["no"] + data["no2"]
            try:
                vals = data[variable].values
            except KeyError:
                warnings.warn(f"{variable} not found in {path_data}. Skipping")
                continue
            ds[variable] = xr.DataArray(
                data=vals,
                dims=["time", "x"],
                coords={"time": (["time"], time_vals), "x": (["x"], [x])},
                attrs={"units": data[variable].keys()[0]},
            )
            ds[variable] = ds[variable].expand_dims(dim={"y": 1})
            ds[variable] = ds[variable].transpose("time", "y", "x")
        compile_data = xr.merge([compile_data, ds])
    # compile_data = compile_data.transpose("time", "y", "x")
    return compile_data


def _concat_data(all_data_paths, resample_freq=None, method="mean"):
    """Loads the data. If inst_resample is not None, tries to interpolate
    to inst_resample frequency. If mean_resample is not None, tries to
    resample using mean values.

    Parameters
    ----------
    all_data_paths: list[str],
        List containing the paths to the data of all the stations.
    resample_freq: str | None
        Frequency of the resample. Uses pandas frequencies. If is None,
        the data will not be resampled.
    method: str {'inst', 'min', 'max', 'mean', 'median'}
        How the data should be resampled. If 'inst', it will be
        interpolated to the desired frequency.

    Returns
    -------
    pd.DataFrame
        Dataframe containing concatenated data, resampled or interpolated
        if required.
    """
    print(all_data_paths[0])
    data = _read_data(all_data_paths[0], resample_freq, method)
    if len(all_data_paths) > 1:
        for d in all_data_paths[1:]:
            print(d)
            data = data.merge(
                _read_data(d, resample_freq, method), how="outer", left_index=True, right_index=True
            )
    return data


def _read_data(data_path, resample_freq=None, method="mean"):
    """Loads the data. If inst_resample is not None, tries to interpolate
    to inst_resample frequency. If mean_resample is not None, tries to
    mean_resample using mean values.

    Parameters
    ----------
    data_path: str,
        str containing the path to the data.
    resample_freq: str | None
        Frequency of the resample. Uses pandas frequencies. If is None,
        the data will not be resampled.
    method: str {'inst', 'min', 'max', 'mean', 'median'}. Default: 'inst'
        How the data should be resampled. If 'inst', it will be
        interpolated to the desired frequency.

    Returns
    -------
    pd.DataFrame
        Dataframe containing data, resampled or interpolated if required.
    """
    data = pd.read_csv(data_path, header=[0, 1], index_col=0, parse_dates=True, comment="#")

    if resample_freq is None:
        return data

    resampled = data.resample(resample_freq)
    if method == "inst":
        sampled_data = resampled.interpolate(method="linear")
    elif method == "min":
        sampled_data = resampled.min()
    elif method == "max":
        sampled_data = resampled.max()
    elif method == "mean":
        sampled_data = resampled.mean()
    elif method == "median":
        sampled_data = resampled.median()
    else:
        warnings.warn(f"method {method} not supported, ignoring.")
        return data

    return sampled_data


def main():
    """If used as a CLI script, it reads and processes everything."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--coordinates", required=True)
    parser.add_argument("-p", "--path", required=True)
    parser.add_argument("-v", "--variables", default="nox")
    parser.add_argument("-r", "--resample_freq", default=None)
    parser.add_argument("-m", "--method", default="mean")
    parser.add_argument("-o", "--output", default="processed.nc")
    args = parser.parse_args()
    processed_data = create_file_with_coords(
        args.coordinates,
        args.path,
        args.variables.split(","),
        args.resample_freq,
        args.method,
    )
    processed_data.to_netcdf(args.output)


if __name__ == "__main__":
    main()
