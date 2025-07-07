import datetime as dt
from glob import glob

import pandas as pd
import xarray as xr


def _parse_metadata(value):
    """Parse metadata to possible values.

    Parameters
    ----------
    value : str
        str to parse

    Returns
    -------
    int | float | datetime | str
        parsed data
    """
    if value == "":  # Deal with the empty string as a special case
        return value
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        pass
    try:
        return pd.to_datetime(value, format="ISO8601").to_datetime64()
    except ValueError:
        return value.lstrip().rstrip()


def _name_columns(index_number):
    """Creates the names for columns.

    Parameters
    ----------
    index_number : int
        Number originally from the df index

    Returns
    -------
    str
        str as Column {index_number + 1}
    """
    return f"Column {index_number + 1}"


def _rename_and_format(df):
    """Renames each variable with the column number, following the convention of PGN files.
    It formats the data as an xarray.Dataset and adds the x dimension

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe with columns to rename

    Returns
    -------
    xr.Dataset
        Dataset with the renamed data
    """
    df2 = df.rename(_name_columns, axis="columns")
    df2 = df2.rename(columns={"Column 1": "time"})
    df2 = df2.set_index("time")
    print("Calculating pandora hourly averages")
    df2 = df2.resample('h').mean()
    ds = xr.Dataset.from_dataframe(df2).expand_dims("x", axis=1).compute()
    return ds.copy(deep=True)


def _read_pandora_file(file_path):
    """Reads in a Pandora data file and formats it correctly as a xr.Dataset

    Parameters
    ----------
    file_path : str
        String with the path to a single Pandora PGN text file

    Returns
    -------
    xr.Dataset
        Dataset from single file formatted for MELODIES-MONET
    """
    count_line_dividers = 0
    headers = {}
    global_attrs = {
        "history": f"{dt.datetime.now()}: created from _read_pandora_files, pandora_pgn.py"
    }
    data_collection = []
    with open(file_path, encoding="latin-1") as f:
        for line in f:
            line_stripped = line.rstrip()
            if line_stripped.startswith("-----------"):
                count_line_dividers += 1
            elif count_line_dividers == 0:
                attr_name, value = tuple(line_stripped.split(":"))
                value = _parse_metadata(value)
                global_attrs[attr_name] = value
            elif count_line_dividers == 1:
                key, metadata = tuple(line_stripped.split(":"))
                headers[key] = metadata
            elif count_line_dividers == 2:
                data_collection.append(line_stripped.split())
    _df = pd.DataFrame(data_collection)
    times = pd.to_datetime(_df[0], format="ISO8601").dt.tz_localize(None)
    # errors = corece turns non valid strings into NaN
    measurements = _df.loc[:, _df.columns != 0].apply(
        pd.to_numeric, errors="coerce", downcast="float"
    )
    df = pd.concat([times, measurements], axis=1)
    data = _rename_and_format(df)
    data["latitude"] = (("x",), [global_attrs["Location latitude [deg]"]])
    data["latitude"].attrs["units"] = "degrees_north"
    data["longitude"] = (("x",), [global_attrs["Location longitude [deg]"]])
    data["longitude"].attrs["units"] = "degrees_east"
    data.attrs = global_attrs
    for k in headers:
        if k in data:
            data[k].attrs["description"] = headers[k]
        elif k.startswith("From Column"):
            optional_keys = headers[k]
    non_shared_keys = list(set(data.keys()) - set(headers.keys()))
    non_shared_keys.remove("latitude")
    non_shared_keys.remove("longitude")
    if len(non_shared_keys) > 0:
        for k in non_shared_keys:
            data[k].attrs["description"] = optional_keys
    data["siteid"] = (("x",), [data.attrs["Short location name"]])
    data = data.assign_coords({"longitude": data["longitude"], "latitude": data["latitude"]})
    return data


def _merge_global_attrs(ds1, ds2, merged):
    """Merges global attributes of two datasets as a list and adds
    them to the merged dataset inplace

    Parameters
    ----------
    ds1 : xr.Dataset
        First dataset
    ds2 : xr.Dataset
        Second dataset
    merged : xr.Dataset
        merged Dataset, the global attributes will be assigned as
        lists of ds1.attrs + ds2.attrs

    Returns
    -------
    None
    """
    for k in merged.attrs:
        merged.attrs[k] = list([ds1.attrs.get(k, "")]) + list([ds2.attrs.get(k, "")])


def open_mfdataset(path):
    """Opens multiple Pandora PGN files and combines them to
    MELODIES-MONET compatible format.

    Parameters
    ----------
    path: str
        String containing the paths

    Returns
    -------
    xr.Dataset
        Formatted dataset. Should work for MELODIES-MONET.
    """
    print("Importing pandora data. The data will be presented as hourly means.")
    if isinstance(path, str):
        files = sorted(glob(path))
    if isinstance(path, list):
        files = []
        for file in path:
            files = files + list(glob(str(file)))
        files = sorted(files)
    ds = _read_pandora_file(files[0])
    if len(files) > 1:
        for f in files[1:]:
            ds2 = _read_pandora_file(f)
            if ds.attrs["Data file version"] != ds2.attrs["Data file version"]:
                raise Exception("Different data file versions, cannot concatenate")
            if ds.attrs["Short location name"] != ds2.attrs["Short location name"]:
                ds = xr.concat([ds, ds2], dim="x")
            else:
                ds = xr.concat([ds, ds2], dim="time")
            _merge_global_attrs(ds, ds2, ds)
        ds.attrs["history"] = [
            f"{dt.datetime.now()}: open_mfdataset from pandora_pgn.py "
        ] + ds.attrs["history"]
    return ds
