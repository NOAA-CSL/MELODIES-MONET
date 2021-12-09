import os
import sys
import logging
from glob import glob

import numpy as np
import xarray as xr

sys.path.append('../../new_monetio')
from hdfio import hdf_open, hdf_close, hdf_read


def read_dataset(fname, variable_dict):
    """
    Parameters
    __________
    fname : str
        Input file path.

    Returns
    _______
    xarray.Dataset
    """
    print('reading ' + fname)

    ds = xr.Dataset()

    f = hdf_open(fname)
    latitude = hdf_read(f, 'Latitude')
    longitude = hdf_read(f, 'Longitude')
    start_time = hdf_read(f, 'Scan_Start_Time')
    for varname in variable_dict:
        print(varname)
        values = hdf_read(f, varname)
        if 'scale' in variable_dict[varname]:
            values = variable_dict[varname]['scale'] \
                * values
        values[values < 0] = np.nan
        ds[varname] = xr.DataArray(values)
    hdf_close(f)

    return ds


def read_mfdataset(fnames, variable_dict, debug=False):
    """
    Parameters
    __________
    fnames : str
        Regular expression for input file paths.

    Returns
    _______
    xarray.Dataset
    """
    if debug:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.INFO
    logging.basicConfig(stream=sys.stdout, level=logging_level)

    for subpath in fnames.split('/'):
        if '$' in subpath:
            envvar = subpath.replace('$', '')
            envval = os.getenv(envvar)
            if envval is None:
                print('environment variable not defined: ' + subpath)
                exit(1)
            else:
                fnames = fnames.replace(subpath, envval)

    print(fnames)
    files = sorted(glob(fnames))
    for file in files:
        granule = read_dataset(file, variable_dict)
