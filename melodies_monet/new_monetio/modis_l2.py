import os
import sys
from glob import glob
import xarray as xr
sys.path.append('../../new_monetio')
from hdfio import hdf_open, hdf_close, hdf_read


def read_dataset(fname, variable_dict):
    """
    Parameters
    __________
    fname : str
        Input file path.
    """
    print('reading ' + fname)
    f = hdf_open(fname)
    latitude = hdf_read(f, 'Latitude')
    longitude = hdf_read(f, 'Longitude')
    for varname in variable_dict:
        print(varname)
        values = hdf_read(f, varname)
    hdf_close(f)


def read_mfdataset(fnames, variable_dict):
    """
    Parameters
    __________
    fnames : str
        Regular expression for input file paths.

    Returns
    _______
    xarray.Dataset
    """
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
        read_dataset(file, variable_dict)
