import os
from glob import glob
import xarray as xr
import h5py as h5


def read_dataset(fname):
    pass


def read_mfdataset(fnames):
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
    files = glob(fnames)
    for file in files:
        print(file)
