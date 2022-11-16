# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#
import numpy as np
from pandas.api.types import is_float_dtype

def write_analysis_ncf(obj, output_dir='', fn_prefix=None, keep_groups=None, title=''):
    """Function to write netcdf4 files with some compression for floats from an attribute of the
    analysis class (models, obs, paired). Writes the objects within the attribute as separate files.

    Parameters
    ----------
    obj : dict
        Dict containing attribute of the driver analysis class (model, observation or paired).
    output_dir : str
        Directory to save files in.
    fn_prefix : str
        Prefix to add to the group name when saving.
    keep_groups : list
        List of groups that should be saved. Other groups will not be saved.

    Returns
    -------
    None

    """
    import pandas as pd
    import json
    import os.path
    
    if fn_prefix is not None:
        base_name = os.path.join(output_dir,'{prefix}_{groupname}.nc4')
    else:
        base_name = os.path.join(output_dir,'{groupname}.nc4')
    
    if keep_groups is not None:
        groups=[elem for elem in obj.keys() if elem in keep_groups]
    else:
        groups=obj.keys()
    
    for group in groups:
        output_name = base_name.format(prefix=fn_prefix, groupname=group)
        print('Writing:', output_name)
        
        dset=obj[group].obj
        comp = dict(zlib=True, complevel=7)
        encoding = {}
        for i in dset.data_vars.keys():
            # if is_float_dtype(dset[i]):  # (dset[i].dtype != 'object') & (i != 'time') & (i != 'time_local') :
            #     print("Compressing: {}, original_dtype: {}".format(i, dset[i].dtype))
            #     dset[i] = compress_variable(dset[i])
            encoding[i] = comp
        dset.attrs['title'] = title
        dset.attrs['format'] = 'NetCDF-4'
        dset.attrs['date_created'] = pd.to_datetime('today').strftime('%Y-%m-%d')
        dict_json = obj[group].__dict__.copy()
        dict_json.pop('obj')
        dset.attrs['dict_json'] = json.dumps(dict_json, indent = 4) 
        dset.attrs['group_name'] = group
        dset.to_netcdf(output_name)

def write_ncf(dset, output_name, title=''):
    """Function to write netcdf4 files with some compression for floats

    Parameters
    ----------
    dset : type
        Description of parameter `dset`.
    output_name : type
        Description of parameter `output_name`.

    Returns
    -------
    type
        Description of returned object.

    """
    import pandas as pd

    print('Writing:', output_name)
    comp = dict(zlib=True, complevel=7)
    encoding = {}
    for i in dset.data_vars.keys():
        if is_float_dtype(dset[i]):  # (dset[i].dtype != 'object') & (i != 'time') & (i != 'time_local') :
            print("Compressing: {}, original_dtype: {}".format(i, dset[i].dtype))
            dset[i] = compress_variable(dset[i])
        encoding[i] = comp
    dset.attrs['title'] = title
    dset.attrs['format'] = 'NetCDF-4'
    dset.attrs['date_created'] = pd.to_datetime('today').strftime('%Y-%m-%d')
    dset.to_netcdf(output_name, encoding=encoding)


def compute_scale_and_offset(mn, mx, n, dtype=np.float32):
    """Calculates the scale and offset to be used for a variable

    Parameters
    ----------
    mn : float
        minimum value.
    mx : float
        maximum value.
    n : number of bits
        default is 32bit.
    dtype : numpy dtype
        default is numpy.float32.

    Returns
    -------
    type
        Description of returned object.

    """
    """
    min is the minimum of the values
    max is the maximum of the values
    n is the integer bit length  (ie 32 for np.int32 or 16 for np.int16)
    """
    # stretch/compress data to the available packed range
    scale_factor = (mx - mn) / (2 ** n - 1)
    # translate the range to be symmetric about zero
    add_offset = mn + 2 ** (n - 1) * scale_factor
    return (scale_factor.astype(dtype), add_offset.astype(dtype))


def pack_value(values, scale_factor, offset, dtype):
    """Values to pack the array with scale factors from a float to integers

    Parameters
    ----------
    values : type
        Description of parameter `values`.
    scale_factor : type
        Description of parameter `scale_factor`.
    offset : type
        Description of parameter `offset`.
    dtype : type
        Description of parameter `dtype`.

    Returns
    -------
    type
        Description of returned object.

    """
    return ((values - offset) / scale_factor).astype(dtype)


def get_min_max(da):
    """Function to return the maximum and minimum value

    Parameters
    ----------
    da : type
        Description of parameter `da`.

    Returns
    -------
    type
        Description of returned object.

    """
    return (da.min().compute(), da.max().compute())


def compress_variable(da):
    """Function to compress a variable from a float to integer and adds netcdf attributes for CF convention.

    Parameters
    ----------
    da : type
        Description of parameter `da`.

    Returns
    -------
    type
        Description of returned object.

    """
    da = da.fillna(-1)
    mn, mx = get_min_max(da)
    scale_factor, offset = compute_scale_and_offset(mn, mx, 32, dtype=da.dtype)
    da.data = pack_value(da, scale_factor, offset, dtype=np.int32).data
    da.attrs['scale_factor'] = scale_factor.values
    da.attrs['add_offset'] = offset.values
    da.attrs['_FillValue'] = -1
    da.attrs['missing_value'] = -1
    return da

def write_pkl(obj, output_name):
    """Function to write a pickle file from an attribute of the analysis class (models, obs, paired)

    Parameters
    ----------
    obj : type
        Description of parameter `obj`.
    output_name : str
        Description of parameter `output_name`.
        
    Returns
    -------
    None

    """
    from joblib import dump
    
    print('Writing:', output_name)
    with open(output_name, 'wb') as outp:
        dump(obj, outp)
