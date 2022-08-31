# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#


def read_pkl(filename):
    """Function to read a pickle file containing part of the analysis class (models, obs, paired)

    Parameters
    ----------
    filename : type
        Description of parameter `filename`.
        
    Returns
    -------
    obj : type
        Description of returned object.

    """
    from joblib import load
    
    print('Reading:', filename)
    with open(filename, 'rb') as file:
        obj = load(file)
        
    return obj

def read_grouped_ncf(filename,xr_kws={}):
    """Function to read netcdf4 files containing a part of the analysis class (models, obs, paired).

    Parameters
    ----------
    filename : str
        Description of parameter `filename`.
    xr_kws : optional
            Additional keyword arguments for xr.open_dataset()
        
    Returns
    -------
    ds_dict : type
        Dict containing xarray datasets for each group.

    """
    import netCDF4
    import xarray as xr
    
    print('Reading:', filename)
    
    # Get a list of netCDF groups in the file
    ncf = netCDF4.Dataset(filename,mode='r')
    groups = ncf.groups.keys()
    ncf.close()
    
    obj = {}
    for group in groups:
        obj[group] = xr.open_dataset(filename,group=group,**xr_kws)
        
    return obj

def xarray_to_class(class_type,group_ds):
    """Remake dict containing driver class instances from dict of xarray datasets. Dict of xarray datasets must contain 
    global attribute that contains json formatted class attributes.

    Parameters
    ----------
    class_type : str
        One of 'model', 'pair' or 'observation'
    group_ds : dict
        dict containing xarray datasets from read_grouped_ncf.

    Returns
    -------
    class_dict
    """
    import json
    from melodies_monet import driver
    
    class_dict = {}
    for group in group_ds.keys():
        if class_type == 'pair':
            c=driver.pair()
        elif class_type == 'model':
            c=driver.model()
        elif class_type == 'observation':
            c=driver.observation()

        obj_dict = json.loads(group_ds[group].attrs['dict_json'])
        for attr in obj_dict.keys():
            setattr(c, attr, obj_dict[attr])
        c.obj = group_ds[group]
        class_dict[group]=c

    return class_dict
