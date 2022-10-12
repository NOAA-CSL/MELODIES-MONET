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

def read_analysis_ncf(filenames,xr_kws={}):
    """Function to read netcdf4 files containing an object within an attribute of a part of the
    analysis class (models, obs, paired). For example, a single model/obs pairing or a single model. 
    If the object is saved in multiple files, the function will merge the files.

    Parameters
    ----------
    filenames : str or iterable
        Description of parameter `filename`.
    xr_kws : optional
        Additional keyword arguments for xr.open_dataset()
        
    Returns
    -------
    ds_out : type
        Xarray dataset containing merged files.

    """
    import xarray as xr
    
    if len(filenames)==1:
        ds_out = xr.open_dataset(filenames[0],**xr_kws)
        
    elif len(filenames)>1:
        for count, file in enumerate(filenames):
            print('Reading:', file)

            if count==0:
                ds_out = xr.open_dataset(file,**xr_kws)
            else:
                ds_append = xr.open_dataset(file,**xr_kws)
                ds_out = xr.merge([ds_out,ds_append])
            
    return ds_out

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
