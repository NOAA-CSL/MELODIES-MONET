# Reader for FV3-RAQMS, a chemistry variant of UFS run at UW-Madison.

# Reader for FV3-RAQMS, a chemistry variant of UFS.

import xarray as xr

def open_dataset(fname):
    '''Open a single dataset from fv3-raqms output. Currently set for netcdf file format.
    Parameters
    ----------
    fname : string
        Filename to be opened.
    Returns
    -------
    xarray.Dataset
    '''
    names, netcdf = _ensure_mfdataset_filenames(fname)
    try:
        if netcdf:
            f = xr.open_dataset(names[0])
            f = _fix_grid(f)
            f = _fix_time(f)
        else:
            raise ValueError
    except ValueError:
        print('''File format not recognized. Note that files should be
                preprocessed to netcdf.''')

    return f

def open_mfdataset(fname):
    '''Open a multiple file dataset from fv3-raqms output.
    Parameters
    ----------
    fname : string
        Filenames to be opened
    Returns
    -------
    xarray.Dataset
    '''
    names, netcdf = _ensure_mfdataset_filenames(fname)
    try:
        if netcdf:
            f = xr.open_mfdataset(names, concat_dim='time')
            f = _fix_grid(f)
            f = _fix_time(f)
        else:
            raise ValueError
    except ValueError:
        print('''File format not recognized. Note that files should be in netcdf
                format. Do not mix and match file types.''')

    return f

def _ensure_mfdataset_filenames(fname):
    '''Checks if netcdf dataset
    Parameters
    ----------
    fname : string or list of strings
    Returns
    -------
    type
    '''
    from glob import glob
    from numpy import sort
    import six

    if isinstance(fname, six.string_types):
        names = sort(glob(fname))
    else:
        names = sort(fname)
    netcdfs = [True for i in names if 'nc' in i]
    netcdf = False
    if len(netcdfs) >= 1:
        netcdf = True
    return names, netcdf

def _fix_grid(f):
    from  numpy import meshgrid
    f = f.rename({'grid_yt':'lat','grid_xt':'lon'})
    lat = f.lat.values
    lon = f.lon.values
    lon[(lon > 180)] -= 360
    lon,lat = meshgrid(lon,lat)
    f = f.rename({'lat':'y','lon':'x','grid_zt':'z'})
    f['longitude']=(('y','x'),lon)
    f['latitude']=(('y','x'),lat)
    f = f.set_coords(['latitude','longitude'])
    del lon,lat
    return f

def _fix_time(f):
    import pandas as pd
    dtstr = f.time
   
    date = pd.to_datetime(dtstr,unit='D',origin=pd.Timestamp('2018-12-31'))
    f['Time'] = (('time',),date)
    
    f = f.set_index({'time':'Time'})
    del date,dtstr
    return f