""" CESM File Reader """
import xarray as xr
from numpy import array, concatenate
from pandas import Series, to_datetime


def open_mfdataset(fname,
                   earth_radius=6370000,
                   convert_to_ppb=True,
                   var_list = ['O3','NO','NO2','lat','lon'],
                   scrip_file = '',
                   **kwargs):
    """Method to open multiple (or single) CESM SE netcdf files.
       This method extends the xarray.open_mfdataset functionality
       It is the main method called by the driver. Other functions defined 
       in this file are internally called by open_mfdataset and are preceeded
       by an underscore (e.g. _get_latlon).

    Parameters
    ----------
    fname : string or list
        fname is the path to the file or files.  It will accept wildcards in
        strings as well.
    earth_radius : float
        The earth radius used for map projections
    convert_to_ppb : boolean
        If true the units of the gas species will be converted to ppbV
        and units of aerosols to ug m^-3
    var_list : string or list
        List of variables to load from the CESM file. Default is to load ozone (O3) and PM2.5 (PM25).
    scrip_file: string
        Scrip file path for unstructured grid output


    Returns
    -------
    xarray.DataSet


    """
    # check that the files are netcdf format
    names, netcdf = _ensure_mfdataset_filenames(fname)
    
    # open the dataset using xarray
    try:
        if netcdf:
            dset_load = xr.open_mfdataset(fname, combine='by_coords', concat_dim='time', **kwargs)
        else:
            raise ValueError
    except ValueError:
        print('''File format not recognized. Note that files should be in netcdf
                format. Do not mix and match file types.''')
    
    # To keep lat & lon variables in the dataset
    if 'lat' not in var_list:
        var_list.append( 'lat' )
    if 'lon' not in var_list:
        var_list.append( 'lon' )

    
    # ===========================
    # Process the loaded data
    #extract variables of choice
    dset = dset_load.get(var_list)
    #rename altitude variable to z for monet use
    dset = dset.rename({'lev':'z'})
    #re-order so surface is associated with the first vertical index
    dset = dset.sortby('z',ascending=False)
    # ===========================
    
        
    # Make sure this dataset has unstructured grid
    dset.attrs["mio_has_unstructured_grid"] = True
    dset.attrs["mio_scrip_file"] = scrip_file
    
    
    # convert units
    if convert_to_ppb:
        for i in dset.variables:
            if 'units' in dset[i].attrs:
                # convert all gas species from mol/mol to ppbv
                if 'mol/mol' in dset[i].attrs['units']:
                    dset[i] *= 1e09
                    dset[i].attrs['units'] = 'ppbV'
                # convert 'kg/m3 to \mu g/m3 '
                elif 'kg/m3' in dset[i].attrs['units']:
                    dset[i] *= 1e09
                    dset[i].attrs['units'] = '$\mu g m^{-3}$'
                    

    # dset_scrip = xr.open_dataset( scrip_file )                
    # return dset, dset_scrip
    return dset


"""
-----------------------------------------
Below are internal functions to this file
-----------------------------------------
"""

def _ensure_mfdataset_filenames(fname):
    '''Checks if dataset in netcdf format
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