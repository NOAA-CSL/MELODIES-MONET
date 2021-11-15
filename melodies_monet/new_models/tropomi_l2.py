'''Readers for TROPOMI Level 2 data''' 


def read_dataset(fname,convert_units=True):
    """Method to open TROPOMI CO netcdf files.
    Parameters
    ----------
    fname : string
        fname is the path to the file.
    convert_units : boolean
        If true the units of the chemical species will be converted to mol/cm^2.
    Returns
    -------
    xarray.DataSet
    """
    
    from satpy.scene import Scene
    
    scn = Scene(reader = 'tropomi_l2', filenames = [fname])
    scn.load(['carbonmonoxide_total_column','column_averaging_kernel','pressure_levels','time_utc','qa_value'])
    
    satdat = scn.to_xarray_dataset()
    satdat['carbonmonoxide_total_column'] *= 6.02214e19
    satdat['carbonmonoxide_total_column'] *= 1.0e-18
    ndat,nl = satdat.longitude.shape
    satdat = satdat.rename({'carbonmonoxide_total_column':'CO'})
    satdat = _fix_sat_dims(satdat,ndat,nl)
    return satdat

def read_mfdataset(fnames,convert_units=True):
    """Method to open TROPOMI CO netcdf files.
    Parameters
    ----------
    fnames : list
        fnames is the list of paths to the files.
    convert_units : boolean
        If true the units of the chemical species will be converted to mol/cm^2.
    Returns
    -------
    xarray.DataSet
    """
    
    from satpy.scene import Scene
    
    scn = Scene(reader = 'tropomi_l2', filenames = fnames)
    
    scn.load(['carbonmonoxide_total_column','column_averaging_kernel','pressure_levels','time_utc','qa_value'])
    
    satdat = scn.to_xarray_dataset()
    satdat['carbonmonoxide_total_column'] *= 6.02214e19
    satdat['carbonmonoxide_total_column'] *= 1.0e-18
    ndat,nl = satdat.longitude.shape
    satdat = satdat.rename({'carbonmonoxide_total_column':'CO'})
    satdat = _fix_sat_dims(satdat,ndat,nl)
    return satdat

def _fix_sat_dims(data_array,ndat,nl):  
    '''Internal helper function for correcting dimension name issues with TROPOMI L2 data and adds 2d time field'''
    import numpy as np
    stop = data_array.end_time
    start = data_array.start_time
    
    td = (stop-start)/ndat
    del stop
    time = np.asarray([start+i*td for i in range(ndat)])
    del start
    time = np.full((ndat,nl),np.expand_dims(time,axis=1))
    data_array['time'] = (('y','x'),time)
    data_array = data_array.set_coords(['time'])
    del time
    try:
        data_array['pressure_levels'] /= 100. 
    except: pass
    data_array = data_array.rename({'layer':'z'})
    return data_array