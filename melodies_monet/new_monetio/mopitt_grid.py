""" MOPITT gridded data File reader """
import pandas as pd
import xarray as xr
import numpy as np
import datetime
from pathlib import Path                   # System agnostic paths
import h5py                                # For loading he5 files
import glob
import os



def getStartTime(filename):
    """Method to read the time in MOPITT level 3 hdf files.

    Parameters
    ----------
    filename : string or list
        filename is the path to the file

    Returns
    -------
    startTime """
    
    structure ='/HDFEOS/ADDITIONAL/FILE_ATTRIBUTES'
    #print("READING FILE " + inFileName)
    fName = os.path.basename(filename)

    try:
        inFile = h5py.File(filename,'r')
    except:
        print("ERROR: CANNOT OPEN " + filename)
        return 0
   
    grp = inFile[structure]
    k = grp.attrs
    startTimeBytes = k.get("StartTime",default="hello")
    startTime = pd.to_datetime(startTimeBytes[0], unit='s', origin='1993-01-01 00:00:00')
    #print("******************", startTime)
   
    try:
        inFile.close()
    except:
        print("ERROR CANNOT CLOSE " + filename)
        return 0
   
    return startTime


def loadAndExtractGridHDF(filename,varname):
    
    """Method to open MOPITT gridded hdf files.
    Masks data that is missing (truns into np.nan).

    Parameters
    ----------
    filename : string
        filename is the path to the file
    varname : string
        The variable to load from the MOPITT file

    Returns
    -------
    xarray.DataSet """
    
    he5_load = h5py.File(filename, mode='r')
    lat = he5_load["/HDFEOS/GRIDS/MOP03/Data Fields/Latitude"][:]
    lon = he5_load["/HDFEOS/GRIDS/MOP03/Data Fields/Longitude"][:]
    alt = he5_load["/HDFEOS/GRIDS/MOP03/Data Fields/Pressure2"][:]
    alt_short = he5_load["/HDFEOS/GRIDS/MOP03/Data Fields/Pressure"][:]
    
    #LAT-LON variables
    if varname=='column':
        data_loaded = he5_load["/HDFEOS/GRIDS/MOP03/Data Fields/RetrievedCOTotalColumnDay"][:]
    elif varname=='apriori_col':
        data_loaded = he5_load["/HDFEOS/GRIDS/MOP03/Data Fields/APrioriCOTotalColumnDay"][:]
    elif varname=='apriori_surf':
        data_loaded = he5_load["/HDFEOS/GRIDS/MOP03/Data Fields/APrioriCOSurfaceMixingRatioDay"][:]
    elif varname=='pressure_surf':
        data_loaded = he5_load["/HDFEOS/GRIDS/MOP03/Data Fields/SurfacePressureDay"][:]
    #LAT-LON-ALT variables    
    elif varname=='ak_col':
        data_loaded = he5_load["/HDFEOS/GRIDS/MOP03/Data Fields/TotalColumnAveragingKernelDay"][:]
    elif varname=='apriori_prof':
        data_loaded = he5_load["/HDFEOS/GRIDS/MOP03/Data Fields/APrioriCOMixingRatioProfileDay"][:]
    
    #DEBEG
    #print(data_loaded.shape)

    # create xarray DataArray
    if (varname=='column' or varname=='apriori_col'
        or varname=='apriori_surf'or varname=='pressure_surf'):
        datarr_new = xr.DataArray(data_loaded, dims=["lon","lat"], coords=[lon,lat])
    elif (varname=='ak_col'):
        datarr_new = xr.DataArray(data_loaded, dims=["lon","lat","alt"], coords=[lon,lat,alt])
    elif (varname=='apriori_prof'):
        datarr_new = xr.DataArray(data_loaded, dims=["lon","lat","alt"], coords=[lon,lat,alt_short])
    
    # missing value -> nan
    ds_masked = datarr_new.where(datarr_new != -9999.)
    he5_load.close()
    
    return ds_masked


def readMOPITTvar(files, varname):
    """Loop through files to open the MOPITT level 3 data.

    Parameters
    ----------
    files : string or list of strings
        The full path to the file or files. Can take a file template with wildcard (*) symbol.
    varname : string
        The variable to load from the MOPITT file

    Returns
    -------
    xarray.DataSet """
    count = 0
    filelist = sorted(glob.glob(files, recursive=False))
    
    for filename in filelist:
        print(filename)
        data = loadAndExtractGridHDF(filename, varname)
        time = getStartTime(filename)
        data = data.expand_dims(axis=0, time=[time])
        if count == 0:
            data_array = data
            count += 1
        else:
            data_array = xr.concat([data_array, data], 'time')
            
    return data_array