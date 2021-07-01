""" RAPchem File reader """
import xarray as xr
from numpy import array, concatenate
from pandas import Series, to_datetime

from monet.monet_accessor import _dataset_to_monet
from wrf import getvar, ALL_TIMES, extract_global_attrs
from netCDF4 import Dataset
import dask.array as da

def can_do(index):
    if index.max():
        return True
    else:
        return False

def open_mfdataset(fname,
                   earth_radius=6370000,
                   convert_to_ppb=True,
                   drop_duplicates=False,
                   **kwargs):
    """Method to open RAP-chem netcdf files.

    Parameters
    ----------
    fname : string or list
        fname is the path to the file or files.  It will accept hot keys in
        strings as well.
    earth_radius : float
        The earth radius used for the map projection
    convert_to_ppb : boolean
        If true the units of the gas species will be converted to ppbV
        and units of aerosols to ug m^-3

    Returns
    -------
    xarray.DataSet


    """
    wrflist = []
    for files in fname:
        wrflist.append(Dataset(files))
    
    var_list = ['o3','no','no2','so2','co','PM2_5_DRY']
               #'pres','tk']
   
    var_wrf_list = []
    for var in var_list:
        if var == 'pres': #Insert special versions.
            var_wrf = getvar(wrflist,var,timeidx=ALL_TIMES,method='cat',squeeze=False,units='Pa')
        else:
            var_wrf = getvar(wrflist,var,timeidx=ALL_TIMES,method='cat',squeeze=False)
        var_wrf_list.append(var_wrf)

    dset = xr.merge(var_wrf_list)
  
    #sum nox
    dset = add_lazy_nox(dset)

    #Add global attributes needed
    a_truelat1 = extract_global_attrs(wrflist[0],'TRUELAT1')    
    a_truelat2 = extract_global_attrs(wrflist[0],'TRUELAT2')
    a_moad_cen_lat = extract_global_attrs(wrflist[0],'MOAD_CEN_LAT')
    a_stand_lon = extract_global_attrs(wrflist[0],'STAND_LON')
    a_map_proj = extract_global_attrs(wrflist[0],'MAP_PROJ')
    a_cen_lat = extract_global_attrs(wrflist[0],'CEN_LAT')
    a_cen_lon = extract_global_attrs(wrflist[0],'CEN_LON')
    dset = dset.assign_attrs(a_truelat1)
    dset = dset.assign_attrs(a_truelat2)
    dset = dset.assign_attrs(a_moad_cen_lat)
    dset = dset.assign_attrs(a_stand_lon)
    dset = dset.assign_attrs(a_map_proj)
    dset = dset.assign_attrs(a_cen_lat)
    dset = dset.assign_attrs(a_cen_lon)

    for i in dset.variables:
        dset[i] = dset[i].assign_attrs(a_truelat1)
        dset[i] = dset[i].assign_attrs(a_truelat2)
        dset[i] = dset[i].assign_attrs(a_moad_cen_lat)
        dset[i] = dset[i].assign_attrs(a_stand_lon)
        dset[i] = dset[i].assign_attrs(a_map_proj)
        dset[i] = dset[i].assign_attrs(a_cen_lat)
        dset[i] = dset[i].assign_attrs(a_cen_lon)
    #convert to monet format
    dset = _dataset_to_monet(dset)

    #Time is already in correct format but need to rename
    dset = dset.rename({'Time': 'time'})
    dset = dset.rename(dict(bottom_top='z'))
    grid = ioapi_grid_from_dataset_wrf(dset, earth_radius=earth_radius)
    dset = dset.assign_attrs({'proj4_srs': grid})
    for i in dset.variables:
        dset[i] = dset[i].assign_attrs({'proj4_srs': grid})
    
    # convert all gas species to ppbv
    if convert_to_ppb:
        for i in dset.variables:
            if 'units' in dset[i].attrs:
                if 'ppmv' in dset[i].attrs['units']:
                    dset[i] *= 1000.
                    dset[i].attrs['units'] = 'ppbV'
    # convert 'ug/kg-dryair -> ug/m3'
    for i in dset.variables:
        if 'units' in dset[i].attrs:
            if 'ug/kg-dryair' in dset[i].attrs['units']:
                # ug/kg -> ug/m3 using dry air density
                dset[i] = dset[i]*dset['pressure']/dset['temp']/287.05535 
                dset[i].attrs['units'] = '$\mu g m^{-3}$'

    #assign mapping table for airnow
    dset = _predefined_mapping_tables(dset)
    return dset 

def ioapi_grid_from_dataset_wrf(ds, earth_radius=6370000):
    """SGet the IOAPI projection out of the file into proj4.
    Parameters
    ----------
    ds : type
        Description of parameter `ds`.
    earth_radius : type
        Description of parameter `earth_radius`.
    Returns
    -------
    type
        Description of returned object.
    """

    pargs = dict()
    pargs['lat_1'] = ds.TRUELAT1
    pargs['lat_2'] = ds.TRUELAT2
    pargs['lat_0'] = ds.MOAD_CEN_LAT
    pargs['lon_0'] = ds.STAND_LON
    #pargs['center_lon'] = ds.XCENT
    #pargs['x0'] = ds.XORIG
    #pargs['y0'] = ds.YORIG
    pargs['r'] = earth_radius
    proj_id = ds.MAP_PROJ
    if proj_id == 6:
        # Cylindrical Equidistant -RAP
        #Assign lat_ts = 0 as true lat are missing value in RAP-chem
        p4 = '+proj=eqc +lat_ts=0' \
             '+lat_0={lat_0} +lon_0={lon_0} ' \
             '+x_0=0 +y_0=0 +datum=WGS84 +units=m +a={r} +b={r}'
        p4 = p4.format(**pargs)
    elif proj_id == 1:
        # lambert - WRF
        #Described here: https://fabienmaussion.info/2018/01/06/wrf-projection/
        p4 = '+proj=lcc +lat_1={lat_1} +lat_2={lat_2} ' \
             '+lat_0={lat_0} +lon_0={lon_0} ' \
             '+x_0=0 +y_0=0 +datum=WGS84 +units=m +a={r} +b={r}'
        p4 = p4.format(**pargs)    
    #elif proj_id == 4:
    #    # Polar stereo
    #    p4 = '+proj=stere +lat_ts={lat_1} +lon_0={lon_0} +lat_0=90.0' \
    #         '+x_0=0 +y_0=0 +a={r} +b={r}'
    #    p4 = p4.format(**pargs)
    #elif proj_id == 3:
    #    # Mercator
    #    p4 = '+proj=merc +lat_ts={lat_1} ' \
    #         '+lon_0={center_lon} ' \
    #         '+x_0={x0} +y_0={y0} +a={r} +b={r}'
    #    p4 = p4.format(**pargs)
    else:
        raise NotImplementedError('IOAPI proj not implemented yet: '
                                  '{}'.format(proj_id))
    # area_def = _get_ioapi_pyresample_area_def(ds)
    return p4  # , area_def

def _get_keys(d):
    keys = Series([i for i in d.data_vars.keys()])
    return keys


def add_lazy_pm25(d):
    """Short summary.
    #Note this is set up for RAP-chem mechanism edit this for your version WRF-chem mechanism.

    Parameters
    ----------
    d : type
        Description of parameter `d`.

    Returns
    -------
    type
        Description of returned object.

    """
    keys = _get_keys(d)
    allvars = Series(['so4aj','so4ai','nh4aj','nh4ai','no3aj','no3ai','naaj','naai',
        'claj','clai','poa0j','poa0i','poa1j','poa1i','poa2j','poa2i','poa3j','poa3i','asoa0j','asoa0i',
        'asoa1j','asoa1i','asoa2j','asoa2i','asoa3j','asoa3i','bsoa1j','bsoa1i','bsoa2j','bsoa2i',
        'bsoa3j','bsoa3i','ecj','eci','p25j','p25i'])
    index = allvars.isin(keys)
    #RHS add extra check that if not all variables in output you do not print PM2.5
    if can_do(index) and all(index):
        newkeys = allvars.loc[index]
        d['PM25_calc'] = add_multiple_lazy(d, newkeys)
        d['PM25_calc'] = d['PM25_calc'].assign_attrs({
            'units': 'ug/kg-dryair','name': 'PM2.5_calc','long_name': 'PM2.5_calc_monet'})
    return d

def add_lazy_nox(d):
    keys = _get_keys(d)
    allvars = Series(['no', 'no2'])
    index = allvars.isin(keys)
    #RHS add extra check that if not all variables do not print
    if can_do(index) and all(index):
        newkeys = allvars.loc[index]
        d['NOx'] = add_multiple_lazy(d, newkeys)
        d['NOx'] = d['NOx'].assign_attrs({'units': 'ppmv','name': 'NOx', 'long_name': 'NOx'})
    return d


def add_multiple_lazy(dset, variables, weights=None):
    from numpy import ones
    if weights is None:
        weights = ones(len(variables))
    else:
        weights = weights.values
    variables = variables.values
    new = dset[variables[0]].copy() * weights[0]
    for i, j in zip(variables[1:], weights[1:]):
        new = new + dset[i] * j
    return new


def _predefined_mapping_tables(dset):
#for now I adapt this to only output what I need and that is to to_airnow so can use in Patrick's scripts
    """Predefined mapping tables for different observational parings used when
        combining data.

    Returns
    -------
    dictionary
        A dictionary of to map to.

    """
    to_airnow = {
        'OZONE': 'o3',
        'PM2.5': 'PM2_5_DRY',
        'CO': 'CO',
        'NOX': 'NOx',
        'SO2': 'SO2',
        'NO': 'NO',
        'NO2': 'NO2',
    }
    dset = dset.assign_attrs({'mapping_tables_to_airnow': to_airnow})
    return dset
