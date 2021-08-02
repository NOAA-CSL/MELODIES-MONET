""" RRFS-CMAQ File Reader """
import xarray as xr
from numpy import array, concatenate
from pandas import Series, to_datetime
import numpy as np


def can_do(index):
    if index.max():
        return True
    else:
        return False

def open_mfdataset(fname, earth_radius=6370000, convert_to_ppb=True, drop_duplicates=False, **kwargs):
    """Method to open RFFS-CMAQ dyn* netcdf files.

    Parameters
    ----------
    fname : string or list
        fname is the path to the file or files.  It will accept hot keys in
        strings as well. This will work for 1 or multiple files.
    earth_radius : float
        The earth radius used for the map projection
    convert_to_ppb : boolean
        If true the units of the gas species will be converted to ppbV

    Returns
    -------
    xarray.DataSet


    """

    # open the dataset using xarray
    dset = xr.open_mfdataset(fname, concat_dim='time', **kwargs)
    
    #Standardize some variable names
    dset = dset.rename({'grid_yt': 'y','grid_xt': 'x','pfull': 'z',
                        'phalf': 'z_i', #Interface pressure levels
                        'lon': 'longitude','lat': 'latitude',
                        'tmp': 'temperature_k'}) #standard temperature (kelvin)

    #Calculate pressure. This has to go before sorting because ak and bk 
    #are not sorted as they are in attributes
    dset['pres_pa'] = _calc_pressure(dset)

    #Adjust pressure levels for all models such that the surface is first.
    dset = dset.sortby('z', ascending=False)
    dset = dset.sortby('z_i', ascending=False)
    
    #Note these pressure and altitude calcs have to go after resorting.
    dset['geohgt_m'] = _calc_hgt(dset)
    
    #Set coordinates
    dset = dset.reset_index(['x','y','z','z_i'])
    dset = dset.reset_coords()
    dset = dset.set_coords(['latitude','longitude','pres_pa','geohgt_m'])
    #time is already set as a coordinate correctly.
    
    #Need to adjust units before summing for aerosols
    # convert all gas species to ppbv
    if convert_to_ppb:
        for i in dset.variables:
            if 'units' in dset[i].attrs:
                if 'ppmv' in dset[i].attrs['units']:
                    dset[i] *= 1000.0
                    dset[i].attrs['units'] = 'ppbV'

    # convert 'ug/kg to ug/m3'
    for i in dset.variables:
        if 'units' in dset[i].attrs:
            if 'ug/kg' in dset[i].attrs['units']:
                # ug/kg -> ug/m3 using dry air density
                dset[i] = dset[i]*dset['pres_pa']/dset['temperature_k']/287.05535
                dset[i].attrs['units'] = '$\mu g m^{-3}$'
                
    # add lazy diagnostic variables
    dset = add_lazy_pm25(dset)
    dset = add_lazy_pm10(dset)
    dset = add_lazy_pm_course(dset)
    dset = add_lazy_clf(dset)
    dset = add_lazy_naf(dset)
    #dset = add_lazy_caf(dset) Leave out for now. Missing a variable and need to check.
    dset = add_lazy_noy(dset)
    dset = add_lazy_nox(dset)
    dset = add_lazy_no3f(dset)
    dset = add_lazy_nh4f(dset)
    dset = add_lazy_so4f(dset)

    # Change the times to pandas format
    dset['time'] = dset.indexes['time'].to_datetimeindex(unsafe=True)
    #Turn off warning for now. This is just because the model is in julian time

    return dset

def _get_keys(d):
    keys = Series([i for i in d.data_vars.keys()])
    return keys


def add_lazy_pm25(d):
    """Short summary.

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
    allvars = Series(concatenate([aitken, accumulation, coarse]))
    weights = Series(concatenate([np.ones(len(aitken)),
                                  np.ones(len(accumulation)),
                                  np.full(len(coarse),0.2)]))
    
    if 'PM25_TOT' in keys.to_list():
        d['PM25'] = d['PM25_TOT']
    else:
        index = allvars.isin(keys)
        if can_do(index):
            newkeys = allvars.loc[index]
            newweights = weights.loc[index]
            d['PM25_calc'] = add_multiple_lazy(d, newkeys, weights=newweights)
            d['PM25_calc'] = d['PM25_calc'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'PM2.5', 'long_name': 'PM2.5 calculated by MONET'})
    return d


def add_lazy_pm10(d):
    keys = _get_keys(d)
    allvars = Series(concatenate([aitken, accumulation, coarse]))
    if 'PMC_TOT' in keys.to_list():
        d['PM10'] = d['PMC_TOT']
    else:
        index = allvars.isin(keys)
        if can_do(index):
            newkeys = allvars.loc[index]
            d['PM10_calc'] = add_multiple_lazy(d, newkeys)
            d['PM10_calc'] = d['PM10_calc'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'PM10', 'long_name': 'Particulate Matter < 10 microns'})
    return d


def add_lazy_pm_course(d):
    keys = _get_keys(d)
    allvars = Series(coarse)
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        d['PM_COURSE'] = add_multiple_lazy(d, newkeys)
        d['PM_COURSE'] = d['PM_COURSE'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'PM_COURSE', 'long_name': 'Course Mode Particulate Matter'})
    return d


def add_lazy_clf(d):
    keys = _get_keys(d)
    allvars = Series(['acli', 'aclj', 'aclk'])
    weights = Series([1, 1, 0.2])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['CLf'] = add_multiple_lazy(d, newkeys, weights=neww)
        d['CLf'] = d['CLf'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'CLf', 'long_name': 'Fine Mode particulate Cl'})
    return d


def add_lazy_caf(d): #Missing ACAI leave out for now.
    keys = _get_keys(d)
    allvars = Series(['ACAI', 'ACAJ', 'ASEACAT', 'ASOIL', 'ACORS'])
    weights = Series([1, 1, 0.2 * 32.0 / 1000.0, 0.2 * 83.8 / 1000.0, 0.2 * 56.2 / 1000.0])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['CAf'] = add_multiple_lazy(d, newkeys, weights=neww)
        d['CAf'] = d['CAf'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'CAf', 'long_name': 'Fine Mode particulate CA'})
    return d


def add_lazy_naf(d):
    keys = _get_keys(d)
    allvars = Series(['anai', 'anaj', 'aseacat', 'asoil', 'acors'])
    weights = Series([1, 1, 0.2 * 837.3 / 1000.0, 0.2 * 62.6 / 1000.0, 0.2 * 2.3 / 1000.0])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['NAf'] = add_multiple_lazy(d, newkeys, weights=neww)
        d['NAf'] = d['NAf'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'NAf', 'long_name': 'NAf'})
    return d


def add_lazy_so4f(d):
    keys = _get_keys(d)
    allvars = Series(['aso4i', 'aso4j', 'aso4k'])
    weights = Series([1.0, 1.0, 0.2])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['SO4f'] = add_multiple_lazy(d, newkeys, weights=neww)
        d['SO4f'] = d['SO4f'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'SO4f', 'long_name': 'SO4f'})
    return d


def add_lazy_nh4f(d):
    keys = _get_keys(d)
    allvars = Series(['anh4i', 'anh4j', 'anh4k'])
    weights = Series([1.0, 1.0, 0.2])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['NH4f'] = add_multiple_lazy(d, newkeys, weights=neww)
        d['NH4f'] = d['NH4f'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'NH4f', 'long_name': 'NH4f'})
    return d


def add_lazy_no3f(d):
    keys = _get_keys(d)
    allvars = Series(['ano3i', 'ano3j', 'ano3k'])
    weights = Series([1.0, 1.0, 0.2])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['NO3f'] = add_multiple_lazy(d, newkeys, weights=neww)
        d['NO3f'] = d['NO3f'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'NO3f', 'long_name': 'NO3f'})
    return d


def add_lazy_noy(d):
    keys = _get_keys(d)
    allvars = Series(noy_gas)
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        d['NOy'] = add_multiple_lazy(d, newkeys)
        d['NOy'] = d['NOy'].assign_attrs({'name': 'NOy', 'long_name': 'NOy'})
    return d

def add_lazy_nox(d):
    keys = _get_keys(d)
    allvars = Series(['no', 'no2'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        d['NOx'] = add_multiple_lazy(d, newkeys)
        d['NOx'] = d['NOx'].assign_attrs({'name': 'NOx', 'long_name': 'NOx'})
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
    """Predefined mapping tables for different observational parings used when
        combining data.

    Returns
    -------
    dictionary
        A dictionary of to map to.

    """
    to_improve = {}
    to_nadp = {}
    to_aqs = {
        'OZONE': ['o3'],
        'PM2.5': ['PM25_calc'],
        'CO': ['co'],
        'NOY': ['NOy'],
        'NOX': ['NOx'],
        'SO2': ['so2'],
        'NO': ['no'],
        'NO2': ['no2'],
        'SO4f': ['SO4f'],
        'PM10': ['PM10_calc'],
        'NO3f': ['NO3f'],
        'ECf': ['ECf'],
        'OCf': ['OCf'],
        'ETHANE': ['eth'],
        'BENZENE': ['benzene'],
        'TOLUENE': ['tol'],
        'ISOPRENE': ['isop'],
        'O-XYLENE': ['XYL'], #maybe xylmn?
        'WS': ['WSPD10'], #Need to update
        'TEMP': ['temperature_k'],
        'WD': ['WDIR10'], #Need to update
        'NAf': ['NAf'],
        'MGf': ['AMGJ'],
        'TIf': ['ATIJ'],
        'SIf': ['ASIJ'],
        'Kf': ['Kf'],
        'CAf': ['CAf'], #Not calc now.
        'NH4f': ['NH4f'],
        'FEf': ['AFEJ'],
        'ALf': ['AALJ'],
        'MNf': ['AMNJ'],
    }
    to_airnow = {
        'OZONE': ['o3'],
        'PM2.5': ['PM25_calc'],
        'CO': ['co'],
        'NOY': ['NOy'],
        'NOX': ['NOx'],
        'SO2': ['so2'],
        'NO': ['no'],
        'NO2': ['no2'],
        'SO4f': ['SO4f'],
        'PM10': ['PM10_calc'],
        'NO3f': ['NO3f'],
        'ECf': ['ECf'],
        'OCf': ['OCf'],
        'ETHANE': ['eth'],
        'BENZENE': ['benzene'],
        'TOLUENE': ['tol'],
        'ISOPRENE': ['isop'],
        'O-XYLENE': ['XYL'], #Need to fix
        'WS': ['WSPD10'], #Need to fix
        'TEMP': ['temperature_k'],
        'WD': ['WDIR10'], #Need to find.
        'NAf': ['NAf'],
        'MGf': ['AMGJ'],
        'TIf': ['ATIJ'],
        'SIf': ['ASIJ'],
        'Kf': ['Kf'],
        'CAf': ['CAf'],
        'NH4f': ['NH4f'],
        'FEf': ['AFEJ'],
        'ALf': ['AALJ'],
        'MNf': ['AMNJ'],
    }
    to_crn = {'SUR_TEMP': ['TEMPG'], 'T_HR_AVG': ['TEMP2'], 'SOLARAD': ['RGRND'], 'SOIL_MOISTURE_5': ['SOIM1'], 'SOIL_MOISTURE_10': ['SOIM2']} #Need to look into.
    to_aeronet = {}
    to_cems = {}
    mapping_tables = {
        'improve': to_improve,
        'aqs': to_aqs,
        'airnow': to_airnow,
        'crn': to_crn,
        'cems': to_cems,
        'nadp': to_nadp,
        'aeronet': to_aeronet,
    }
    dset = dset.assign_attrs({'mapping_tables': mapping_tables})
    return dset


# Arrays for different gasses and pm groupings
accumulation = array(
    [
        'aalj',
        'aalk1j',
        'aalk2j',
        'abnz1j',
        'abnz2j',
        'abnz3j',
        'acaj',
        'aclj',
        'aecj',
        'afej',
        'aiso1j',
        'aiso2j',
        'aiso3j',
        'akj',
        'amgj',
        'amnj',
        'anaj',
        'anh4j',
        'ano3j',
        'aolgaj',
        'aolgbj',
        'aorgcj',
        'aothrj',
        'apah1j',
        'apah2j',
        'apah3j',
        'asij',
        'aso4j',
        'asqtj',
        'atij',
        'atol1j',
        'atol2j',
        'atol3j',
        'atrp1j',
        'atrp2j',
        'axyl1j',
        'axyl2j',
        'axyl3j',
    ] #Had to remove 'APNCOMJ','APOCJ', 'AORGAJ','AORGPAJ','AORGBJ',double check this list?
)
aitken = array(['acli', 'aeci', 'anai', 'anh4i', 'ano3i', 'aothri', 'aso4i']) #Had to remove 'APNCOMI','APOCI','AORGAI','AORGPAI', and 'AORGBI' double check this?
coarse = array(['aclk', 'acors', 'anh4k', 'ano3k', 'aseacat', 'aso4k', 'asoil'])
noy_gas = array(['no', 'no2', 'no3', 'n2o5', 'hono', 'hno3', 'pan', 'panx', 'pna', 'intr', 'ntr1', 'ntr2', 'cron', 'opan']) #RHS updated to include 3 NTR's. Need to compare and make sure have all of them. Removed CRN2, CRNO, and CRPX. Not in output are they in mech?
pec = array(['aeci', 'aecj'])
pso4 = array(['aso4i', 'aso4j'])
pno3 = array(['ano3i', 'ano3j'])
pnh4 = array(['anh4i', 'anh4j'])
pcl = array(['acli', 'aclj'])
poc = array(
    [
        'aothri',
        'atol1j',
        'atol2j',
        'atol3j',
        'atrp1j',
        'atrp2j',
        'axyl1j',
        'axyl2j',
        'axyl3j',
        'aolgaj',
        'aolgbj',
        'aorgcj',
        'aothrj',
        'apah1j',
        'apah2j',
        'apah3j',
        'asqtj',
        'aiso1j',
        'aiso2j',
        'aiso3j',
        'aalk1j',
        'aalk2j',
        'abnz1j',
        'abnz2j',
        'abnz3j',
    ] #Had to remove some. Double check this later? 'APNCOMI','APOCI','AORGAI','AORGPAI','AORGBI',
        #'AORGAJ','AORGPAJ','AORGBJ','APNCOMJ','APOCJ', 'AORGAI','AORGAJ','AORGPAI','AORGPAJ',
        #'AORGBI','AORGBJ',
)
minerals = array(['aalj', 'acaj', 'afej', 'akj', 'amgj', 'amnj', 'anaj', 'atij', 'asij'])

def _calc_hgt(f):
    """Calculates the geopotential height in m from the variables hgtsfc and delz.
    Note: To use this function the delz value needs to go from surface to TOA in vertical.
    Parameters
    ----------
    f : xarray.DataSet
        the NEMSIO opened object.  Can be lazily loaded.
    Returns
    -------
    xr.DataArray
        Geoptential height with varialbes, coordinates and variable attributes.
    """
    sfc = f.hgtsfc
    dz = f.delz
    z = dz + sfc
    z = z.rolling(z=len(f.z), min_periods=1).sum()
    z.name = 'geohgt'
    z.attrs['long_name'] = 'Geopotential Height'
    z.attrs['units'] = 'm'
    return z


def _calc_pressure(dset):
    """Calculate the pressure in pa form presss and ak and bk constants.
    
    Interface pressures are calculated by:
    phalf(k) = a(k) + surfpres * b(k)
    
    Mid layer pressures are calculated by:
    pfull(k) = (phalf(k+1)-phalf(k))/log(phalf(k+1)/phalf(k))
    
    Parameters
    ----------
    dset : xarray.Dataset
        nemsio dataset opened
    Returns
    -------
    xarray.DataArray
        Description of returned object.
    """
    pres = dset.dpres.copy().load() #Have to load into memory here so can assign levels.
    srfpres = dset.pressfc.copy()
    for k in range(len(dset.z)):
        pres_2 = dset.ak[k+1] + srfpres * dset.bk[k+1]
        pres_1 = dset.ak[k] + srfpres * dset.bk[k]
        pres[:,k,:,:] = (pres_2-pres_1)/np.log(pres_2/pres_1)
    
    pres.name = 'pres_pa'
    pres.attrs['units'] = 'pa'
    pres.attrs['long_name'] = 'Mid Layer Pressure'
    return pres