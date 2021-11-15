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

def open_mfdataset(fname, 
                   convert_to_ppb=True, 
                   mech='cb6r3_ae6_aq',
                   var_list = None,
                   fname_pm25=None,
                   surf_only=False,
                   **kwargs):
    #Like WRF-chem add var list that just determines whether to calculate sums or not to speed this up.
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

    #Get dictionary of summed species for the mechanism of choice.
    dict_sum = dict_species_sums(mech=mech)
    
    if var_list is not None:
        #Read in only a subset of variables and only do calculations if needed.
        list_calc_sum = [] 
        for var_sum in ['PM25', 'PM25_wopc', 'PM10', 'noy_gas', 'noy_aer', 
                        'nox', 'pm25_cl', 'pm25_ec', 'pm25_ca', 'pm25_na',
                        'pm25_nh4', 'pm25_no3', 'pm25_so4', 'pm25_om']:
            if var_sum in var_list:
                if var_sum == 'PM25':
                    var_list.extend(dict_sum['aitken'])
                    var_list.extend(dict_sum['accumulation'])
                    var_list.extend(dict_sum['coarse'])
                elif var_sum == 'PM25_wopc':
                    var_list.extend(dict_sum['aitken'])
                    var_list.extend(dict_sum['accumulation_wopc'])
                    var_list.extend(dict_sum['coarse'])
                elif var_sum == 'PM10':
                    var_list.extend(dict_sum['aitken'])
                    var_list.extend(dict_sum['accumulation'])
                    var_list.extend(dict_sum['coarse'])
                else:
                    var_list.extend(dict_sum[var_sum])
                var_list.remove(var_sum)
                list_calc_sum.append(var_sum)
        #append the other needed species.
        var_list.append('lat')
        var_list.append('lon')
        var_list.append('phalf')
        var_list.append('tmp')
        var_list.append('pressfc')
        var_list.append('dpres')
        var_list.append('hgtsfc')
        var_list.append('delz')
        
        #Remove duplicates just in case:
        var_list = list( dict.fromkeys(var_list) )
        #If variables in pm25 files are included remove these as these are not in the main file
        #And will be added later.
        for pm25_var in ['PM25_TOT','PM25_TOT_NSOM','PM25_EC','PM25_NH4',
                         'PM25_NO3','PM25_SO4','PM25_OC','PM25_OM']:
            if pm25_var in var_list:
                var_list.remove(pm25_var)
    
        # open the dataset using xarray
        dset = xr.open_mfdataset(fname, concat_dim='time', **kwargs)[var_list]
    else:
        #Read in all variables and do all calculations.
        dset = xr.open_mfdataset(fname, concat_dim='time', **kwargs)
        list_calc_sum = ['PM25', 'PM25_wopc', 'PM10', 'noy_gas', 'noy_aer', 
                        'nox', 'pm25_cl', 'pm25_ec', 'pm25_ca', 'pm25_na',
                        'pm25_nh4', 'pm25_no3', 'pm25_so4', 'pm25_om']
    
    if fname_pm25 is not None:
        #Add the processed pm2.5 species.
        dset_pm25 = xr.open_mfdataset(fname_pm25, concat_dim='time', **kwargs)
        dset_pm25 = dset_pm25.drop(labels=['lat','lon','pfull']) #Drop duplicate variables so can merge. 
        #Slight differences in pfull value between the files, but I assume that these still represent the
        #same pressure levels from the model dynf* files.
        #Attributes are formated differently in pm25 file so remove attributes and use those from dynf* files.
        dset_pm25.attrs = {}
        dset = dset.merge(dset_pm25)
        
    #Standardize some variable names
    dset = dset.rename({'grid_yt': 'y','grid_xt': 'x','pfull': 'z',
                        'phalf': 'z_i', #Interface pressure levels
                        'lon': 'longitude','lat': 'latitude',
                        'tmp': 'temperature_k', #standard temperature (kelvin)
                        'pressfc': 'surfpres_pa','dpres': 'dp_pa', #Change names so standard surfpres_pa and dp_pa
                        'hgtsfc': 'surfalt_m','delz':'dz_m'}) #Optional, but when available include altitude info

    #Calculate pressure. This has to go before sorting because ak and bk 
    #are not sorted as they are in attributes
    dset['pres_pa_mid'] = _calc_pressure(dset)

    #Adjust pressure levels for all models such that the surface is first.
    dset = dset.sortby('z', ascending=False)
    dset = dset.sortby('z_i', ascending=False)
    
    #Note this altitude calcs needs to always go after resorting.
    #Altitude calculations are all optional, but for each model add values that are easy to calculate.
    dset['alt_msl_m_full'] = _calc_hgt(dset)
    dset['dz_m'] = dset['dz_m']*-1. #Change to positive values.
    
    #Set coordinates
    dset = dset.reset_index(['x','y','z','z_i'],drop=True) #For now drop z_i no variables use it. 
    dset['latitude'] = dset['latitude'].isel(time=0)
    dset['longitude'] = dset['longitude'].isel(time=0)
    dset = dset.reset_coords()
    dset = dset.set_coords(['latitude','longitude'])  
    
    #These sums and units are quite expensive and memory intensive, 
    #so add option to shrink dataset to just surface when needed
    if surf_only == True:
        dset=dset.isel(z=0).expand_dims('z',axis=1)
        
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
                dset[i] = dset[i]*dset['pres_pa_mid']/dset['temperature_k']/287.05535
                dset[i].attrs['units'] = '$\mu g m^{-3}$'
    
    # add lazy diagnostic variables
    # Note that because there are so many species to sum. Summing the aerosols is slowing down the code.
    if 'PM25' in list_calc_sum:
        dset = add_lazy_pm25(dset,dict_sum)
    if 'PM25_wopc' in list_calc_sum:
        dset = add_lazy_pm25_wopc(dset,dict_sum)
    if 'PM10' in list_calc_sum:    
        dset = add_lazy_pm10(dset,dict_sum)
    if 'noy_gas' in list_calc_sum:
        dset = add_lazy_noy_g(dset,dict_sum)
    if 'noy_aer' in list_calc_sum:
        dset = add_lazy_noy_a(dset,dict_sum)
    if 'nox' in list_calc_sum:
        dset = add_lazy_nox(dset,dict_sum)
    if 'pm25_cl' in list_calc_sum:
        dset = add_lazy_cl_pm25(dset,dict_sum)
    if 'pm25_ec' in list_calc_sum:
        dset = add_lazy_ec_pm25(dset,dict_sum)
    if 'pm25_ca' in list_calc_sum:
        dset = add_lazy_ca_pm25(dset,dict_sum)
    if 'pm25_na' in list_calc_sum:
        dset = add_lazy_na_pm25(dset,dict_sum)
    if 'pm25_nh4' in list_calc_sum:
        dset = add_lazy_nh4_pm25(dset,dict_sum)
    if 'pm25_no3' in list_calc_sum:
        dset = add_lazy_no3_pm25(dset,dict_sum)
    if 'pm25_so4' in list_calc_sum:
        dset = add_lazy_so4_pm25(dset,dict_sum)
    if 'pm25_om' in list_calc_sum:
        dset = add_lazy_om_pm25(dset,dict_sum)
    # Change the times to pandas format
    dset['time'] = dset.indexes['time'].to_datetimeindex(unsafe=True)
    #Turn off warning for now. This is just because the model is in julian time

    return dset

def _get_keys(d):
    keys = Series([i for i in d.data_vars.keys()])
    return keys

def add_lazy_pm25(d,dict_sum):
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
    allvars = Series(concatenate([dict_sum['aitken'], dict_sum['accumulation'], dict_sum['coarse']]))
    weights = Series(concatenate([np.ones(len(dict_sum['aitken'])),
                                  np.ones(len(dict_sum['accumulation'])),
                                  np.full(len(dict_sum['coarse']),0.2)]))
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        newweights = weights.loc[index]
        d['PM25'] = add_multiple_lazy2(d, newkeys, weights=newweights)
        d['PM25'] = d['PM25'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'PM2.5', 
                                            'long_name': 'PM2.5 calculated by MONET assuming coarse mode 20%'})
    return d

def add_lazy_pm25_wopc(d,dict_sum):
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
    allvars = Series(concatenate([dict_sum['aitken'], dict_sum['accumulation_wopc'], dict_sum['coarse']]))
    weights = Series(concatenate([np.ones(len(dict_sum['aitken'])),
                                  np.ones(len(dict_sum['accumulation_wopc'])),
                                  np.full(len(dict_sum['coarse']),0.2)]))
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        newweights = weights.loc[index]
        d['PM25_wopc'] = add_multiple_lazy2(d, newkeys, weights=newweights)
        d['PM25_wopc'] = d['PM25_wopc'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'PM2.5_wopc', 
                                            'long_name': 'PM2.5 calculated by MONET assuming coarse mode 20% excluding apcsoj'})
    return d

def add_lazy_pm10(d,dict_sum):
    keys = _get_keys(d)
    allvars = Series(concatenate([dict_sum['aitken'], dict_sum['accumulation'], dict_sum['coarse']]))
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        d['PM10'] = add_multiple_lazy2(d, newkeys)
        d['PM10'] = d['PM10'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'PM10', 
                                            'long_name': 'Particulate Matter < 10 microns'})
    return d
                      
def add_lazy_noy_g(d,dict_sum):
    keys = _get_keys(d)
    allvars = Series(dict_sum['noy_gas'])
    weights = Series(dict_sum['noy_gas_weight'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        newweights = weights.loc[index]
        d['noy_gas'] = add_multiple_lazy2(d, newkeys, weights=newweights)
        d['noy_gas'] = d['noy_gas'].assign_attrs({'name': 'noy_gas', 'long_name': 'NOy gases'})
    return d                      

def add_lazy_noy_a(d,dict_sum):
    keys = _get_keys(d)
    allvars = Series(dict_sum['noy_aer'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        d['noy_aer'] = add_multiple_lazy2(d, newkeys)
        d['noy_aer'] = d['noy_aer'].assign_attrs({'units': '$\mu g m^{-3}$','name': 'noy_aer', 
                                              'long_name': 'NOy aerosol'})
    return d
                      
def add_lazy_nox(d,dict_sum):
    keys = _get_keys(d)
    allvars = Series(dict_sum['nox'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        d['nox'] = add_multiple_lazy2(d, newkeys)
        d['nox'] = d['nox'].assign_attrs({'name': 'nox', 'long_name': 'nox'})
    return d

def add_lazy_cl_pm25(d,dict_sum):
    keys = _get_keys(d)
    allvars = Series(dict_sum['pm25_cl'])
    weights = Series(dict_sum['pm25_cl_weight'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['pm25_cl'] = add_multiple_lazy2(d, newkeys, weights=neww)
        d['pm25_cl'] = d['pm25_cl'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'pm25_cl', 
                                                  'long_name': 'PM2.5 CL assuming coarse mode 20%'})
    return d

def add_lazy_ec_pm25(d,dict_sum):
    keys = _get_keys(d)
    allvars = Series(dict_sum['pm25_ec'])
    weights = Series(dict_sum['pm25_ec_weight'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['pm25_ec'] = add_multiple_lazy2(d, newkeys, weights=neww)
        d['pm25_ec'] = d['pm25_ec'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'pm25_ec', 
                                                  'long_name': 'PM2.5 EC assuming coarse mode 20%'})
    return d

def add_lazy_ca_pm25(d,dict_sum):
    keys = _get_keys(d)
    allvars = Series(dict_sum['pm25_ca'])
    weights = Series(dict_sum['pm25_ca_weight'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['pm25_ca'] = add_multiple_lazy2(d, newkeys, weights=neww)
        d['pm25_ca'] = d['pm25_ca'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'pm25_ca', 
                                                  'long_name': 'PM2.5 CA assuming coarse mode 20%'})
    return d


def add_lazy_na_pm25(d,dict_sum):
    keys = _get_keys(d)
    allvars = Series(dict_sum['pm25_na'])
    weights = Series(dict_sum['pm25_na_weight'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['pm25_na'] = add_multiple_lazy2(d, newkeys, weights=neww)
        d['pm25_na'] = d['pm25_na'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'pm25_na', 
                                                  'long_name': 'PM2.5 NA assuming coarse mode 20%'})
    return d

def add_lazy_nh4_pm25(d,dict_sum):
    keys = _get_keys(d)
    allvars = Series(dict_sum['pm25_nh4'])
    weights = Series(dict_sum['pm25_nh4_weight'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['pm25_nh4'] = add_multiple_lazy2(d, newkeys, weights=neww)
        d['pm25_nh4'] = d['pm25_nh4'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'pm25_nh4', 
                                                  'long_name': 'PM2.5 NH4 assuming coarse mode 20%'})
    return d

def add_lazy_no3_pm25(d,dict_sum):
    keys = _get_keys(d)
    allvars = Series(dict_sum['pm25_no3'])
    weights = Series(dict_sum['pm25_no3_weight'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['pm25_no3'] = add_multiple_lazy2(d, newkeys, weights=neww)
        d['pm25_no3'] = d['pm25_no3'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'pm25_no3', 
                                                  'long_name': 'PM2.5 NO3 assuming coarse mode 20%'})
    return d

def add_lazy_so4_pm25(d,dict_sum):
    keys = _get_keys(d)
    allvars = Series(dict_sum['pm25_so4'])
    weights = Series(dict_sum['pm25_so4_weight'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['pm25_so4'] = add_multiple_lazy2(d, newkeys, weights=neww)
        d['pm25_so4'] = d['pm25_so4'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'pm25_so4', 
                                                  'long_name': 'PM2.5 SO4 assuming coarse mode 20%'})
    return d

def add_lazy_om_pm25(d,dict_sum):
    keys = _get_keys(d)
    allvars = Series(dict_sum['pm25_om'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        d['pm25_om'] = add_multiple_lazy2(d, newkeys)
        d['pm25_om'] = d['pm25_om'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'pm25_om', 
                                                  'long_name': 'PM2.5 OM'})
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

def add_multiple_lazy2(dset, variables, weights=None):
    
    dset2 = dset[variables.values]
    if weights is not None:
        for i, j in zip(variables.values, weights.values):
            dset2[i] = dset2[i] * j
    
    new = dset2.to_array().sum('variable')
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
        'PM2.5': ['PM25'],
        'CO': ['co'],
        'NOY': ['NOy'],
        'NOX': ['NOx'],
        'SO2': ['so2'],
        'NO': ['no'],
        'NO2': ['no2'],
    }
    to_airnow = {
        'OZONE': ['o3'],
        'PM2.5': ['PM25'],
        'CO': ['co'],
        'NOY': ['NOy'],
        'NOX': ['NOx'],
        'SO2': ['so2'],
        'NO': ['no'],
        'NO2': ['no2'],
    }
    to_crn = {} 
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

#For the different mechanisms, just update these arrays as needed. 

def dict_species_sums(mech):
    if mech == 'cb6r3_ae6_aq':
        sum_dict = {}
        # Arrays for different gasses and pm groupings
        sum_dict.update({'accumulation': ['aso4j', 'ano3j', 'anh4j', 'anaj', 'aclj', 'aecj', 'aothrj',
                                                'afej', 'asij', 'atij', 'acaj', 'amgj', 'amnj', 'aalj', 
                                                'akj', 'alvpo1j', 'asvpo1j', 'asvpo2j', 'asvpo3j', 'aivpo1j',
                                                'axyl1j', 'axyl2j', 'axyl3j', 'atol1j', 'atol2j', 'atol3j',
                                                'abnz1j', 'abnz2j', 'abnz3j', 'aiso1j', 'aiso2j', 'aiso3j',
                                                'atrp1j', 'atrp2j', 'asqtj', 'aalk1j', 'aalk2j', 'apah1j',
                                                'apah2j', 'apah3j', 'aorgcj', 'aolgbj', 'aolgaj', 'alvoo1j',
                                                'alvoo2j','asvoo1j','asvoo2j','asvoo3j','apcsoj']})
        sum_dict.update({'accumulation_wopc': ['aso4j', 'ano3j', 'anh4j', 'anaj', 'aclj', 'aecj', 'aothrj',
                                                'afej', 'asij', 'atij', 'acaj', 'amgj', 'amnj', 'aalj', 
                                                'akj', 'alvpo1j', 'asvpo1j', 'asvpo2j', 'asvpo3j', 'aivpo1j',
                                                'axyl1j', 'axyl2j', 'axyl3j', 'atol1j', 'atol2j', 'atol3j',
                                                'abnz1j', 'abnz2j', 'abnz3j', 'aiso1j', 'aiso2j', 'aiso3j',
                                                'atrp1j', 'atrp2j', 'asqtj', 'aalk1j', 'aalk2j', 'apah1j',
                                                'apah2j', 'apah3j', 'aorgcj', 'aolgbj', 'aolgaj', 'alvoo1j',
                                                'alvoo2j','asvoo1j','asvoo2j','asvoo3j']})
        sum_dict.update({'aitken' : ['aso4i','ano3i','anh4i','anai','acli', 'aeci','aothri',
                                          'alvpo1i','asvpo1i','asvpo2i','alvoo1i','alvoo2i','asvoo1i','asvoo2i']})
        sum_dict.update({'coarse' : ['asoil','acors', 'aseacat', 'aclk', 'aso4k', 'ano3k', 'anh4k']})
        sum_dict.update({'noy_gas' : ['no', 'no2', 'no3', 'n2o5', 'hono', 'hno3', 'pna',
                         'cron', 'clno2', 'pan', 'panx','opan', 'ntr1', 'ntr2','intr']}) 
        sum_dict.update({'noy_gas_weight' : [1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]})
        sum_dict.update({'noy_aer' : ['ano3i', 'ano3j', 'ano3k']}) #Need to confirm here if there is a size cutoff for noy obs?
        sum_dict.update({'nox' : ['no', 'no2']})              
        sum_dict.update({'pm25_cl' : ['acli', 'aclj','aclk']})
        sum_dict.update({'pm25_cl_weight' : [1, 1, 0.2]})
        sum_dict.update({'pm25_ec' : ['aeci', 'aecj']})
        sum_dict.update({'pm25_ec_weight' : [1, 1]})
        sum_dict.update({'pm25_na' : ['anai', 'anaj','aseacat', 'asoil', 'acors']})
        sum_dict.update({'pm25_na_weight' : [1, 1, 0.2 * 0.8373, 0.2 * 0.0626, 0.2 * 0.0023]})
        sum_dict.update({'pm25_ca' : ['acaj','aseacat', 'asoil', 'acors']})
        sum_dict.update({'pm25_ca_weight' : [1, 0.2 * 0.0320, 0.2 * 0.0838, 0.2 * 0.0562]})              
        sum_dict.update({'pm25_nh4' : ['anh4i', 'anh4j','anh4k']})
        sum_dict.update({'pm25_nh4_weight' : [1, 1, 0.2]})
        sum_dict.update({'pm25_no3' : ['ano3i', 'ano3j','ano3k']})
        sum_dict.update({'pm25_no3_weight' : [1, 1, 0.2]})
        sum_dict.update({'pm25_so4' : ['aso4i', 'aso4j','aso4k']})
        sum_dict.update({'pm25_so4_weight' : [1, 1, 0.2]})
        sum_dict.update({'pm25_om' : ['alvpo1i','asvpo1i','asvpo2i','alvoo1i','alvoo2i','asvoo1i','asvoo2i',
                                            'alvpo1j', 'asvpo1j', 'asvpo2j', 'asvpo3j', 'aivpo1j',
                                            'axyl1j', 'axyl2j', 'axyl3j', 'atol1j', 'atol2j', 'atol3j',
                                            'abnz1j', 'abnz2j', 'abnz3j', 'aiso1j', 'aiso2j', 'aiso3j',
                                            'atrp1j', 'atrp2j', 'asqtj', 'aalk1j', 'aalk2j', 'apah1j',
                                            'apah2j', 'apah3j', 'aorgcj', 'aolgbj', 'aolgaj', 'alvoo1j',
                                            'alvoo2j','asvoo1j','asvoo2j','asvoo3j','apcsoj']})              

    else:
        raise NotImplementedError('Mechanism not supported, update rrfs_cmaq.py file in MONETIO')
    
    return sum_dict
        

def _calc_hgt(f):
    """Calculates the geopotential height in m from the variables hgtsfc and delz.
    Note: To use this function the delz value needs to go from surface to TOA in vertical.
    Because we are adding the height of each grid box these are really grid top values
    Parameters
    ----------
    f : xarray.DataSet
        the NEMSIO opened object.  Can be lazily loaded.
    Returns
    -------
    xr.DataArray
        Geoptential height with varialbes, coordinates and variable attributes.
    """
    sfc = f.surfalt_m.load()
    dz = f.dz_m.load()*-1. #These are negative in RRFS-CMAQ, but you resorted and are adding from the surface, so make them positive.
    dz[:,0,:,:] = dz[:,0,:,:] + sfc #Add the surface altitude to the first model level only
    z = dz.rolling(z=len(f.z), min_periods=1).sum()
    z.name = 'alt_msl_m_full'
    z.attrs['long_name'] = 'Altitude MSL Full Layer in Meters'
    z.attrs['units'] = 'm'
    return z


def _calc_pressure(dset):
    """Calculate the pressure in pa from presss and ak and bk constants.
    
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
    pres = dset.dp_pa.copy().load() #Have to load into memory here so can assign levels.
    srfpres = dset.surfpres_pa.copy().load()
    for k in range(len(dset.z)):
        pres_2 = dset.ak[k+1] + srfpres * dset.bk[k+1]
        pres_1 = dset.ak[k] + srfpres * dset.bk[k]
        pres[:,k,:,:] = (pres_2-pres_1)/np.log(pres_2/pres_1)
    
    pres.name = 'pres_pa_mid'
    pres.attrs['units'] = 'pa'
    pres.attrs['long_name'] = 'Pressure Mid Layer in Pa'
    return pres
