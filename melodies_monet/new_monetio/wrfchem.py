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
                   convert_to_ppb=True,
                   mech = 'racm_esrl_vcp',
                   var_list = ['o3'],
                   surf_only=False,
                   **kwargs):
    """Method to open WRF-chem and RAP-chem netcdf files.

    Parameters
    ----------
    fname : string or list
        fname is the path to the file or files.  It will accept hot keys in
        strings as well.
    convert_to_ppb : boolean
        If true the units of the gas species will be converted to ppbV
    mech: str
        Mechanism to be used for calculating sums. Supported mechanisms 
        include: 'racm_esrl_vcp' and 'redhc'
    var_list: list
        List of variables to include in output. MELODIES-MONET only reads in 
        variables need to plot in order to save on memory and simulation cost
        especially for vertical data
    surf_only: boolean
        Whether to save only surface data to save on memory and computational 
        cost (True) or not (False).

    Returns
    -------
    xarray.DataSet
        WRF-Chem or RAP-Chem model dataset in standard format for use 
        in MELODIES-MONET


    """
    #Get dictionary of summed species for the mechanism of choice.
    dict_sum = dict_species_sums(mech=mech)
    
    list_calc_sum = [] 
    for var_sum in ['noy_gas','noy_aer','nox','pm25_cl','pm25_ec','pm25_na',
                    'pm25_nh4','pm25_no3','pm25_so4','pm25_om']:
        if var_sum in var_list:
            var_list.extend(dict_sum[var_sum])
            var_list.remove(var_sum)
            list_calc_sum.append(var_sum)

    wrflist = []
    for files in fname:
        wrflist.append(Dataset(files))
    
    if surf_only == False:
        #Add some additional defaults needed for aircraft analysis
        #Turn this on also if need to convert aerosols 
        var_list.append('pres')
        var_list.append('height')
        var_list.append('tk')
        var_list.append('height_agl')
        var_list.append('PSFC')
        #need to calculate surface pressure and dp and optionally dz here. 
        #Meng or Jian since you need this for satellite info, I'll have you add these here.
   
    var_wrf_list = []
    for var in var_list:
        if var == 'pres': #Insert special versions.
            var_wrf = getvar(wrflist,var,timeidx=ALL_TIMES,method='cat',squeeze=False,units='Pa')
        elif var == 'height':
            var_wrf = getvar(wrflist,var,timeidx=ALL_TIMES,method='cat',squeeze=False,units='m')
        elif var == 'height_agl':
            var_wrf = getvar(wrflist,var,timeidx=ALL_TIMES,method='cat',squeeze=False,units='m') 
        else:
            var_wrf = getvar(wrflist,var,timeidx=ALL_TIMES,method='cat',squeeze=False)
        var_wrf_list.append(var_wrf)

    dset = xr.merge(var_wrf_list)

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
    
    #Convert names to standards used in MONET
    dset = dset.rename({'Time': 'time','south_north':'y',
                        'west_east':'x','XLONG':'longitude',
                        'XLAT':'latitude'})
        
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
    
    # add lazy diagnostic variables
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

    dset = dset.reset_index(['XTIME','datetime'],drop=True)
    if surf_only == False:
        #Reset more variables
        dset = dset.rename({'bottom_top': 'z','temp':'temperature_k',
                            'height':'alt_msl_m_mid','height_agl':'alt_agl_m_mid',
                            'PSFC': 'surfpres_pa','pressure': 'pres_pa_mid'}) #
        dset2 = dset
    else:
        #Expand into z coordinate so that format is consistent.
        dset2 = dset.expand_dims('z',axis=3).copy()
    
    dset2 = dset2.reset_coords()
    dset2 = dset2.set_coords(['latitude','longitude'])
    
    return dset2 

def _get_keys(d):
    """Calculates keys

    Parameters
    ----------
    d : xarray.Dataset
        WRF-Chem model data

    Returns
    -------
    list
        list of keys

    """
    keys = Series([i for i in d.data_vars.keys()])
    return keys

def add_lazy_noy_g(d,dict_sum):
    """Calculates NOy gas

    Parameters
    ----------
    d : xarray.Dataset
        WRF-Chem model data

    Returns
    -------
    xarray.Dataset
        WRF-Chem model data including new NOy gas calculation

    """
    keys = _get_keys(d)
    allvars = Series(dict_sum['noy_gas'])
    weights = Series(dict_sum['noy_gas_weight'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        newweights = weights.loc[index]
        d['noy_gas'] = add_multiple_lazy(d, newkeys, weights=newweights)
        d['noy_gas'] = d['noy_gas'].assign_attrs({'name': 'noy_gas', 'long_name': 'NOy gases'})
    return d                      

def add_lazy_noy_a(d,dict_sum):
    """Calculates NOy aerosol

    Parameters
    ----------
    d : xarray.Dataset
        WRF-Chem model data

    Returns
    -------
    xarray.Dataset
        WRF-Chem model data including new NOy aerosol calculation

    """
    keys = _get_keys(d)
    allvars = Series(dict_sum['noy_aer'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        d['noy_aer'] = add_multiple_lazy(d, newkeys)
        d['noy_aer'] = d['noy_aer'].assign_attrs({'units': '$\mu g m^{-3}$','name': 'noy_aer', 
                                              'long_name': 'NOy aerosol'})
    return d
                      
def add_lazy_nox(d,dict_sum):
    """Calculates NOx

    Parameters
    ----------
    d : xarray.Dataset
        WRF-Chem model data

    Returns
    -------
    xarray.Dataset
        WRF-Chem model data including new NOx calculation

    """
    keys = _get_keys(d)
    allvars = Series(dict_sum['nox'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        d['nox'] = add_multiple_lazy(d, newkeys)
        d['nox'] = d['nox'].assign_attrs({'name': 'nox', 'long_name': 'nox'})
    return d

def add_lazy_cl_pm25(d,dict_sum):
    """Calculates CL

    Parameters
    ----------
    d : xarray.Dataset
        WRF-Chem model data

    Returns
    -------
    xarray.Dataset
        WRF-Chem model data including new CL calculation

    """
    keys = _get_keys(d)
    allvars = Series(dict_sum['pm25_cl'])
    weights = Series(dict_sum['pm25_cl_weight'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['pm25_cl'] = add_multiple_lazy(d, newkeys, weights=neww)
        d['pm25_cl'] = d['pm25_cl'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'pm25_cl', 
                                                  'long_name': 'PM2.5 CL assuming coarse mode 20%'})
    return d

def add_lazy_ec_pm25(d,dict_sum):
    """Calculates EC

    Parameters
    ----------
    d : xarray.Dataset
        WRF-Chem model data

    Returns
    -------
    xarray.Dataset
        WRF-Chem model data including new EC calculation

    """
    keys = _get_keys(d)
    allvars = Series(dict_sum['pm25_ec'])
    weights = Series(dict_sum['pm25_ec_weight'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['pm25_ec'] = add_multiple_lazy(d, newkeys, weights=neww)
        d['pm25_ec'] = d['pm25_ec'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'pm25_ec', 
                                                  'long_name': 'PM2.5 EC assuming coarse mode 20%'})
    return d

def add_lazy_na_pm25(d,dict_sum):
    """Calculates NA

    Parameters
    ----------
    d : xarray.Dataset
        WRF-Chem model data

    Returns
    -------
    xarray.Dataset
        WRF-Chem model data including new NA calculation

    """
    keys = _get_keys(d)
    allvars = Series(dict_sum['pm25_na'])
    weights = Series(dict_sum['pm25_na_weight'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['pm25_na'] = add_multiple_lazy(d, newkeys, weights=neww)
        d['pm25_na'] = d['pm25_na'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'pm25_na', 
                                                  'long_name': 'PM2.5 NA assuming coarse mode 20%'})
    return d

def add_lazy_nh4_pm25(d,dict_sum):
    """Calculates NH4

    Parameters
    ----------
    d : xarray.Dataset
        WRF-Chem model data

    Returns
    -------
    xarray.Dataset
        WRF-Chem model data including new NH4 calculation

    """
    keys = _get_keys(d)
    allvars = Series(dict_sum['pm25_nh4'])
    weights = Series(dict_sum['pm25_nh4_weight'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['pm25_nh4'] = add_multiple_lazy(d, newkeys, weights=neww)
        d['pm25_nh4'] = d['pm25_nh4'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'pm25_nh4', 
                                                  'long_name': 'PM2.5 NH4 assuming coarse mode 20%'})
    return d

def add_lazy_no3_pm25(d,dict_sum):
    """Calculates NO3 particulate

    Parameters
    ----------
    d : xarray.Dataset
        WRF-Chem model data

    Returns
    -------
    xarray.Dataset
        WRF-Chem model data including new NO3 particulate calculation

    """
    keys = _get_keys(d)
    allvars = Series(dict_sum['pm25_no3'])
    weights = Series(dict_sum['pm25_no3_weight'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['pm25_no3'] = add_multiple_lazy(d, newkeys, weights=neww)
        d['pm25_no3'] = d['pm25_no3'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'pm25_no3', 
                                                  'long_name': 'PM2.5 NO3 assuming coarse mode 20%'})
    return d

def add_lazy_so4_pm25(d,dict_sum):
    """Calculates SO4

    Parameters
    ----------
    d : xarray.Dataset
        WRF-Chem model data

    Returns
    -------
    xarray.Dataset
        WRF-Chem model data including new SO4 calculation

    """
    keys = _get_keys(d)
    allvars = Series(dict_sum['pm25_so4'])
    weights = Series(dict_sum['pm25_so4_weight'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['pm25_so4'] = add_multiple_lazy(d, newkeys, weights=neww)
        d['pm25_so4'] = d['pm25_so4'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'pm25_so4', 
                                                  'long_name': 'PM2.5 SO4 assuming coarse mode 20%'})
    return d

def add_lazy_om_pm25(d,dict_sum):
    """Calculates OM

    Parameters
    ----------
    d : xarray.Dataset
        WRF-Chem model data

    Returns
    -------
    xarray.Dataset
        WRF-Chem model data including new OM calculation

    """
    keys = _get_keys(d)
    allvars = Series(dict_sum['pm25_om'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        d['pm25_om'] = add_multiple_lazy(d, newkeys)
        d['pm25_om'] = d['pm25_om'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'pm25_om', 
                                                  'long_name': 'PM2.5 OM'})
    return d
                       
def add_multiple_lazy(dset, variables, weights=None):
    """Sums variables

    Parameters
    ----------
    d : xarray.Dataset
        WRF-Chem model data
    variables : series
        series of variables
    variables : series
        series of weights to apply to each variable during the sum 

    Returns
    -------
    xarray.Dataarray
        Weighted sum of all specified variables

    """
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
#for now this is mechanism independent.
    """Predefined mapping tables for different observational parings used when
        combining data.

    Returns
    -------
    dictionary
        dictionary defining default mapping tables

    """
    to_airnow = {
        'OZONE': 'o3',
        'PM2.5': 'PM2_5_DRY',
        'PM10': 'PM10',
        'CO': 'co',
        'SO2': 'so2',
        'NO': 'no',
        'NO2': 'no2',
    }
    dset = dset.assign_attrs({'mapping_tables_to_airnow': to_airnow})
    return dset

def dict_species_sums(mech):
    """Predefined mapping tables for different observational parings used when
        combining data.
    
    Parameters
    ----------
    mech : string
        mechanism name
        
    Returns
    -------
    dictionary
        dictionary defining the variables to sum based on the specified mechanism

    """
    if mech == 'racm_esrl_vcp':
        sum_dict = {}
        # Arrays for different gasses and pm groupings
        sum_dict.update({'noy_gas' : ['hno3','no', 'no2', 'no3', 'pan', 'tpan', 'hono', 'hno4',
                         'onit', 'n2o5', 'ison', 'nald','mpan']}) 
        sum_dict.update({'noy_gas_weight' : [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1]})
        sum_dict.update({'noy_aer' : ['no3ai', 'no3aj']}) #Need to confirm here if there is a size cutoff for noy obs?
        sum_dict.update({'nox' : ['no', 'no2']})              
        sum_dict.update({'pm25_cl' : ['clai', 'claj']})
        sum_dict.update({'pm25_cl_weight' : [1, 1]})
        sum_dict.update({'pm25_ec' : ['eci', 'ecj']})
        sum_dict.update({'pm25_ec_weight' : [1, 1]})
        sum_dict.update({'pm25_na' : ['naai', 'naaj']})
        sum_dict.update({'pm25_na_weight' : [1, 1]})            
        sum_dict.update({'pm25_nh4' : ['nh4ai', 'nh4aj']})
        sum_dict.update({'pm25_nh4_weight' : [1, 1]})
        sum_dict.update({'pm25_no3' : ['no3ai', 'no3aj']})
        sum_dict.update({'pm25_no3_weight' : [1, 1]})
        sum_dict.update({'pm25_so4' : ['so4ai', 'so4aj']})
        sum_dict.update({'pm25_so4_weight' : [1, 1]})
        sum_dict.update({'pm25_om' : ['asoa1i', 'asoa1j', 'asoa2i', 'asoa2j','asoa3i', 'asoa3j',
                                     'asoa4i', 'asoa4j','bsoa1i', 'bsoa1j','bsoa2i', 'bsoa2j',
                                     'bsoa3i', 'bsoa3j','bsoa4i', 'bsoa4j','orgpai','orgpaj']}) 
    elif mech == 'redhc':
        sum_dict = {}
        # Arrays for different gasses and pm groupings
        sum_dict.update({'noy_gas' : ['hno3','no', 'no2', 'no3', 'pan', 'ho2no2',
                         'onit', 'n2o5']}) 
        sum_dict.update({'noy_gas_weight' : [1, 1, 1, 1, 1, 1, 1, 2]})
        sum_dict.update({'noy_aer' : ['no3ai', 'no3aj']}) #Need to confirm here if there is a size cutoff for noy obs?
        sum_dict.update({'nox' : ['no', 'no2']})              
        sum_dict.update({'pm25_cl' : ['clai', 'claj']})
        sum_dict.update({'pm25_cl_weight' : [1, 1]})
        sum_dict.update({'pm25_ec' : ['eci', 'ecj']})
        sum_dict.update({'pm25_ec_weight' : [1, 1]})
        sum_dict.update({'pm25_na' : ['naai', 'naaj']})
        sum_dict.update({'pm25_na_weight' : [1, 1]})            
        sum_dict.update({'pm25_nh4' : ['nh4ai', 'nh4aj']})
        sum_dict.update({'pm25_nh4_weight' : [1, 1]})
        sum_dict.update({'pm25_no3' : ['no3ai', 'no3aj']})
        sum_dict.update({'pm25_no3_weight' : [1, 1]})
        sum_dict.update({'pm25_so4' : ['so4ai', 'so4aj']})
        sum_dict.update({'pm25_so4_weight' : [1, 1]})
        sum_dict.update({'pm25_om' : ['asoa0j','asoa0i','asoa1i', 'asoa1j', 'asoa2i', 'asoa2j',
                                      'asoa3i', 'asoa3j','bsoa1i', 'bsoa1j','bsoa2i', 'bsoa2j',
                                      'bsoa3i', 'bsoa3j','poa0j', 'poa0i', 'poa1j', 'poa1i', 
                                      'poa2j', 'poa2i', 'poa3j', 'poa3i']})

    else:
        raise NotImplementedError('Mechanism not supported, update rrfs_cmaq.py file in MONETIO')
    
    return sum_dict
