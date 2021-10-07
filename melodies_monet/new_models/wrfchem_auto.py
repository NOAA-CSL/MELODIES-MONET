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
                   mech = 'racm_esrl_vcp',
                   var_list = ['o3'],
                   vert = True,
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
    
    if vert == True:
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
    
    #Calculate projection
    grid = ioapi_grid_from_dataset_wrf(dset, earth_radius=earth_radius)
    dset = dset.assign_attrs({'proj4_srs': grid})
    for i in dset.variables:
        dset[i] = dset[i].assign_attrs({'proj4_srs': grid})
        
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
    if vert == True:
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

def add_lazy_noy_g(d,dict_sum):
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
    keys = _get_keys(d)
    allvars = Series(dict_sum['nox'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        d['nox'] = add_multiple_lazy(d, newkeys)
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
        d['pm25_cl'] = add_multiple_lazy(d, newkeys, weights=neww)
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
        d['pm25_ec'] = add_multiple_lazy(d, newkeys, weights=neww)
        d['pm25_ec'] = d['pm25_ec'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'pm25_ec', 
                                                  'long_name': 'PM2.5 EC assuming coarse mode 20%'})
    return d

def add_lazy_na_pm25(d,dict_sum):
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
        A dictionary of to map to.

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
