""" CMAQ File Reader """
import xarray as xr
from numpy import array, concatenate
from pandas import Series, to_datetime

from monetio import grids
#from ..grids import get_ioapi_pyresample_area_def, grid_from_dataset


def can_do(index):
    if index.max():
        return True
    else:
        return False

def open_mfdataset(fname, 
                   earth_radius=6370000, 
                   convert_to_ppb=True, 
                   drop_duplicates=False, 
                   mech='cb6r3_ae6_aq',
                   var_list = ['O3'],
                   fname_vert=None,
                   fname_surf=None,
                   concatenate_forecasts=False,
                   **kwargs):
    """Method to open CMAQ IOAPI netcdf files.

    Parameters
    ----------
    fname : string or list
        fname is the path to the file or files.  It will accept hot keys in
        strings as well.
    earth_radius : float
        The earth radius used for the map projection
    convert_to_ppb : boolean
        If true the units of the gas species will be converted to ppbV
    mech: str
        Mechanism to be used for calculating sums. Only one mechanism is 
        included now. But plan to adapt this like the rrfs_cmaq.py code to have 
        default mechanisms. 
    var_list: list
        List of variables to include in output. MELODIES-MONET only reads in 
        variables need to plot in order to save on memory and simulation cost
        especially for vertical data
    fname_vert: string
        Path to the metcro3d file in CMAQ. This file is needed to calculate 
        vertical metrics of grid.
    fname_surf: string
        Path to the metcro2d file in CMAQ. This file is needed to calculate 
        vertical metrics of grid.
    concatentate_forecasts: boolean
        Whether you need to concatenate_forecasts (True) or not (False). Time 
        in the default CMAQ model data are not in a format to facilitate 
        concatenation by xarray.

    Returns
    -------
    xarray.DataSet
        CMAQ model dataset in standard format for use in MELODIES-MONET


    """
    #For CMAQ, times are not in file to concatenate on for different forecast periods,
    #so need to read in one file at a time, calc times, and then merge together.
    if concatenate_forecasts == True:
        dset_list = []
        for file in fname:
            # open the dataset using xarray
            dset = xr.open_mfdataset(file, concat_dim='TSTEP', **kwargs)
            # get the times
            dset = _get_times(dset, drop_duplicates=drop_duplicates)
            dset_list.append(dset)
        dset = xr.concat(dset_list,'time')  
    else:
        # open the dataset using xarray
        dset = xr.open_mfdataset(fname, concat_dim='TSTEP', **kwargs)
        # get the times
        dset = _get_times(dset, drop_duplicates=drop_duplicates)
        
    # add lazy diagnostic variables
    if 'PM25' in var_list:
        dset = add_lazy_pm25(dset)
    if 'PM10' in var_list:    
        dset = add_lazy_pm10(dset)
    if 'PM_COURSE' in var_list:
        dset = add_lazy_pm_course(dset)
    if 'CLf' in var_list:
        dset = add_lazy_clf(dset)
    if 'NAf' in var_list:
        dset = add_lazy_naf(dset)
    if 'CAf' in var_list:
        dset = add_lazy_caf(dset)
    if 'NOy' in var_list:
        dset = add_lazy_noy(dset)
    if 'NOx' in var_list:
        dset = add_lazy_nox(dset)
    if 'NO3f' in var_list:
        dset = add_lazy_no3f(dset)
    if 'NH4f' in var_list:
        dset = add_lazy_nh4f(dset)
    if 'SO4f' in var_list:
        dset = add_lazy_so4f(dset)
    if 'RH' in var_list:
        dset = add_lazy_rh(dset)

    # get the grid information
    grid = grids.grid_from_dataset(dset, earth_radius=earth_radius)
    area_def = grids.get_ioapi_pyresample_area_def(dset, grid)
    # assign attributes for dataset and all DataArrays
    dset = dset.assign_attrs({'proj4_srs': grid})
    for i in dset.variables:
        dset[i] = dset[i].assign_attrs({'proj4_srs': grid})
        for j in dset[i].attrs:
            dset[i].attrs[j] = dset[i].attrs[j].strip()
        # dset[i] = dset[i].assign_attrs({'area': area_def})
    # dset = dset.assign_attrs(area=area_def)

    # get the lat lon
    dset = _get_latlon(dset, area_def)
    
    #Remove TFLAG so can merge and make sure dimensions consistent with other models.
    dset = dset.drop(labels='TFLAG')
    
    if fname_vert is not None:
        dset_vert = xr.open_mfdataset(fname_vert, concat_dim='TSTEP', **kwargs)
        dset_vert = _get_times(dset_vert, drop_duplicates=drop_duplicates)
        dset_vert = dset_vert.drop(labels='TFLAG')
        dset = dset.merge(dset_vert)
        #Rename variables 
        dset = dset.rename({'TA':'temperature_k',
                            'PRES':'pres_pa_mid',
                            'PRES-F_lvl':'pres_pa_full',
                            'ZH':'alt_agl_m_mid',
                            'ZF':'alt_agl_m_full'})
    if fname_surf is not None:
        dset_surf = xr.open_mfdataset(fname_surf, concat_dim='TSTEP', **kwargs)
        dset_surf = _get_times(dset_surf, drop_duplicates=drop_duplicates)
        dset_surf = dset_surf.drop(labels='TFLAG').squeeze()
        dset = dset.merge(dset_surf)
        dset = dset.rename({'PRSFC':'surfpres_pa'})
        
    # get Predefined mapping tables for observations
    # d set = _predefined_mapping_tables(dset)
    # rename dimensions
    dset = dset.rename({'COL': 'x', 'ROW': 'y', 'LAY': 'z'})
    
    # convert all gas species to ppbv
    if convert_to_ppb:
        for i in dset.variables:
            if 'units' in dset[i].attrs:
                if 'ppmV' in dset[i].attrs['units']:
                    dset[i] *= 1000.0
                    dset[i].attrs['units'] = 'ppbV'

    # convert 'micrograms to \mu g'
    for i in dset.variables:
        if 'units' in dset[i].attrs:
            if 'micrograms' in dset[i].attrs['units']:
                dset[i].attrs['units'] = '$\mu g m^{-3}$'

    return dset


def _get_times(d, drop_duplicates):
    """Convert time from the default CMAQ model results to standard format

    Parameters
    ----------
    d : xarray.Dataset
        CMAQ model data
    drop_duplicates: boolean
        Whether to drop duplicates (True) or not (False)

    Returns
    -------
    xarray.Dataset
        CMAQ model data with time in standard format.

    """
    idims = len(d.TFLAG.dims)
    if idims == 2:
        tflag1 = Series(d['TFLAG'][:, 0]).astype(str).str.zfill(7)
        tflag2 = Series(d['TFLAG'][:, 1]).astype(str).str.zfill(6)
    else:
        tflag1 = Series(d['TFLAG'][:, 0, 0]).astype(str).str.zfill(7)
        tflag2 = Series(d['TFLAG'][:, 0, 1]).astype(str).str.zfill(6)
    date = to_datetime([i + j for i, j in zip(tflag1, tflag2)], format='%Y%j%H%M%S')
    if drop_duplicates:
        indexdates = Series(date).drop_duplicates(keep='last').index.values
        d = d.isel(TSTEP=indexdates)
        d['TSTEP'] = date[indexdates]
    else:
        d['TSTEP'] = date
    return d.rename({'TSTEP': 'time'})


def _get_latlon(dset, area):
    """Calculates the lat and lons from the pyreample.geometry.AreaDefinition

    Parameters
    ----------
    dset : xarray.Dataset
        CMAQ model data

    Returns
    -------
    xarray.Dataset
        CMAQ model data including the latitude and longitude in standard 
        format.

    """
    lon, lat = area.get_lonlats()
    dset['longitude'] = xr.DataArray(lon[::-1, :], dims=['ROW', 'COL'])
    dset['latitude'] = xr.DataArray(lat[::-1, :], dims=['ROW', 'COL'])
    dset = dset.assign_coords(longitude=dset.longitude, latitude=dset.latitude)
    return dset


def _get_keys(d):
    """Calculates keys

    Parameters
    ----------
    d : xarray.Dataset
        CMAQ model data

    Returns
    -------
    list
        list of keys

    """
    keys = Series([i for i in d.data_vars.keys()])
    return keys


def add_lazy_pm25(d):
    """Calculates PM2.5 sum. 20% of coarse mode is included in PM2.5 sum.

    Parameters
    ----------
    d : xarray.Dataset
        CMAQ model data

    Returns
    -------
    xarray.Dataset
        CMAQ model data including new PM2.5 calculation

    """
    keys = _get_keys(d)
    allvars = Series(concatenate([aitken, accumulation, coarse]))
    weights = Series(
        [
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
        ]
    )
    if 'PM25_TOT' in keys.to_list():
        d['PM25'] = d['PM25_TOT']
    else:
        index = allvars.isin(keys)
        if can_do(index):
            newkeys = allvars.loc[index]
            newweights = weights.loc[index]
            d['PM25'] = add_multiple_lazy(d, newkeys, weights=newweights)
            d['PM25'] = d['PM25'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'PM2.5', 'long_name': 'PM2.5'})
    return d


def add_lazy_pm10(d):
    """Calculates PM10 sum.

    Parameters
    ----------
    d : xarray.Dataset
        CMAQ model data

    Returns
    -------
    xarray.Dataset
        CMAQ model data including new PM10 calculation

    """
    keys = _get_keys(d)
    allvars = Series(concatenate([aitken, accumulation, coarse]))
    if 'PMC_TOT' in keys.to_list():
        d['PM10'] = d['PMC_TOT']
    else:
        index = allvars.isin(keys)
        if can_do(index):
            newkeys = allvars.loc[index]
            d['PM10'] = add_multiple_lazy(d, newkeys)
            d['PM10'] = d['PM10'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'PM10', 'long_name': 'Particulate Matter < 10 microns'})
    return d


def add_lazy_pm_course(d):
    """Calculates sum of particulate coarse mode.

    Parameters
    ----------
    d : xarray.Dataset
        CMAQ model data

    Returns
    -------
    xarray.Dataset
        CMAQ model data including new PM_COURSE calculation

    """
    keys = _get_keys(d)
    allvars = Series(coarse)
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        d['PM_COURSE'] = add_multiple_lazy(d, newkeys)
        d['PM_COURSE'] = d['PM_COURSE'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'PM_COURSE', 'long_name': 'Course Mode Particulate Matter'})
    return d


def add_lazy_clf(d):
    """Calculates sum of particulate Cl.

    Parameters
    ----------
    d : xarray.Dataset
        CMAQ model data

    Returns
    -------
    xarray.Dataset
        CMAQ model data including new CLf calculation

    """
    keys = _get_keys(d)
    allvars = Series(['ACLI', 'ACLJ', 'ACLK'])
    weights = Series([1, 1, 0.2])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['CLf'] = add_multiple_lazy(d, newkeys, weights=neww)
        d['CLf'] = d['CLf'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'CLf', 'long_name': 'Fine Mode particulate Cl'})
    return d


def add_lazy_caf(d):
    """Calculates sum of particulate Ca.

    Parameters
    ----------
    d : xarray.Dataset
        CMAQ model data

    Returns
    -------
    xarray.Dataset
        CMAQ model data including new CAf calculation

    """
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
    """Calculates sum of particulate Na.

    Parameters
    ----------
    d : xarray.Dataset
        CMAQ model data

    Returns
    -------
    xarray.Dataset
        CMAQ model data including new NAf calculation

    """
    keys = _get_keys(d)
    allvars = Series(['ANAI', 'ANAJ', 'ASEACAT', 'ASOIL', 'ACORS'])
    weights = Series([1, 1, 0.2 * 837.3 / 1000.0, 0.2 * 62.6 / 1000.0, 0.2 * 2.3 / 1000.0])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['NAf'] = add_multiple_lazy(d, newkeys, weights=neww)
        d['NAf'] = d['NAf'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'NAf', 'long_name': 'NAf'})
    return d


def add_lazy_so4f(d):
    """Calculates sum of particulate SO4.

    Parameters
    ----------
    d : xarray.Dataset
        CMAQ model data

    Returns
    -------
    xarray.Dataset
        CMAQ model data including new SO4f calculation

    """
    keys = _get_keys(d)
    allvars = Series(['ASO4I', 'ASO4J', 'ASO4K'])
    weights = Series([1.0, 1.0, 0.2])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['SO4f'] = add_multiple_lazy(d, newkeys, weights=neww)
        d['SO4f'] = d['SO4f'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'SO4f', 'long_name': 'SO4f'})
    return d


def add_lazy_nh4f(d):
    """Calculates sum of particulate NH4.

    Parameters
    ----------
    d : xarray.Dataset
        CMAQ model data

    Returns
    -------
    xarray.Dataset
        CMAQ model data including new NH4f calculation

    """
    keys = _get_keys(d)
    allvars = Series(['ANH4I', 'ANH4J', 'ANH4K'])
    weights = Series([1.0, 1.0, 0.2])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['NH4f'] = add_multiple_lazy(d, newkeys, weights=neww)
        d['NH4f'] = d['NH4f'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'NH4f', 'long_name': 'NH4f'})
    return d


def add_lazy_no3f(d):
    """Calculates sum of particulate NO3.

    Parameters
    ----------
    d : xarray.Dataset
        CMAQ model data

    Returns
    -------
    xarray.Dataset
        CMAQ model data including new NO3f calculation

    """
    keys = _get_keys(d)
    allvars = Series(['ANO3I', 'ANO3J', 'ANO3K'])
    weights = Series([1.0, 1.0, 0.2])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        neww = weights.loc[index]
        d['NO3f'] = add_multiple_lazy(d, newkeys, weights=neww)
        d['NO3f'] = d['NO3f'].assign_attrs({'units': '$\mu g m^{-3}$', 'name': 'NO3f', 'long_name': 'NO3f'})
    return d


def add_lazy_noy(d):
    """Calculates sum of NOy.

    Parameters
    ----------
    d : xarray.Dataset
        CMAQ model data

    Returns
    -------
    xarray.Dataset
        CMAQ model data including new NOy calculation

    """
    keys = _get_keys(d)
    allvars = Series(noy_gas)
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        d['NOy'] = add_multiple_lazy(d, newkeys)
        d['NOy'] = d['NOy'].assign_attrs({'name': 'NOy', 'long_name': 'NOy'})
    return d


def add_lazy_rh(d):
    # keys = Series([i for i in d.variables])
    # allvars = Series(['TEMP', 'Q', 'PRES'])
    # index = allvars.isin(keys)
    # if can_do(index):
    #     import atmos
    #     data = {
    #         'T': dset['TEMP'][:].compute().values,
    #         'rv': dset['Q'][:].compute().values,
    #         'p': dset['PRES'][:].compute().values
    #     }
    #     d['NOx'] = add_multiple_lazy(d, newkeys)
    #     d['NOx'] = d['NOx'].assign_attrs({'name': 'NOx', 'long_name': 'NOx'})
    return d


def add_lazy_nox(d):
    """Calculates sum of NOx.

    Parameters
    ----------
    d : xarray.Dataset
        CMAQ model data

    Returns
    -------
    xarray.Dataset
        CMAQ model data including new NOx calculation

    """
    keys = _get_keys(d)
    allvars = Series(['NO', 'NOX'])
    index = allvars.isin(keys)
    if can_do(index):
        newkeys = allvars.loc[index]
        d['NOx'] = add_multiple_lazy(d, newkeys)
        d['NOx'] = d['NOx'].assign_attrs({'name': 'NOx', 'long_name': 'NOx'})
    return d


def add_multiple_lazy(dset, variables, weights=None):
    """Sums variables

    Parameters
    ----------
    d : xarray.Dataset
        CMAQ model data
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
    """Predefined mapping tables for different observational parings used when
        combining data.

    Returns
    -------
    dictionary
        dictionary defining default mapping tables

    """
    to_improve = {}
    to_nadp = {}
    to_aqs = {
        'OZONE': ['O3'],
        'PM2.5': ['PM25'],
        'CO': ['CO'],
        'NOY': ['NOy'],
        'NOX': ['NOx'],
        'SO2': ['SO2'],
        'NOX': ['NOx'],
        'NO': ['NO'],
        'NO2': ['NO2'],
        'SO4f': ['SO4f'],
        'PM10': ['PM10'],
        'NO3f': ['NO3f'],
        'ECf': ['ECf'],
        'OCf': ['OCf'],
        'ETHANE': ['ETHA'],
        'BENZENE': ['BENZENE'],
        'TOLUENE': ['TOL'],
        'ISOPRENE': ['ISOP'],
        'O-XYLENE': ['XYL'],
        'WS': ['WSPD10'],
        'TEMP': ['TEMP2'],
        'WD': ['WDIR10'],
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
    to_airnow = {
        'OZONE': ['O3'],
        'PM2.5': ['PM25'],
        'CO': ['CO'],
        'NOY': ['NOy'],
        'NOX': ['NOx'],
        'SO2': ['SO2'],
        'NOX': ['NOx'],
        'NO': ['NO'],
        'NO2': ['NO2'],
        'SO4f': ['SO4f'],
        'PM10': ['PM10'],
        'NO3f': ['NO3f'],
        'ECf': ['ECf'],
        'OCf': ['OCf'],
        'ETHANE': ['ETHA'],
        'BENZENE': ['BENZENE'],
        'TOLUENE': ['TOL'],
        'ISOPRENE': ['ISOP'],
        'O-XYLENE': ['XYL'],
        'WS': ['WSPD10'],
        'TEMP': ['TEMP2'],
        'WD': ['WDIR10'],
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
    to_crn = {'SUR_TEMP': ['TEMPG'], 'T_HR_AVG': ['TEMP2'], 'SOLARAD': ['RGRND'], 'SOIL_MOISTURE_5': ['SOIM1'], 'SOIL_MOISTURE_10': ['SOIM2']}
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
        'AALJ',
        'AALK1J',
        'AALK2J',
        'ABNZ1J',
        'ABNZ2J',
        'ABNZ3J',
        'ACAJ',
        'ACLJ',
        'AECJ',
        'AFEJ',
        'AISO1J',
        'AISO2J',
        'AISO3J',
        'AKJ',
        'AMGJ',
        'AMNJ',
        'ANAJ',
        'ANH4J',
        'ANO3J',
        'AOLGAJ',
        'AOLGBJ',
        'AORGCJ',
        'AOTHRJ',
        'APAH1J',
        'APAH2J',
        'APAH3J',
        'APNCOMJ',
        'APOCJ',
        'ASIJ',
        'ASO4J',
        'ASQTJ',
        'ATIJ',
        'ATOL1J',
        'ATOL2J',
        'ATOL3J',
        'ATRP1J',
        'ATRP2J',
        'AXYL1J',
        'AXYL2J',
        'AXYL3J',
        'AORGAJ',
        'AORGPAJ',
        'AORGBJ',
    ]
)
aitken = array(['ACLI', 'AECI', 'ANAI', 'ANH4I', 'ANO3I', 'AOTHRI', 'APNCOMI', 'APOCI', 'ASO4I', 'AORGAI', 'AORGPAI', 'AORGBI'])
coarse = array(['ACLK', 'ACORS', 'ANH4K', 'ANO3K', 'ASEACAT', 'ASO4K', 'ASOIL'])
noy_gas = array(['NO', 'NO2', 'NO3', 'N2O5', 'HONO', 'HNO3', 'PAN', 'PANX', 'PNA', 'NTR', 'CRON', 'CRN2', 'CRNO', 'CRPX', 'OPAN'])
pec = array(['AECI', 'AECJ'])
pso4 = array(['ASO4I', 'ASO4J'])
pno3 = array(['ANO3I', 'ANO3J'])
pnh4 = array(['ANH4I', 'ANH4J'])
pcl = array(['ACLI', 'ACLJ'])
poc = array(
    [
        'AOTHRI',
        'APNCOMI',
        'APOCI',
        'AORGAI',
        'AORGPAI',
        'AORGBI',
        'ATOL1J',
        'ATOL2J',
        'ATOL3J',
        'ATRP1J',
        'ATRP2J',
        'AXYL1J',
        'AXYL2J',
        'AXYL3J',
        'AORGAJ',
        'AORGPAJ',
        'AORGBJ',
        'AOLGAJ',
        'AOLGBJ',
        'AORGCJ',
        'AOTHRJ',
        'APAH1J',
        'APAH2J',
        'APAH3J',
        'APNCOMJ',
        'APOCJ',
        'ASQTJ',
        'AISO1J',
        'AISO2J',
        'AISO3J',
        'AALK1J',
        'AALK2J',
        'ABNZ1J',
        'ABNZ2J',
        'ABNZ3J',
        'AORGAI',
        'AORGAJ',
        'AORGPAI',
        'AORGPAJ',
        'AORGBI',
        'AORGBJ',
    ]
)
minerals = array(['AALJ', 'ACAJ', 'AFEJ', 'AKJ', 'AMGJ', 'AMNJ', 'ANAJ', 'ATIJ', 'ASIJ'])
