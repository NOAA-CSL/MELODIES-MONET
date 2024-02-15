# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#
from __future__ import division

from builtins import range

import numpy as np

__author__ = 'barry'


def search_listinlist(array1, array2):
    # find intersections

    s1 = set(array1.flatten())
    s2 = set(array2.flatten())

    inter = s1.intersection(s2)

    index1 = np.array([])
    index2 = np.array([])
    # find the indexes in array1
    for i in inter:
        index11 = np.where(array1 == i)
        index22 = np.where(array2 == i)
        index1 = np.concatenate([index1[:], index11[0]])
        index2 = np.concatenate([index2[:], index22[0]])

    return np.sort(np.int32(index1)), np.sort(np.int32(index2))


def list_contains(list1, list2):
    """Return True if any item in `list1` is also in `list2`."""
    for m in list1:
        for n in list2:
            if n == m:
                return True

    return False


def linregress(x, y):
    import statsmodels.api as sm

    xx = sm.add_constant(x)
    model = sm.OLS(y, xx)
    fit = model.fit()
    b, a = fit.params[0], fit.params[1]
    rsquared = fit.rsquared
    std_err = np.sqrt(fit.mse_resid)
    return a, b, rsquared, std_err


def findclosest(list, value):
    """Return (index, value) of the closest value in `list` to `value`."""
    a = min((abs(x - value), x, i) for i, x in enumerate(list))
    return a[2], a[1]


def _force_forder(x):
    """
    Converts arrays x to fortran order. Returns
    a tuple in the form (x, is_transposed).
    """
    if x.flags.c_contiguous:
        return (x.T, True)
    else:
        return (x, False)


def kolmogorov_zurbenko_filter(df, col, window, iterations):
    """KZ filter implementation
        series is a pandas series
        window is the filter window m in the units of the data (m = 2q+1)
        iterations is the number of times the moving average is evaluated
        """
    df.index = df.time_local
    z = df.copy()
    for i in range(iterations):
        z.index = z.time_local
        z = z.groupby('siteid')[col].rolling(
            window, center=True, min_periods=1).mean().reset_index().dropna()
    df = df.reset_index(drop=True)
    return df.merge(z, on=['siteid', 'time_local'])


def wsdir2uv(ws, wdir):
    from numpy import pi, sin, cos
    u = -ws * sin(wdir * pi / 180.)
    v = -ws * cos(wdir * pi / 180.)
    return u, v


def get_relhum(temp, press, vap):
    # temp:  temperature (K)
    # press: pressure (Pa)
    # vap:   water vapor mixing ratio (kg/kg)
    temp_o = 273.16
    es_vap = 611.0 * np.exp(17.67 * ((temp - temp_o) / (temp - 29.65)))
    ws_vap = 0.622 * (es_vap / press)
    relhum = 100.0 * (vap / ws_vap)
    return relhum


def long_to_wide(df):
    from pandas import Series, merge
    w = df.pivot_table(values='obs',
                       index=['time', 'siteid'],
                       columns='variable').reset_index()
    cols = Series(df.columns)
    g = df.groupby('variable')
    for name, group in g:
        w[name + '_unit'] = group.units.unique()[0]
    # mergeon = hstack((index.values, df.variable.unique()))
    return merge(w, df, on=['siteid', 'time'])


def calc_8hr_rolling_max(df, col=None, window=None):
    df.index = df.time_local
    df_rolling = df.groupby('siteid')[col].rolling(
        window, center=True, win_type='boxcar').mean().reset_index().dropna()
    df_rolling_max = df_rolling.groupby('siteid').resample(
        'D', on='time_local').max().reset_index(drop=True)
    df = df.reset_index(drop=True)
    return df.merge(df_rolling_max, on=['siteid', 'time_local'])


def calc_24hr_ave(df, col=None):
    df.index = df.time_local
    df_24hr_ave = df.groupby('siteid')[col].resample('D').mean().reset_index()
    df = df.reset_index(drop=True)
    return df.merge(df_24hr_ave, on=['siteid', 'time_local'])


def calc_3hr_ave(df, col=None):
    df.index = df.time_local
    df_3hr_ave = df.groupby('siteid')[col].resample('3H').mean().reset_index()
    df = df.reset_index(drop=True)
    return df.merge(df_3hr_ave, on=['siteid', 'time_local'])


def calc_annual_ave(df, col=None):
    df.index = df.time_local
    df_annual_ave = df.groupby('siteid')[col].resample(
        'A').mean().reset_index()
    df = df.reset_index(drop=True)
    return df.merge(df_annual_ave, on=['siteid', 'time_local'])


def get_giorgi_region_bounds(index=None, acronym=None):
    import pandas as pd
    i = [
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
        21, 22
    ]
    acro = [
        'NAU', 'SAU', 'AMZ', 'SSA', 'CAM', 'WNA', 'CNA', 'ENA', 'ALA', 'GRL',
        'MED', 'NEU', 'WAF', 'EAF', 'SAF', 'SAH', 'SEA', 'EAS', 'SAS', 'CAS',
        'TIB', 'NAS'
    ]
    lonmax = [
        155, 155, -34, -40, -83, -103, -85, -60, -103, -10, 40, 40, 22, 52, 52,
        65, 155, 145, 100, 75, 100, 180
    ]
    lonmin = [
        110, 110, -82, -76, -116, -130, -103, -85, -170, -103, -10, -10, -20,
        22, -10, -20, 95, 100, 65, 40, 75, 40
    ]
    latmax = [
        -11, -28, 12, -20, 30, 60, 50, 50, 72, 85, 48, 75, 18, 18, -12, 30, 20,
        50, 30, 50, 50, 70
    ]
    latmin = [
        -28, -45, -20, -56, 10, 30, 30, 25, 60, 50, 30, 48, -12, -12, -35, 18,
        -11, 20, 5, 30, 30, 50
    ]
    df = pd.DataFrame(
        {
            'latmin': latmin,
            'lonmin': lonmin,
            'latmax': latmax,
            'lonmax': lonmax,
            'acronym': acro
        },
        index=i)
    try:
        if index is None and acronym is None:
            print('either index or acronym needs to be supplied')
            print(
                'look here https://web.northeastern.edu/sds/web/demsos/images_002/subregions.jpg'
            )
            raise ValueError
        elif index is not None:
            return df.loc[df.index == index].values.flatten()
        else:
            return df.loc[df.acronym == acronym.upper()].values.flatten()
    except ValueError:
        exit


def get_giorgi_region_df(df):
    df.loc[:, 'GIORGI_INDEX'] = None
    df.loc[:, 'GIORGI_ACRO'] = None
    for i in range(22):
        latmin, lonmin, latmax, lonmax, acro = get_giorgi_region_bounds(
            index=int(i + 1))
        con = (df.longitude <= lonmax) & (df.longitude >= lonmin) & (
            df.latitude <= latmax) & (df.latitude >= latmin)
        df.loc[con, 'GIORGI_INDEX'] = i + 1
        df.loc[con, 'GIORGI_ACRO'] = acro
    return df


def get_epa_region_bounds(index=None, acronym=None):
    import pandas as pd
    i = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
    acro = [
        'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9', 'R10', 'AK',
        'PR', 'VI'
    ]
    lonmax = [
        -66.8628, -73.8885, -74.8526, -75.4129, -80.5188, -88.7421, -89.1005,
        -96.438, -109.0475, -111.0471, -129.99, -65.177765, -64.26384
    ]
    lonmin = [
        -73.7272, -79.7624, -83.6753, -91.6589, -97.2304, -109.0489, -104.0543,
        -116.0458, -124.6509, -124.7305, -169.9146, -67.289886, -64.861221
    ]
    latmax = [
        47.455, 45.0153, 42.5167, 39.1439, 49.3877, 37.0015, 43.5008, 48.9991,
        42.0126, 49.0027, 71.5232, 18.520551, 18.751244
    ]
    latmin = [
        40.9509, 38.8472, 36.5427, 24.3959, 36.9894, 25.8419, 35.9958, 36.9949,
        31.3325, 41.9871, 52.5964, 17.904834, 18.302014
    ]
    df = pd.DataFrame(
        {
            'latmin': latmin,
            'lonmin': lonmin,
            'latmax': latmax,
            'lonmax': lonmax,
            'acronym': acro
        },
        index=i)
    try:
        if index is None and acronym is None:
            print('either index or acronym needs to be supplied')
            print(
                'Look here for more information: https://www.epa.gov/enviro/epa-regional-kml-download',
                'https://gist.github.com/jakebathman/719e8416191ba14bb6e700fc2d5fccc5'
            )
            raise ValueError
        elif index is not None:
            return df.loc[df.index == index].values.flatten()
        else:
            return df.loc[df.acronym == acronym.upper()].values.flatten()
    except ValueError:
        exit


def get_epa_region_df(df):
    df.loc[:, 'EPA_INDEX'] = None
    df.loc[:, 'EPA_ACRO'] = None
    for i in range(13):
        latmin, lonmin, latmax, lonmax, acro = get_epa_region_bounds(
            index=int(i + 1))
        con = (df.longitude <= lonmax) & (df.longitude >= lonmin) & (
            df.latitude <= latmax) & (df.latitude >= latmin)
        df.loc[con, 'EPA_INDEX'] = i + 1
        df.loc[con, 'EPA_ACRO'] = acro
    return df

def resample_stratify(da, levels, vertical, axis=1,interpolation='linear',extrapolation='nan'):
    import stratify
    import xarray as xr

    result = stratify.interpolate(levels, vertical.chunk().data, da.chunk().data, axis=axis,
                                 interpolation = interpolation,extrapolation = extrapolation)
    dims = da.dims
    out = xr.DataArray(result, dims=dims)
    for i in dims:
        if i != "z":
            out[i] = da[i]
    out.attrs = da.attrs.copy()
    if len(da.coords) > 0:
        for i in da.coords:
            if i != "z":
                out.coords[i] = da.coords[i]
    return out

def vert_interp(ds_model,df_obs,var_name_list):
    import xarray as xr
    from pandas import merge_asof, Series

    var_out_list = []
    for var_name in var_name_list:
        if var_name == 'pressure_model':
            out = resample_stratify(ds_model[var_name],sorted(ds_model.pressure_obs.squeeze().values,reverse=True),
                                      ds_model['pressure_model'],axis=1,
                                      interpolation='linear',extrapolation='nan')
            #Use extrapolation nan for the pressure so that later you can assign the nan values to the pressure_obs value 
            #instead of the midpoint of the edge model cells. This is needed for the pairing later on.
        else:
            out = resample_stratify(ds_model[var_name],sorted(ds_model.pressure_obs.squeeze().values,reverse=True),
                                  ds_model['pressure_model'],axis=1,
                                  interpolation='linear',extrapolation='nearest')
        out.name = var_name
        var_out_list.append(out)

    df_model = xr.merge(var_out_list).to_dataframe().reset_index()
    df_model.fillna({'pressure_model':df_model.pressure_obs},inplace=True)
    df_model.drop(labels=['x','y','z','pressure_obs','time_obs'], axis=1, inplace=True)
    df_model.rename(columns={'pressure_model':'pressure_obs'}, inplace=True)

    final_df_model = merge_asof(df_obs, df_model, 
                            by=['latitude', 'longitude', 'pressure_obs'], 
                            on='time', direction='nearest')

    return final_df_model

def mobile_and_ground_pair(ds_model,df_obs, var_name_list):
    import xarray as xr
    from pandas import merge_asof, Series
    
    var_out_list = []
    # Extract just the surface level data from correct model variables
    # if there is a z dimension, extract the surface, otherwise assume data is at surface and issue warning
    if 'z' in ds_model.dims:
        for var_name in var_name_list:
            out = ds_model[var_name].isel(z=0)
            out.name = var_name
            var_out_list.append(out)
    else:
        print('WARNING: No z dimension in model, assuming all data at surface.')
        for var_name in var_name_list:
            out = ds_model[var_name]
            out.name = var_name
            var_out_list.append(out)
    
    df_model = xr.merge(var_out_list).to_dataframe().reset_index()
    df_model.drop(labels=['x','y','time_obs'], axis=1, inplace=True)

    final_df_model = merge_asof(df_obs, df_model, 
                            by=['latitude', 'longitude'], 
                            on='time', direction='nearest')

    return final_df_model

def find_obs_time_bounds(files=[],time_var=None):
    """Function to read a series of ict files and print a list of min and max times for each.

    Parameters
    ----------
    files : str or iterable
        str or list of str containing filenames that should be read.
        
    time_var : str
        Optional, variable name that should be assumed to be time when reading aircaft csv files.

    Returns
    -------
    bounds : dict
        Dict containing time bounds for each file.

    """
    import os 
    import monetio as mio
    import xarray as xr
    
    if isinstance(files,str):
        files = [files]
    
    bounds = {}
    for file in files:
        _, extension = os.path.splitext(files[0])
        try:
            if extension in {'.nc', '.ncf', '.netcdf', '.nc4'}:
                obs = xr.open_dataset(file)
            elif extension in ['.ict', '.icartt']:
                obs = mio.icartt.add_data(file)
            elif extension in ['.csv']:
                from .read_util import read_aircraft_obs_csv
                obs = read_aircraft_obs_csv(filename=file,time_var=time_var)
            else:
                raise ValueError(f'extension {extension!r} currently unsupported')
        except Exception as e:
            print('something happened opening file:', e)
            return
        
        time_min = obs['time'].min()
        time_max = obs['time'].max()
        
        print('For {}, time bounds are, Min: {}, Max: {}'.format(file,time_min,time_max))
        bounds[file] = {'Min':time_min,'Max':time_max}

        del obs
        
    return bounds

def loop_pairing(control,file_pairs_yaml='',file_pairs={},save_types=['paired']):
    
    """Function to loop over sets of pairings and save them out as multiple netcdf files.

    Parameters
    ----------
    control: str
        str containing path to control file.
    file_pairs : dict (optional)
        Dict containing filenames for obs and models. This should be specified if file_pairs_yaml is not. 
        Example: file_pairs = {'0722':{'model':{'wrfchem_v4.2':'/wrk/users/charkins/melodies-monet_data/wrfchem/run_CONUS_fv19_BEIS_1.0xISO_RACM_v4.2.2_racm_berk_vcp_noI_phot/0722/*'},
                      'obs':{'firexaq':'/wrk/d2/rschwantes/obs/firex-aq/R1/10s_merge/firexaq-mrg10-dc8_merge_20190722_R1.ict'}},
            '0905':{'model':{'wrfchem_v4.2':'/wrk/users/charkins/melodies-monet_data/wrfchem/run_CONUS_fv19_BEIS_1.0xISO_RACM_v4.2.2_racm_berk_vcp_noI_phot_soa/0905/*'},
                    'obs':{'firexaq':'/wrk/d2/rschwantes/obs/firex-aq/R1/10s_merge/firexaq-mrg10-dc8_merge_20190905_R1.ict'}}
            }
    file_pairs_yaml: str
        str containing path to a yaml file with file pairings. 
        Example below
        '0722':
          'model':
            'wrfchem_v4.2':'/wrk/users/charkins/melodies-monet_data/wrfchem/run_CONUS_fv19_BEIS_1.0xISO_RACM_v4.2.2_racm_berk_vcp_noI_phot/0722/*'
          'obs':
            'firexaq':'/wrk/d2/rschwantes/obs/firex-aq/R1/10s_merge/firexaq-mrg10-dc8_merge_20190722_R1.ict'
        '0905':
           'model':
            'wrfchem_v4.2':'/wrk/users/charkins/melodies-monet_data/wrfchem/run_CONUS_fv19_BEIS_1.0xISO_RACM_v4.2.2_racm_berk_vcp_noI_phot_soa/0905/*'
          'obs':
            'firexaq':'/wrk/d2/rschwantes/obs/firex-aq/R1/10s_merge/firexaq-mrg10-dc8_merge_20190905_R1.ict'
    save_types: list (optional)
        List containing the types of data to save to netcdf. Can include any of 'paired', 'models', and 'obs'
    Returns
    -------
    None

    """
    
    from melodies_monet import driver
    
    if file_pairs_yaml:
        import yaml
        with open(file_pairs_yaml, 'r') as stream:
            file_pairs = yaml.safe_load(stream)
    
    for file in file_pairs.keys():
    
        an = driver.analysis()
        an.control=control
        an.read_control()
    
        for model in an.control_dict['model']:
            an.control_dict['model'][model]['files'] = file_pairs[file]['model'][model]
        for obs in an.control_dict['obs']:
            an.control_dict['obs'][obs]['filename'] = file_pairs[file]['obs'][obs]
        
        an.control_dict['analysis']['save']={}
        an.save={}
        for t in save_types:
            an.control_dict['analysis']['save'][t]={'method':'netcdf','prefix':file,'data':'all'}
            an.save[t]={'method':'netcdf','prefix':file,'data':'all'}
        
        an.open_models()
        an.open_obs()
        an.pair_data()
        
        # TEMPORARY sanitize variable names so can be saved to netcdf
        for pair in an.paired.keys():
            datavars = list(an.paired[pair].obj.variables)
            contains_slash = [name for name in datavars if '/' in name]
            contains_and = [name for name in datavars if '&' in name]
            an.paired[pair].obj = an.paired[pair].obj.drop(contains_slash)
            an.paired[pair].obj = an.paired[pair].obj.drop(contains_and)
        
        an.save_analysis()

