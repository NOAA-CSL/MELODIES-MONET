# SPDX-License-Identifier: Apache-2.0
#
# File started by Maggie Bruckner. 
# Contains satellite specific pairing operators
import numpy as np
import pandas as pd
import xarray as xr

import warnings
import logging

numba_logger = logging.getLogger("numba")
numba_logger.setLevel(logging.WARNING)

def vertical_regrid(input_press, input_values, output_press):
    '''
    This function uses interp1d to regrid vertical layers in a 3D array
    
    Function requires:
        input_press = input pressure levels in hPa and same dimensions as input_values (lon, lat, alt)
        input_values = Dataarray of input values to be regridded (lon, lat, alt)
        output_press = output pressure levels in hPa, dimensions are the same as input values, except for the altitude (lon, lat, newalt)
        
    Function Returns:
        regrid_array = the data regridded to the new pressure levels

    '''
    from scipy import interpolate
    
    out_array = np.full_like(output_press,np.nan)
    for y in range (input_press.shape[0]):
        # Longitude values
        for x in range (input_press.shape[1]):
            xx = input_press[y,x,:]
            yy = input_values[y,x,:]
            xnew = output_press[y,x,:]
            f = interpolate.interp1d(xx, yy, fill_value="extrapolate")

            out_array[y,x,:] = f(xnew)
    return out_array

def mod_to_overpasstime(modobj,opass_tms,partial_col=None):
    '''
    Interpolate model to satellite overpass time.

    Parameters
    ----------

    modobj : xarray.Dataset
        model data
    opass_tms : pandas.DatetimeIndex
        satellite overpass local time
    partial_col : str
        variable to calculate partial columns for
    Output
    ------
    outmod : xarray.Dataset 
        revised model data at local overpass time
    '''
    from .tools import calc_geolocaltime
    
    nst, = opass_tms.shape

    # initialize local time as variable
    modobj['localtime'] = calc_geolocaltime(modobj)

    # initialize new model object with satellite datetimes
    outmod = []

    for ti in np.arange(nst):
        # Apply filter to select model data within +/- 1 output time step of the overpass time
        tempmod = modobj.where(np.abs(modobj['localtime'] - opass_tms[ti].to_datetime64()) < (modobj.time[1] - modobj.time[0]))
        
        # determine factors for linear interpolation in time
        tfac = 1 - (np.abs(tempmod['localtime'] - opass_tms[ti].to_datetime64())/(modobj.time[1] - modobj.time[0]))
        tempmod = tempmod.drop_vars('localtime')
        # Carry out time interpolation
        ## Note regarding current behavior: will only carry out time interpolation if at least 2 model timesteps
        outmod.append((tfac*tempmod).sum(dim='time', min_count=2,keep_attrs=True))
    #print(outmod)
    outmod = xr.concat(outmod,dim='time')
    outmod['time'] = (['time'],opass_tms)
    
    if partial_col:
        from melodies_monet.util.tools import calc_partialcolumn        
        outmod[f'{partial_col}_col'] = calc_partialcolumn(outmod,var=partial_col)
        
    return outmod

def check_timestep(model_data,obs_data):
    ''' When pairing to level 3 data, model data may need to be aggregated to observation timestep.
        This function checks if the model data and observation data have the same timestep. Model data 
        is aggregated to observation timestep. Assumes level 3 data has a monthly or daily timestep and 
        that the model data is higher frequency or same frequency.
    '''

    # check if l3 is daily
    timestep = xr.infer_freq(obs_data.time.dt.round('D'))
    # if not daily, check if l3 is monthly
    if timestep != 'D':
        timestep = xr.infer_freq(pd.to_datetime(obs_data.time.dtstrftime('%Y-%m')))
    if timestep == 'D' or timestep == 'MS':
        print('Aggregating model to observation timestep')
        return model_data.resample(time=timestep).mean()
    else:
        print('Timestep check and model resample failed')
        raise

def mopitt_l3_pairing(model_data,obs_data,co_ppbv_varname,global_model=True):
    ''' Calculate model CO column, with MOPITT averaging kernel applied.
    '''
    try:
        import xesmf as xe
    except ImportError:
        print('satellite_utilities: xesmf module not found')
        raise
    
    ## Check if obs are monthly or daily
    if obs_data.attrs['monthly']:
        # if obs_data is monthly, take monthly mean of model data
        model_obstime = model_data.resample(time='MS').mean()
        filtstr = '%Y-%m'
    elif not obs_data.attrs['monthly']:
        # obs_data is daily, so model and obs seem to be on same time step
        model_obstime = model_data
        filtstr = '%Y-%m-%d'
    else:
        # check frequency of model data 
        # Should not get here.
        tstep = xr.infer_freq(model_data.time.dt.round('D'))
        if tstep == 'MS' or tstep == 'M':
            model_obstime = model_data
            filtstr = '%Y-%m'
        else:
            print('Time resolution of model data and MOPITT data is incompatible')
            raise
        
    # initialize regridder for horizontal interpolation 
    # from model grid to MOPITT grid
    grid_adjust = xe.Regridder(model_obstime[['latitude','longitude']],obs_data[['lat','lon']],
                               'bilinear',periodic=global_model,unmapped_to_nan=True)
    co_model_regrid = grid_adjust(model_obstime[co_ppbv_varname])
    pressure_model_regrid = grid_adjust(model_obstime['pres_pa_mid']/100.)
    
    # enforce dimension order as (time,lat,lon,z)
    co_model_regrid = co_model_regrid.transpose('time','lon','lat','z')
    pressure_model_regrid = pressure_model_regrid.transpose('time','lon','lat','z')
    
    # vertical regrid of model to satellite
    co_regrid = xr.full_like(obs_data['pressure'], np.nan)
    # MEB: loop over time outside of regrid lowers memory usage
    for t in range(obs_data.time.size):
        obs_day = obs_data.time[t].dt.strftime(filtstr)
        co_regrid[t] = vertical_regrid(pressure_model_regrid.sel(time=obs_day).values.squeeze(), 
                                       co_model_regrid.sel(time=obs_day).values.squeeze(), 
                                       obs_data['pressure'][t].values)
    
    # apply AK
    ## log apriori and model data
    log_ap = np.log10(obs_data['apriori_prof'])
    log_mod = np.log10(co_regrid)
    diff_arr = log_mod-log_ap
    ## smooth/apply ak
    smoothed = obs_data['apriori_col'] + (obs_data['ak_col']*diff_arr).sum(dim='alt', min_count=1)
    
    # Add variable name to smoothed model dataarray, combine with obs_data
    smoothed = smoothed.rename(co_ppbv_varname+'_column_model')
    ds = xr.merge([smoothed,obs_data.copy(deep=True)]) 
    
    # Apply scaling to drop scientific notation (x10^{18} molec/cm2 instead of molec/cm2)
    ##  Taylor plot doesn't work if don't do this.
    ds[co_ppbv_varname+'_column_model'] /= 1e18
    ds[co_ppbv_varname+"_column_model"] = ds[co_ppbv_varname+'_column_model'].assign_attrs(units='$10^{18} molec./cm^{2}$')
    ds['column'] /= 1e18
    ds["column"] = ds['column'].assign_attrs(units='$10^{18} molec./cm^{2}$')
    
    # rename dims from lon/lat to x/y for consistency with other datasets
    ds = ds.rename_dims({'lat':'x','lon':'y'})
    # Makde lat/lon coordinates 2d
    lat_2d,lon_2d = np.meshgrid(ds.lat,ds.lon)
    ds['latitude'] = (['y','x'],lat_2d)
    ds['longitude'] = (['y','x'],lon_2d)
    ds = ds.reset_coords().set_coords(['latitude','longitude','time','alt'])
    return ds    

def omps_l3_daily_o3_pairing(model_data,obs_data,ozone_ppbv_varname):
    '''Calculate model ozone column from model ozone profile in ppbv. Move data from model grid 
        to 1x1 degree OMPS L3 data grid. Following data grid matching, take daily mean for model data.
    '''
    try:
        import xesmf as xe
    except ImportError:
        print('satellite_utilities: xesmf module not found')
        raise
    
    # factor for converting ppbv profiles to DU column
    # also requires conversion of dp from Pa to hPa
    du_fac = 1.0e-5*6.023e23/28.97/9.8/2.687e19
    column = (du_fac*(model_data['dp_pa']/100.)*model_data[ozone_ppbv_varname]).sum('z')
    
    # initialize regrid and apply to column data
    grid_adjust = xe.Regridder(model_data[['latitude','longitude']],obs_data[['latitude','longitude']],'bilinear',periodic=True)
    mod_col_obsgrid = grid_adjust(column)
    # Aggregate time-step to daily means
    daily_mean = mod_col_obsgrid.resample(time='1D').mean()
    # change dimension name for date to time
    daily_mean = daily_mean.rename(ozone_ppbv_varname)

    return xr.merge([daily_mean,obs_data])

def space_and_time_pairing(model_data,obs_data,pair_variables):
    '''Bilinear spatial and temporal satellite pairing code. 
    Assumes model data has (time,pressure,latitude,longitude) dimensions.
    Assumes observation data contains fields named time, pressure, latiutde, and longitude.
    
    
    *** need to make setup work for surface/1z fields, as some pairing requires surface pressure field *** 
    '''
    try:
        import xesmf as xe
    except ImportError:
        print('satellite_utilities: xesmf module not found')
        raise
    mod_nf,mod_nz,mod_nx,mod_ny = model_data[pair_variables[0]].shape # assumes model data is structured (time,z,lon,lat). lon/lat dimension order likely unimportant
    # obs_nz = obs_data['pressure'].shape # assumes 1d pressure field in observation set
    obs_nx,obs_ny = obs_data['longitude'].shape # assumes 2d lat/lon fields in observation set
    # initialize dictionary and arrays for interpolated model data
    ds = {i:np.zeros((mod_nz,obs_nx,obs_ny)) for i in pair_variables}
    
    # loop over model time steps
    for f in range(mod_nf):
        
        # set index for observation data less than 1 model timestep from working model file.
        tindex = np.where(np.abs(obs_data.time - model_data.time[f]) <= (model_data.time[1]-model_data.time[0]))[0]
        
        # if there is observation data within the selected time range, proceed with pairing
        if len(tindex):
            # initialize spatial regridder (model lat/lon to satellite swath lat/lon)
            # dimensions of new variables will be (time, z, satellite_x, satellite_y)
            regridr = xe.Regridder(model_data.isel(time=f),obs_data[['latitude','longitude']].sel(x=tindex),'bilinear') # standard bilinear spatial regrid. 
            
            # regrid for each variable in pair_variables
            for j in pair_variables:
                interm_var = regridr(model_data[j][f])
                
                # apply  time interpolation
                if f == (mod_nf-1):
                #    print('last')
                    t2 = np.where((obs_data.time[tindex] >= model_data.time[f]))[0]
                    ds[j][:,tindex[t2]] = interm_var[:,t2].values

                    tind_2 = np.where((obs_data.time[tindex] < model_data.time[f]) & 
                                      (np.abs(obs_data.time[tindex] - model_data.time[f]) <= (model_data.time[1]-model_data.time[0])))[0]
                    tfac1 = 1-(np.abs(model_data.time[f] - obs_data.time[tindex[tind_2]])/(model_data.time[1]-model_data.time[0]))

                    ds[j][:,tindex[tind_2]] += np.expand_dims(tfac1.values,axis=1)*interm_var[:,tind_2].values
                
                elif f == (0):
                #    print('first')
                    t2 = np.where((obs_data.time[tindex] <= model_data.time[f]))[0]
                    ds[j][:,tindex[t2],:] = interm_var[:,t2].values
                    
                    tind_2 = np.where((obs_data.time[tindex] > model_data.time[f]) & 
                                      (np.abs(obs_data.time[tindex] - model_data.time[f]) <= (model_data.time[1]-model_data.time[0])))[0]
                    tfac1 = 1-(np.abs(model_data.time[f] - obs_data.time[tindex[tind_2]])/(model_data.time[1]-model_data.time[0]))

                    ds[j][:,tindex[tind_2],:] += np.expand_dims(tfac1.values,axis=1)*interm_var[:,tind_2,:].values
                   
                else:


                    tfac1 = 1-(np.abs(model_data.time[f] - obs_data.time[tindex])/(model_data.time[1]-model_data.time[0]))
                    
                    ds[j][:,tindex,:] += np.expand_dims(tfac1.values,axis=1)*interm_var.values
    return ds

def omps_nm_pairing(model_data,obs_data,ozone_ppbv_varname):
    'Pairs model ozone mixing ratio with OMPS nadir mapper retrievals. Calculates column without applying apriori'
 
    print('pairing without applying averaging kernel')

    if len(ozone_ppbv_varname) != 1:
        print('ozone_ppbv_varname has more than one entry')

    
    du_fac = 1.0e-5*6.023e23/28.97/9.8/2.687e19 # conversion factor; moves model from ppbv to dobson
    pair_variables = ['dp_pa',ozone_ppbv_varname[0]]
    paired_ds = space_and_time_pairing(model_data,obs_data,pair_variables)
    
    # calculate ozone column, no averaging kernel or apriori applied.
    col = np.nansum(du_fac*(paired_ds['dp_pa']/100.)*paired_ds[ozone_ppbv_varname[0]],axis=0) # new dimensions will be (satellite_x, satellite_y)
    ds = xr.Dataset({ozone_ppbv_varname[0]: (['time','y'],col),
                     'ozone_column':(['time','y'],obs_data.ozone_column.values)
                               },
                    coords={
                        'longitude':(['time','y'],obs_data['longitude'].values),
                        'latitude':(['time','y'],obs_data['latitude'].values),
                        'time':(['time'],obs_data.time.values),
                    })    

    return ds
                                                                            
                                                                            

def calculate_omps_dp(swath_data,model_surface_pressure):
    ''' Calculates OMPS vertical layer thickness in hPa. Takes into account model surface pressure

    Parameters
    ----------
    swath_data : xr.Dataset
        OMPS data. Must include pressure
    model_surface_pressure : xr.DataArray
        surface pressure from model in Pa
    Returns
    -------
    xr.DataArray
        OMPS vertical layer thickness in hPa
    
    '''
    dp_swath_hPa = xr.full_like(swath_data.apriori,np.nan)
    
    down = swath_data.pressure.roll(z=-1)
    up = swath_data.pressure.roll(z=1)
    down[-1] = 0
    
    dp_swath_hPa[:,:,:] = ((up-swath_data['pressure'])/2+(swath_data['pressure']-down)/2) #.values
    dp_swath_hPa[0,:,:] = ((swath_data['pressure']-down)[0]/2) + (model_surface_pressure/100.-swath_data['pressure'][0])
    dp_swath_hPa[-1,:,:] = ((up-swath_data['pressure'])[-1]/2+(swath_data['pressure']-down)[-1]).values
    
    return dp_swath_hPa                                                                                  

def omps_nm_pairing_apriori(model_data,obs_data,ozone_ppbv_varname):
    'Pairs model ozone mixing ratio data with OMPS nm. Applies satellite apriori column to model observations.'
    try:
        import xesmf as xe
        import stratify
    except ImportError:
        print('satellite_utilities: xesmf module not found')
        raise

    du_fac = 1.0e-5*6.023e23/28.97/9.8/2.687e19 # conversion factor; moves model from ppbv to dobson
    
    print('pairing with averaging kernel application')

    paired_ds = xr.Dataset({ozone_ppbv_varname[0]: (['time','y'],np.zeros_like(obs_data.ozone_column.values)),
                 'ozone_column':(['time','y'],obs_data.ozone_column.values)
                           },
                coords={
                    'longitude':(['time','y'],obs_data['longitude'].values),
                    'latitude':(['time','y'],obs_data['latitude'].values),
                    'time':(['time'],obs_data.time.values),
                })
    # use model datestamps to loop over
    ## once omps reader is brought in-line with newer readers, follow similar structure to tropomi with ak application
    model_dates = model_data['time'].dt.floor("D")
    obs_dates = obs_data['time'].dt.floor('D')
    for d,day in enumerate(model_dates):
        obs_day = obs_data.where( obs_dates == day,drop = True).transpose('z','x','y')
        if obs_day.sizes['x'] == 0:
            warnings.warn(f'Observations does not contain data for {day}, skipping.')
            continue
        model_day = model_data.isel(time=d)
        
        # horizontal regridding of model to obs
        rgd_fn = xe.Regridder(model_day[['latitude','longitude']], obs_day[['latitude','longitude']],method="bilinear", periodic=True)
        mod_regrid = rgd_fn(model_day)

        # vertical regridding
        obs_day['dp_hPa'] = calculate_omps_dp(obs_day,mod_regrid['surfpres_pa'])
        mod_regrid = stratify.interpolate(obs_day.pressure*100,mod_regrid['pres_pa_mid'],mod_regrid[ozone_ppbv_varname[0]].values,axis=0)
        # partial column calc
        mod_o3_pcol = mod_regrid*obs_day['dp_hPa']*du_fac
        # apply ak and calculate total columns
        mod_o3_totcol = obs_day['apriori'].sum('z') + (obs_day['layer_efficiency'] * (mod_o3_pcol - obs_day['apriori'])).sum('z')

        paired_ds[ozone_ppbv_varname[0]].loc[dict(time=obs_day.time)] = mod_o3_totcol
  
    return paired_ds