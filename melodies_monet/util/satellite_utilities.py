# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#
# File started by Maggie Bruckner. 
# Contains satellite specific pairing operators
import numpy as np
from datetime import datetime,timedelta

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

def check_timestep(model_data,obs_data):
    ''' When pairing to level 3 data, model data may need to be aggregated to observation timestep.
        This function checks if the model data and observation data have the same timestep. Model data 
        is aggregated to observation timestep. Assumes level 3 data has a monthly or daily timestep and 
        that the model data is higher frequency or same frequency.
    '''
    import xarray as xr
    import pandas as pd
    
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

def mopitt_l3_pairing(model_data,obs_data,co_ppbv_varname):
    ''' Calculate model CO column, with MOPITT averaging kernel applied.
    '''
    import xarray as xr
    try:
        import xesmf as xe
    except ImportError as e:
        print('satellite_utilities: xesmf module not found')
        raise
    
    # Aggregate time-step, if needed
    ## Check if same number of timesteps:
    if obs_data.time.size == model_data.time.size:
        model_obstime = model_data
    elif obs_data.time.size < model_data.time.size and obs_data.time.size >2:
        model_obstime = check_timestep(model_data,obs_data)
    elif obs_data.time.size < model_data.time.size and obs_data.time.size < 3:
        print('Model data and obs data timesteps do not match, and there are not enough observation timesteps to infer the spacing from.')
        raise
    elif obs_data.time.size > mod_data.time.size:
        print('Observation data appears to be a finer time resolution than model data')
        raise
    # initialize regridder for horizontal interpolation 
    # from model grid to MOPITT grid
    grid_adjust = xe.Regridder(model_obstime[['latitude','longitude']],obs_data[['lat','lon']],'bilinear')
    co_model_regrid = grid_adjust(model_obstime[co_ppbv_varname])
    pressure_model_regrid = grid_adjust(model_obstime['pres_pa_mid']/100.)
    
    # enforce dimension order as (time,lat,lon,z)
    co_model_regrid = co_model_regrid.transpose('time','lon','lat','z')
    pressure_model_regrid = pressure_model_regrid.transpose('time','lon','lat','z')
    
    # vertical regrid of model to satellite
    co_regrid = xr.full_like(obs_data['pressure'], np.nan)
    # MEB: loop over time outside of regrid lowers memory usage
    for t in range(obs_data.time.size):
        co_regrid[t] = vertical_regrid(pressure_model_regrid[t].values, co_model_regrid[t].values, obs_data['pressure'][t].values)
    
    # apply AK
    ## log apriori and model data
    log_ap = np.log10(obs_data['apriori_prof'])
    log_mod = np.log10(co_regrid)
    diff_arr = log_mod-log_ap
    ## smooth/apply ak
    smoothed = obs_data['apriori_col'] + (obs_data['ak_col']*diff_arr).sum(dim='alt')
    
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
        import xarray as xr
        import xesmf as xe
    except ImportError as e:
        print('satellite_utilities: xesmf module not found')
        raise
    
    # factor for converting ppbv profiles to DU column
    # also requires conversion of dp from Pa to hPa
    du_fac = 1.0e-5*6.023e23/28.97/9.8/2.687e19
    column = (du_fac*(model_data['dp_pa']/100.)*model_data[ozone_ppbv_varname]).sum('z')
    
    # initialize regrid and apply to column data
    grid_adjust = xe.Regridder(model_data[['latitude','longitude']],obs_data[['latitude','longitude']],'bilinear')
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
    except ImportError as e:
        print('satellite_utilities: xesmf module not found')
        raise
    mod_nf,mod_nz,mod_nx,mod_ny = model_data[pair_variables[0]].shape # assumes model data is structured (time,z,lon,lat). lon/lat dimension order likely unimportant
    obs_nz = obs_data['pressure'].shape # assumes 1d pressure field in observation set
    obs_nx,obs_ny = obs_data['longitude'].shape # assumes 2d lat/lon fields in observation ser
    
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
    import xarray as xr
    import pandas as pd
 
    print('pairing without applying averaging kernel')

    if len(ozone_ppbv_varname) != 1:
        print('ozone_ppbv_varname has more than one entry')

    
    du_fac = 1.0e-5*6.023e23/28.97/9.8/2.687e19 # conversion factor; moves model from ppbv to dobson
    pair_variables = ['dp_pa',ozone_ppbv_varname]
    paired_ds = space_and_time_pairing(model_data,obs_data,pair_variables)
    
    # calculate ozone column, no averaging kernel or apriori applied.
    col = np.nansum(du_fac*(paired_ds['dp_pa']/100.)*paired_ds['o3vmr'],axis=0) # new dimensions will be (satellite_x, satellite_y)
    ds = xr.Dataset({ozone_ppbv_varname[0]: (['time','y'],col),
                     'ozone_column':(['time','y'],obs_data.ozone_column.values)
                               },
                    coords={
                        'longitude':(['time','y'],obs_data['longitude'].values),
                        'latitude':(['time','y'],obs_data['latitude'].values),
                        'time':(['time'],obs_data.time.values),
                    })    

    return ds
                                                                            
                                                                            

def omps_nm_pairing_apriori(model_data,obs_data,ozone_ppbv_varname):
    'Pairs model ozone mixing ratio data with OMPS nm. Applies satellite apriori column to model observations.'
    import xarray as xr
    import pandas as pd
    try:
        import xesmf as xe
    except ImportError as e:
        print('satellite_utilities: xesmf module not found')
        raise

    du_fac = 1.0e-5*6.023e23/28.97/9.8/2.687e19 # conversion factor; moves model from ppbv to dobson
    
    print('pairing with averaging kernel application')
                     
    # Grab necessary shape information
    nf,nz_m,nx_m,ny_m = model_data[ozone_ppbv_varname[0]].shape
    nx,ny = obs_data.ozone_column.shape
    ## initialize intermediates for use in calcluating column
    pressure_temp = np.zeros((nz_m,nx,ny))
    ozone_temp = np.zeros((nz_m,nx,ny))
    sfc = np.zeros((nx,ny))
    ## loop over model time steps
    for f in range(nf):
        
        tindex = np.where(np.abs(obs_data.time - model_data.time[f]) <= (model_data.time[1]-model_data.time[0]))[0]
        if len(tindex):
            # regrid spatially (model lat/lon to satellite swath lat/lon)
            regridr = xe.Regridder(model_data.isel(time=f),obs_data[['latitude','longitude']].sel(x=tindex),'bilinear')
            regrid_oz = regridr(model_data[ozone_ppbv_varname[0]][f])
            regrid_p = regridr(model_data['pres_pa_mid'][f]) # this one should be pressure variable (for the interpolation).
            sfp = regridr(model_data['surfpres_pa'][f])
            # fixes for observations before/after model time range.
            if f == (nf-1):
                t2 = np.where((obs_data.time[tindex] >= model_data.time[f]))[0]
                ozone_temp[:,tindex[t2],:] = regrid_oz[:,t2,:].values
                pressure_temp[:,tindex[t2],:] = regrid_p[:,t2,:].values
                sfc[t2,:] = sfp[t2,:].values 
                tind_2 = np.where((obs_data.time[tindex] < model_data.time[f]) & 
                                  (np.abs(obs_data.time[tindex] - model_data.time[f]) <= (model_data.time[1]-model_data.time[0])))[0]
                tfac1 = 1-(np.abs(model_data.time[f] - obs_data.time[tindex[tind_2]])/(model_data.time[1]-model_data.time[0]))

                ozone_temp[:,tindex[tind_2],:] += np.expand_dims(tfac1.values,axis=1)*regrid_oz[:,tind_2,:].values
                pressure_temp[:,tindex[tind_2],:] += np.expand_dims(tfac1.values,axis=1)*regrid_p[:,tind_2,:].values
                sfc[tindex[tind_2],:] += np.expand_dims(tfac1.values,axis=1)*sfp[tind_2,:].values
            elif f == 0:
                t2 = np.where((obs_data.time[tindex] <= model_data.time[f]))[0]
                ozone_temp[:,tindex[t2],:] = regrid_oz[:,t2,:].values
                pressure_temp[:,tindex[t2],:] = regrid_p[:,t2,:].values
                sfc[tindex[t2],:] = sfp[t2,:].values 
                tind_2 = np.where((obs_data.time[tindex] > model_data.time[f]) & 
                                  (np.abs(obs_data.time[tindex] - model_data.time[f]) <= (model_data.time[1]-model_data.time[0])))[0]
                tfac1 = 1-(np.abs(model_data.time[f] - obs_data.time[tindex[tind_2]])/(model_data.time[1]-model_data.time[0]))
                ozone_temp[:,tindex[tind_2],:] += np.expand_dims(tfac1.values,axis=1)*regrid_oz[:,tind_2,:].values
                pressure_temp[:,tindex[tind_2],:] += np.expand_dims(tfac1.values,axis=1)*regrid_p[:,tind_2,:].values
                sfc[tind_2,:] += np.expand_dims(tfac1.values,axis=1)*sfp[tind_2,:].values
            else:
                tfac1 = 1-(np.abs(model_data.time[f] - obs_data.time[tindex])/(model_data.time[1]-model_data.time[0]))
                ozone_temp[:,tindex,:] += np.expand_dims(tfac1.values,axis=1)*regrid_oz.values
                pressure_temp[:,tindex,:] += np.expand_dims(tfac1.values,axis=1)*regrid_p.values
                sfc[tindex,:] += np.expand_dims(tfac1.values,axis=1)*sfp.values
    # Interpolate model data to satellite pressure levels
    from wrf import interplevel
    # note: for interpolation in pressure coordinates to work, z dimension must be such that the smallest 
    # pressure is on the bottom. With Melodies-Monet model datasets, this requires flipping the z dimension 
    # as the model readers are set up to ensure the surface is at index 0. 
    ozone_satp = interplevel(ozone_temp[::-1],pressure_temp[::-1]/100.,obs_data.pressure,missing=np.nan)
    ozone_satp = ozone_satp.values
    
    ozone_satp[np.isnan(ozone_satp)] = 0
    oz = np.zeros_like(obs_data.ozone_column.values)
    
    nl,n1,n2 = ozone_satp.shape
    
    # delta pressure calculation for satellite pressure midlevels
    p = obs_data.pressure.values
    shift_down = np.roll(p,-1)
    shift_down[-1] =0

    shift_up = np.roll(p,1)
    band = (shift_up-p)/2+(p-shift_down)/2
   
    band[0] = (p-shift_down)[0]/2

    band[-1] = (shift_up-p)[-1]/2 + (p-shift_down)[-1]
    for i in range(nl):
        
        if i != 0:
            dp = band[i]
        else:
            sfc[sfc == 0] = np.nan
            dp = np.abs(sfc/100. - obs_data.pressure[i].values) + band[i]

        add = du_fac*dp*ozone_satp[i]
        eff = obs_data.layer_efficiency[:,:,i].values
        ap = obs_data.apriori[:,:,i].values
        oz = oz + ap*(1-eff) + (eff)*(add)
 
    ds = xr.Dataset({ozone_ppbv_varname[0]: (['time','y'],oz),
                     'ozone_column':(['time','y'],obs_data.ozone_column.values)
                               },
                    coords={
                        'longitude':(['time','y'],obs_data['longitude'].values),
                        'latitude':(['time','y'],obs_data['latitude'].values),
                        'time':(['time'],obs_data.time.values),
                    })
    return ds
