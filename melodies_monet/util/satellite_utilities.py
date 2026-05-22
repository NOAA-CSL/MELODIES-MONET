# SPDX-License-Identifier: Apache-2.0
#
# File started by Maggie Bruckner. 
# Contains satellite specific pairing operators
import numpy as np
import pandas as pd
import xarray as xr

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

    nst, = opass_tms.shape
    # nmt, = modobj.time.shape
    # ny,nx = modobj.longitude.shape
    
    # Determine local time offset
    local_utc_offset = (modobj['longitude']/15).round().astype('timedelta64[h]')
    # initialize local time as variable
    modobj['localtime'] = modobj['time'] + local_utc_offset

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
                                                                            
                                                                            

def omps_nm_pairing_apriori(model_data,obs_data,ozone_ppbv_varname):
    'Pairs model ozone mixing ratio data with OMPS nm. Applies satellite apriori column to model observations.'
    try:
        import xesmf as xe
    except ImportError:
        print('satellite_utilities: xesmf module not found')
        raise

    du_fac = 1.0e-5*6.023e23/28.97/9.8/2.687e19 # conversion factor; moves model from ppbv to dobson
    
    print('pairing with averaging kernel application')
                     
    # Grab necessary shape information
    nf,nz_m,nx_m,ny_m = model_data[ozone_ppbv_varname[0]].shape
    nx,ny = obs_data.ozone_column.shape
    ## initialize intermediates for use in calculating column
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


##new code for mapting satellite OMPS NO2 into model grid (5/22/2026; nazrul)

def _standardize_omps_swath_dims_new(obs_data):
    """
    Ensure OMPS swath dimensions are consistent.
    """
    ds = obs_data

    if 'xtrack' in ds.dims and 'y' not in ds.dims:
        ds = ds.rename({'xtrack': 'y'})

    if 'Latitude' in ds and 'latitude' not in ds:
        ds = ds.rename({'Latitude': 'latitude'})
    if 'Longitude' in ds and 'longitude' not in ds:
        ds = ds.rename({'Longitude': 'longitude'})

    return ds


def _get_valid_model_time_indices(model_data, obs_swath):
    import numpy as np

    swath_hours = obs_swath["utc_hour"].values.astype(int)
    model_hours = model_data["time_utc_hour"].values.astype(int)

    valid_idx = []

    for i, mh in enumerate(model_hours):
        diff = np.abs(swath_hours - mh)
        diff = np.minimum(diff, 24 - diff)

        if np.any(diff <= 1):
            valid_idx.append(i)

    return valid_idx


def _select_swath_for_model_hour(obs_swath, model_hour):
    import numpy as np

    swath_hours = obs_swath["utc_hour"].values.astype(int)

    diff = np.abs(swath_hours - model_hour)
    diff = np.minimum(diff, 24 - diff)

    swath_idx = np.where(diff <= 1)[0]

    if len(swath_idx) == 0:
        return None

    return obs_swath.isel(time=swath_idx)


def _compute_pixel_level_satellite_adjustment(model_data_t,obs_swath,obs_no2_var="no2_totalcolumn"):
    """
    Use model vertical shape factor to adjust satellite column using AK + prior.
    model_data_t should contain only ONE model time.
    """

    import numpy as np
    import xesmf as xe
    from scipy import interpolate

    obs_swath = _standardize_omps_swath_dims_new(obs_swath)

    # Regrid model profile to satellite pixels
    regridder = xe.Regridder(
        model_data_t[["latitude", "longitude"]],
        obs_swath[["latitude", "longitude"]],
        "nearest_s2d",
        reuse_weights=False
    )

    no2_layer_swath = regridder(model_data_t["no2_layer"])
    pres_mid_swath = regridder(model_data_t["pres_pa_mid"])

    # Remove single model time dimension if present
    #no2_layer_swath = no2_layer_swath.squeeze()
    #pres_mid_swath = pres_mid_swath.squeeze()

    PressureLevel = obs_swath["PressureLevel"].values * 100.0 # hPa to Pa
    sat_presmid_pa = 0.5 * (PressureLevel[:, :, :-1] +
                            PressureLevel[:, :, 1:])

    shp_prior = obs_swath["NO2_ShapeFactor"].values
    ak_obs = obs_swath["AveragingKernel"].values

    ntime, nxtrack, nlay_sat = sat_presmid_pa.shape
    shp_mod = np.full_like(shp_prior, np.nan, dtype=np.float32)

    # Compute model shape factor per satellite pixel
    for i in range(ntime):
        for j in range(nxtrack):

            pm = pres_mid_swath[:, i, j].values
            nm = no2_layer_swath[:, i, j].values

            valid = (pm > 0) & (nm > 0) & np.isfinite(pm) & np.isfinite(nm)

            if valid.sum() < 2:
                continue

            try:
                spl = interpolate.splrep(np.log10(pm[valid]), nm[valid])
                interp_prof = interpolate.splev(
                    np.log10(sat_presmid_pa[i, j, :]),
                    spl
                )

                interp_prof = np.where(interp_prof > 0, interp_prof, np.nan)

                tot = np.nansum(interp_prof)

                if tot > 0:
                    shp_mod[i, j, :] = interp_prof / tot

            except Exception:
                continue

    # AK correction
    PTROP = 15000.0 # Pa = 150 hPa

    mask_total = np.ones_like(sat_presmid_pa, dtype=bool)
    mask_tropo = sat_presmid_pa > PTROP
    mask_strato = sat_presmid_pa <= PTROP

    term = (ak_obs - 1.0) * (shp_prior - shp_mod)

    ratio_total = 1.0 + np.nansum(term * mask_total, axis=2)
    ratio_tropo = 1.0 + np.nansum(term * mask_tropo, axis=2)
    ratio_strato = 1.0 + np.nansum(term * mask_strato, axis=2)

    if obs_no2_var == "no2_totalcolumn":
        ratio_use = ratio_total
    elif obs_no2_var == "no2_tropocolumn":
        ratio_use = ratio_tropo
    else:
        ratio_use = ratio_strato

    obs_raw = obs_swath[obs_no2_var].values
    obs_revised = np.where(
        ratio_use > 0,
        obs_raw * ratio_use,
        np.nan
    )

    return obs_raw, obs_revised

def omps_l2_no2_pairing_apriori_new(model_data,
                                    obs_swath,
                                    model_var_list,
                                    obs_no2_var="no2_totalcolumn"):

    import numpy as np
    import xarray as xr
    import gc

    valid_idx = _get_valid_model_time_indices(model_data, obs_swath)
    print(valid_idx)
    if len(valid_idx) == 0:
        raise ValueError("No model time within ±1 hour of swath UTC hour")

    model_hours_all = model_data["time_utc_hour"].values.astype(int)
    print(model_hours_all)
    # Use model grid
    lat_mod = model_data["latitude"].isel(time=0).values \
        if "time" in model_data["latitude"].dims else model_data["latitude"].values

    lon_mod = model_data["longitude"].isel(time=0).values \
        if "time" in model_data["longitude"].dims else model_data["longitude"].values

    grid_shape = lat_mod.shape

    # Final accumulated 2D grids
    sat_raw_sum = np.zeros(grid_shape, dtype=np.float32)
    sat_rev_sum = np.zeros(grid_shape, dtype=np.float32)
    model_sum = np.zeros(grid_shape, dtype=np.float32)
    n_obs_total = np.zeros(grid_shape, dtype=np.float32)

    for it in valid_idx:
        print(f"START {it}", flush=True)
        mh = int(model_hours_all[it])
        print(f"Processing model time index {it}, UTC hour {mh}",flush=True)

        obs_sub = _select_swath_for_model_hour(obs_swath, mh)

        if obs_sub is None:
            continue

        # Select one model time only; removes time dimension
        model_data_t = model_data.isel(time=it)
        print(f"BEFORE AK {it}", flush=True)
        # Pixel-level model-shape / AK adjustment
        obs_raw, obs_revised = _compute_pixel_level_satellite_adjustment(
            model_data_t,
            obs_sub,
            obs_no2_var
        )
        print(f"DONE AK {it}", flush=True)
        sat_lat = obs_sub["latitude"].values
        sat_lon = obs_sub["longitude"].values
        sat_lon = np.where(sat_lon < 0, sat_lon + 360, sat_lon)

        lon_mod_use = np.where(lon_mod < 0, lon_mod + 360, lon_mod)

        model_col = model_data_t[model_var_list[0]].values

        # Temporary 2D grid for this time window
        sat_raw_sum_t = np.zeros(grid_shape, dtype=np.float32)
        sat_rev_sum_t = np.zeros(grid_shape, dtype=np.float32)
        count_t = np.zeros(grid_shape, dtype=np.float32)

        for i in range(sat_lat.shape[0]):
            for j in range(sat_lat.shape[1]):

                if not np.isfinite(obs_revised[i, j]):
                    continue

                lat_val = sat_lat[i, j]
                lon_val = sat_lon[i, j]

                if not np.isfinite(lat_val) or not np.isfinite(lon_val):
                    continue

                # Nearest model grid cell
                iy = np.argmin(np.abs(lat_mod[:, 0] - lat_val))
                ix = np.argmin(np.abs(lon_mod_use[0, :] - lon_val))

                sat_raw_sum_t[iy, ix] += obs_raw[i, j]
                sat_rev_sum_t[iy, ix] += obs_revised[i, j]
                count_t[iy, ix] += 1
        print(f"DONE PIXEL LOOP {it}", flush=True)
        valid_cell = count_t > 0

        # Accumulate satellite values into final 2D grid
        sat_raw_sum[valid_cell] += sat_raw_sum_t[valid_cell]
        sat_rev_sum[valid_cell] += sat_rev_sum_t[valid_cell]

        # Model value repeated for cells that have satellite obs
        model_sum[valid_cell] += model_col[valid_cell] * count_t[valid_cell]

        n_obs_total[valid_cell] += count_t[valid_cell]
        print(f"DONE ACCUMULATION {it}",flush=True)
        del obs_sub, model_data_t
        gc.collect()

    # Final 2D satellite averages
    sat_raw_2d = np.where(n_obs_total > 0,sat_raw_sum / n_obs_total,np.nan)
    sat_rev_2d = np.where(n_obs_total > 0,sat_rev_sum / n_obs_total,np.nan)
    model_2d = np.where(n_obs_total > 0,model_sum / n_obs_total,np.nan)
    ds_out = xr.Dataset(
        {
            f"{obs_no2_var}_model": (
                ["grid_yt", "grid_xt"],
                model_2d.astype(np.float32)
            ),
            f"{obs_no2_var}_sat": (
                ["grid_yt", "grid_xt"],
                sat_rev_2d.astype(np.float32)
            ),
            f"{obs_no2_var}_sat_raw": (
                ["grid_yt", "grid_xt"],
                sat_raw_2d.astype(np.float32)
            ),
            "n_obs": (
                ["grid_yt", "grid_xt"],
                n_obs_total.astype(np.float32)
            ),
        },
        coords={
            "grid_yt": model_data["grid_yt"].values,
            "grid_xt": model_data["grid_xt"].values,
            "latitude": (
                ["grid_yt", "grid_xt"],
                lat_mod
            ),
            "longitude": (
                ["grid_yt", "grid_xt"],
                lon_mod
            ),
        },
        attrs={
            "pairing_method": "pixel_level_AK_adjustment_then_satellite_binning_to_model_grid",
            "time_selection": "model times within +/-1 hour of swath utc_hour",
            "final_dimension": "2D model grid",
        }
    )

    return ds_out

