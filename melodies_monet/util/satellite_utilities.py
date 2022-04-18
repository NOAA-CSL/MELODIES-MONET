# File started by Maggie Bruckner. 
# Contains satellite specific pairing operators
import xesmf as xe
import numpy as np
from datetime import datetime,timedelta

def omps_nm_pairing(model_data,obs_data):
    'Pairs UFS-RAQMS ozone mixing ratio with OMPS nadir mapper retrievals. Calculates column without applying apriori'
    import xarray as xr

    import pandas as pd
    
    print('pairing without applying averaging kernel')

    
    du_fac = 1.0e4*6.023e23/28.97/9.8/2.687e19 # conversion factor; moves model from ppv to dobson
    
    # define satellite lat/lon grid as dataset
    ds_out = xr.Dataset({'lat': (['x','y'], obs_data['latitude'].values),
                         'lon': (['x','y'], obs_data['longitude'].values),
                        }
                       )
    ds_mod = xr.Dataset({'lat':(['x','y'],model_data['latitude'].values),
                         'lon':(['x','y'],model_data['longitude'].values),})
    # regrid spatially (model lat/lon to satellite swath lat/lon)
    # dimensions of new variables will be (time, z, satellite_x, satellite_y)
    regridr = xe.Regridder(model_data,ds_out,'bilinear') # standard bilinear spatial regrid. 
    regrid_oz = regridr(model_data['o3vmr'])
    regrid_dp = regridr(model_data['dpm'])
    nf,nz,nx,ny = regrid_oz.shape

    # calculate ozone column, no averaging kernel or apriori applied.
    col = np.nansum(du_fac*regrid_dp*regrid_oz,axis=1) # new dimensions will be (time, satellite_x, satellite_y)
   
    nf,nx,ny = col.shape
    
    ## initialize dataset holder for final dataset
    oz = np.zeros_like(obs_data.ozone_column.values)
    
    
    ## loop over model time steps
    for f in range(nf):
        tindex = np.where(np.abs(obs_data.time - model_data.time[f]) <= (model_data.time[1]-model_data.time[0]))[0]
        #tfac1 = 1-np.abs(model_data.time[f] - obs_data.time[tindex])/(model_data.time[1]-model_data.time[0])

        #sum_tf[tindex] += tfac1
        # fixes for observations before/after model time range.
        if f == (nf-1):
            print('last')
            tindex = np.where((obs_data.time >= model_data.time[f]))[0]
            oz[tindex,:] = col[f][tindex,:]#.values
            #sum_tf[tindex] += 1
            
            tind_2 = np.where((obs_data.time < model_data.time[f]) & 
                              (np.abs(obs_data.time - model_data.time[f]) <= (model_data.time[1]-model_data.time[0])))[0]
            tfac1 = 1-(np.abs(model_data.time[f] - obs_data.time[tind_2])/(model_data.time[1]-model_data.time[0]))
            
            oz[tind_2,:] += np.expand_dims(tfac1.values,axis=1)*col[f][tind_2,:]
            #sum_tf[tind_2] += tfac1
        elif f == (0):
            print('first')
            tindex = np.where((obs_data.time <= model_data.time[f]))[0]
            oz[tindex,:] = col[f][tindex,:]#.values
            #sum_tf[tindex] += 1
            tind_2 = np.where((obs_data.time > model_data.time[f]) & 
                              (np.abs(obs_data.time - model_data.time[f]) <= (model_data.time[1]-model_data.time[0])))[0]
            tfac1 = 1-(np.abs(model_data.time[f] - obs_data.time[tind_2])/(model_data.time[1]-model_data.time[0]))
            
            oz[tind_2,:] += np.expand_dims(tfac1.values,axis=1)*col[f][tind_2,:]
            #sum_tf[tind_2] += tfac1
        else:
            print('not 1 or last')
            tindex = np.where(np.abs(obs_data.time - model_data.time[f]) <= (model_data.time[1]-model_data.time[0]))[0]
            
            tfac1 = 1-(np.abs(model_data.time[f] - obs_data.time[tindex])/(model_data.time[1]-model_data.time[0]))
            #sum_tf[tindex] += tfac1
            oz[tindex,:] += np.expand_dims(tfac1.values,axis=1)*col[f][tindex,:]#.values

    ds = xr.Dataset({'o3vmr': (['x','y'],oz),
                     'ozone_column':(['x','y'],obs_data.ozone_column.values)
                               },
                    coords={
                        'longitude':(['x','y'],obs_data['longitude'].values),
                        'latitude':(['x','y'],obs_data['latitude'].values),
                        'time':(['x'],obs_data.time.values),
                    })    

    return ds

def omps_nm_pairing_apriori(model_data,obs_data):
    'Pairs UFS-RAQMS data with OMPS nm. Applies satellite apriori column to model observations.'

    import xarray as xr

    import pandas as pd
    du_fac = 1.0e4*6.023e23/28.97/9.8/2.687e19 # conversion factor; moves model from ppv to dobson
    
    print('pairing with averaging kernel application')
    # define satellite lat/lon grid as dataset
    ds_out = xr.Dataset({'lat': (['x','y'], obs_data['latitude'].values),
                         'lon': (['x','y'], obs_data['longitude'].values),
                        }
                       )
    # regrid spatially (model lat/lon to satellite swath lat/lon)
    # dimensions of new variables will be (time, z, satellite_x, satellite_y)
    regridr = xe.Regridder(model_data,ds_out,'bilinear') # standard bilinear spatial regrid. 
    regrid_oz = regridr(model_data['o3vmr'])
    regrid_p = regridr(model_data['pdash']) # this one should be pressure variable (for the interpolation).
    sfp = regridr(model_data['sfcp'])
    # Interpolate model data in time. Final dataset will be (satellite_x, satellite_y) dimensions at satellite times.
    nf,nz_m,nx,ny = regrid_oz.shape

    ## initialize intermediates for use in calcluating column
    pressure_temp = np.zeros((nz_m,nx,ny))
    ozone_temp = np.zeros((nz_m,nx,ny))
    sfc = np.zeros((nx,ny))
    ## loop over model time steps
    for f in range(nf):
        tindex = np.where(np.abs(obs_data.time - model_data.time[f]) <= (model_data.time[1]-model_data.time[0]))[0]
        # fixes for observations before/after model time range.
        if f == (nf-1):
            tindex = np.where((obs_data.time >= model_data.time[f]))[0]
            ozone_temp[:,tindex,:] = regrid_oz[f,:,tindex,:].values
            pressure_temp[:,tindex,:] = regrid_p[f,:,tindex,:].values
            sfc[tindex,:] = sfp[f,tindex,:].values 
            tind_2 = np.where((obs_data.time < model_data.time[f]) & 
                              (np.abs(obs_data.time - model_data.time[f]) <= (model_data.time[1]-model_data.time[0])))[0]
            tfac1 = 1-(np.abs(model_data.time[f] - obs_data.time[tind_2])/(model_data.time[1]-model_data.time[0]))
            
            ozone_temp[:,tind_2,:] += np.expand_dims(tfac1.values,axis=1)*regrid_oz[f,:,tind_2,:].values
            pressure_temp[:,tind_2,:] += np.expand_dims(tfac1.values,axis=1)*regrid_p[f,:,tind_2,:].values
            sfc[tind_2,:] += np.expand_dims(tfac1.values,axis=1)*sfp[f,tind_2,:].values
        elif f == 0:
            tindex = np.where((obs_data.time <= model_data.time[f]))[0]
            ozone_temp[:,tindex,:] = regrid_oz[f,:,tindex,:].values
            pressure_temp[:,tindex,:] = regrid_p[f,:,tindex,:].values
            sfc[tindex,:] = sfp[f,tindex,:].values 
            tind_2 = np.where((obs_data.time > model_data.time[f]) & 
                              (np.abs(obs_data.time - model_data.time[f]) <= (model_data.time[1]-model_data.time[0])))[0]
            tfac1 = 1-(np.abs(model_data.time[f] - obs_data.time[tind_2])/(model_data.time[1]-model_data.time[0]))
            ozone_temp[:,tind_2,:] += np.expand_dims(tfac1.values,axis=1)*regrid_oz[f,:,tind_2,:].values
            pressure_temp[:,tind_2,:] += np.expand_dims(tfac1.values,axis=1)*regrid_p[f,:,tind_2,:].values
            sfc[tind_2,:] += np.expand_dims(tfac1.values,axis=1)*sfp[f,tind_2,:].values
        else:
            tfac1 = 1-(np.abs(model_data.time[f] - obs_data.time[tindex])/(model_data.time[1]-model_data.time[0]))
            ozone_temp[:,tindex,:] += np.expand_dims(tfac1.values,axis=1)*regrid_oz[f,:,tindex,:].values
            pressure_temp[:,tindex,:] += np.expand_dims(tfac1.values,axis=1)*regrid_p[f,:,tindex,:].values
            sfc[tindex,:] += np.expand_dims(tfac1.values,axis=1)*sfp[f,tindex,:].values
    # Interpolate model data to satellite pressure levels
    from wrf import interplevel
    ozone_satp = interplevel(ozone_temp,pressure_temp/100.,obs_data.pressure,missing=np.nan)
    ozone_satp = ozone_satp.values
    
    ozone_satp[np.isnan(ozone_satp)] = 0
    oz = np.zeros_like(obs_data.ozone_column.values)
    # flip z dimension. Had to flip obs pressure in interpolation for reasons (to make it work)
    #ozone_satp = ozone_satp[::-1]
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

    ds = xr.Dataset({'o3vmr': (['x','y'],oz),
                     'ozone_column':(['x','y'],obs_data.ozone_column.values)
                               },
                    coords={
                        'longitude':(['x','y'],obs_data['longitude'].values),
                        'latitude':(['x','y'],obs_data['latitude'].values),
                        'time':(['x'],obs_data.time.values),
                    })
    return ds
