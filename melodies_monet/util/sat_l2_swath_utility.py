# SPDX-License-Identifier: Apache-2.0
#

# read all swath data for the time range
# developed for TROPOMI Level2 NO2
#

import numpy as np
import xarray as xr
from datetime import datetime
import xesmf as xe

import logging
numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)

def trp_interp_swatogrd(obsobj, modobj,no2varname='no2'):

    """
    interpolate sat swath to model grid
    
    Parameters
    ------
    obsobj  : satellite swath data
    modobj  : model data (with no2 col calculated)
    
    Output
    ------
    no2_modgrid_avg: Regridded satellite data at model grids for all datetime

    """
    
    # model grids attributes
    nmodt, nz, ny, nx  = modobj[f'{no2varname}_col'].shape # time, z, y, x, no2 columns at molec cm^-2
    
    time   = [datetime.strptime(x,'%Y-%m-%d') for x in obsobj.keys()]
    nobstime  = len(list(obsobj.keys()))

    # daily averaged sat data at model grids
    no2_modgrid_avg=xr.Dataset(data_vars = {
        'nitrogendioxide_tropospheric_column':(["time", "x", "y"],
                                                np.full([nobstime, ny, nx], np.nan, dtype=np.float32)),
        f'{no2varname}trpcol':(["time", "x", "y"], np.full([nobstime, ny, nx], np.nan, dtype=np.float32))
            },
        coords = dict(
            time=time,
            longitude=(["x", "y"], modobj.coords['longitude'].values),
            latitude=(["x", "y"], modobj.coords['latitude'].values)),
        attrs=dict(description="daily tropomi data at model grids"),)

    for nd in range(nobstime):
        days = list(obsobj.keys())[nd]
        # --- model
        # get model no2 trop. columns at 13:00 - 14:00 localtime
        modobj_tm = modobj.sel(time=days.strfime('%Y-%d-%m'))
        
        # intermediate need: model NO2 partial columns for day
        # no2col_satm = np.nanmean(modobj_tm['no2col'].values, axis = 0)
        
        # --- TROPOMI
        # number of swath
        nswath = len(obsobj[days])

        # intermediate array for all swaths
        no2_modgrid_all = np.zeros([ny, nx, nswath], dtype=np.float64)
        tropopause_intermediate = np.zeros([ny, nx, nswath], dtype=np.float64)
        for ns in range(nswath):
            satlon = obsobj[days][ns]['lon']
            satlat = obsobj[days][ns]['lat']
            satno2 = obsobj[days][ns]['nitrogendioxide_tropospheric_column']
            tm5_tropopause = obsobj[days][ns]['troppres']
            # regridding from swath grid to model grids
            grid_in = {'lon':satlon.values, 'lat':satlat.values}

            regridder = xe.Regridder(grid_in, no2_modgrid_avg[['lat','lon']],'bilinear',ignore_degenerate=True,reuse_weights=False)
            
            # regridded no2 trop. columns
            no2_modgrid = regridder(satno2) # , keep_attrs=True
            tpause_modgrid = regridder(tm5_tropopause)
            print('Done with TROPOMI regridding', days, ns)

            #regridder.destroy()
            del regridder
 
            no2_modgrid_all[:,:,ns] = no2_modgrid
            tropopause_intermediate[:,:,ns] = tpause_modgrid
            print(' no2 satellite:', np.nanmin(no2_modgrid), np.nanmax(no2_modgrid))

        # daily averaged no2 trop. columns at model grids
        no2_modgrid_avg['nitrogendioxide_tropospheric_column'][nd,:,:] = np.nanmean(np.where(no2_modgrid_all > 0.0, no2_modgrid_all, np.nan), axis=2)
        # daily averaged TM5 tropopause at model grid
        tm5_tpause = np.nanmean(tropopause_intermediate,axis=2)
        
        # sum model to TM5 tropopause
        no2_modgrid_avg[f'{no2varname}trpcol'][nd, :,:] = modobj_tm[f'{no2varname}_col'].where(modobj_tm['pres_pa_mid'] >= tm5_tpause).sum(dim='z').values.squeeze()
    del(modobj)
    del(obsobj)

    return no2_modgrid_avg


def trp_interp_swatogrd_ak(obsobj, modobj,no2varname='no2'):

    """
    interpolate sat swath to model grid applied with averaging kernel
    
    Parameters
    ------
    obsobj  : satellite swath data
    modobj  : model data (with no2 col calculated)
    
    Output
    ------
    no2_modgrid_avg: Regridded satellite data at model grids for all datetime

    """

    # model grids attributes
    nmodt, nz, ny, nx  = modobj[f'{no2varname}_col'].shape # time, z, y, x, no2 columns at molec cm^-2
    
    time   = [datetime.strptime(x,'%Y-%m-%d') for x in obsobj.keys()]
    nobstime  = len(list(obsobj.keys()))

    # daily averaged sat data at model grids
    no2_modgrid_avg=xr.Dataset(data_vars = {
        'nitrogendioxide_tropospheric_column':(["time", "x", "y"],
                                                np.full([nobstime, ny, nx], np.nan, dtype=np.float32)),
        f'{no2varname}trpcol':(["time", "x", "y"],np.full([nobstime, ny, nx], np.nan, dtype=np.float32))
            },
        coords = dict(
            time=time,
            longitude=(["x", "y"], modobj.coords['longitude'].values),
            latitude=(["x", "y"], modobj.coords['latitude'].values)),
        attrs=dict(description="daily tropomi data at model grids"),)

    # tmpvalue = np.zeros([ny, nx], dtype = np.float64)

    # loop over all days
    for nd in range(nobstime):

        days = time[nd].strftime('%Y-%m-%d')
        # --- model ---
        # get model no2 trop. columns at 13:00 - 14:00 localtime
        try:
            modobj_tm = modobj.sel(time=days)
        except KeyError:
            print(days)
            print('Satellite data was outside available model times')
            continue
        #modobj_tm = modobj.sel(time=days)
        # no2col_satm = modobj_tm[f'{no2varname}_col'].mean(dim='time')
              
        # --- tropomi ---
        # number of swath
        nswath = len(obsobj[days])

        # array for all swaths
        no2_modgrid_all = np.zeros([ny, nx, nswath], dtype=np.float32)
        tropopause_intermediate = np.zeros([ny, nx, nswath], dtype=np.float32)
        for ns in range(nswath):
            working_swath = obsobj[days][ns]     

            grid_sat = {'lon':working_swath['lon'].values, 'lat':working_swath['lat'].values}
            grid_mod= {'lon':modobj.coords['longitude'].values, 'lat':modobj.coords['latitude'].values}


            nysat, nxsat, nzsat = working_swath['averaging_kernel'].shape

            # regridding from model grid to sat grid
            regridder_ms = xe.Regridder(grid_mod, grid_sat,'bilinear',ignore_degenerate=True,reuse_weights=False)
            
            # force model data to put z dimension last for pressure and no2 partial columns
            mod_pres_no2 = modobj_tm[['pres_pa_mid',f'{no2varname}_col']].mean(dim='time')#.transpose('y','x','z')
            #print(mod_pres_no2['no2col'].shape)
            # regridding for model pressure, and no2 vertical columns
            mod_rgd_sat = regridder_ms(mod_pres_no2)
            mod_rgd_sat = mod_rgd_sat.transpose('y','x','z')
            # convert from aks to trop.aks
            working_swath['averaging_kernel'] = working_swath['averaging_kernel'] * working_swath['air_mass_factor_total'] / working_swath['air_mass_factor_troposphere']
            # calculate the revised tamf_mod, and ratio = tamf_mod / tamf_org
            ratio = cal_amf_wrfchem(working_swath['averaging_kernel'], mod_rgd_sat['pres_pa_mid'].values, working_swath['preslev'], working_swath['troppres'], mod_rgd_sat[f'{no2varname}_col'].values,
                                    working_swath['air_mass_factor_troposphere'], grid_sat['lon'], grid_sat['lat'], grid_mod['lon'], grid_mod['lat'])

            # averaing kernel applied done
            satno2 = working_swath['nitrogendioxide_tropospheric_column'] * ratio 
            
            tm5_tropopause = working_swath['troppres'] 

            # regridding from swath grid to model grids
            regridder = xe.Regridder(grid_sat, grid_mod,'bilinear',ignore_degenerate=True,reuse_weights=False,unmapped_to_nan=True)

            # regridded no2 trop. columns
            no2_modgrid = regridder(satno2, keep_attrs=True)
            tpause_modgrid = regridder(tm5_tropopause, keep_attrs=True)
            no2_modgrid_all[:,:,ns] = no2_modgrid[:,:]
            tropopause_intermediate[:,:,ns] = tpause_modgrid

        # daily averaged TM5 tropopause at model grid
        tm5_tpause = np.nanmean(np.where(no2_modgrid_all > 0.0,tropopause_intermediate,np.nan),axis=2)

        # sum model to TM5 tropopause
        no2_modgrid_avg[f'{no2varname}trpcol'][nd, :,:] = modobj_tm[f'{no2varname}_col'].where(modobj_tm['pres_pa_mid'] >= tm5_tpause).sum(dim='z').values.squeeze()

        # daily averaged no2 trop. columns at model grids
        no2_modgrid_avg['nitrogendioxide_tropospheric_column'][nd,:,:] = np.nanmean(np.where(no2_modgrid_all > 0.0, no2_modgrid_all, np.nan), axis=2)

    return no2_modgrid_avg


def cal_amf_wrfchem(scatw, wrfpreslayer, tpreslev, troppres, wrfno2layer_molec, tamf_org, satlon, satlat, modlon, modlat):
    from scipy import interpolate

    nsaty, nsatx, nz    = wrfpreslayer.shape
    nsatz, nsaty, nsatx = tpreslev.shape # mli, update to new dimension


    nume             = np.zeros([nsaty, nsatx], dtype=np.float32)
    deno             = np.zeros([nsaty, nsatx], dtype=np.float32)
    amf_wrfchem      = np.zeros([nsaty, nsatx], dtype=np.float32)
    amf_wrfchem[:,:] = np.nan
    wrfavk           = np.zeros([nsaty, nsatx, nz], dtype = np.float32)
    wrfavk[:,:,:]    = np.nan
    wrfavk_scl       = np.zeros([nsaty, nsatx], dtype=np.float32) 
    preminus         = np.zeros([nsaty, nsatx], dtype=np.float32)
    wrfpreslayer_slc = np.zeros([nsaty, nsatx], dtype=np.float32)
    tmpvalue_sat     = np.zeros([nsaty, nsatx], dtype=np.float32)
    tmpvalue_mod     = np.zeros([nsaty, nsatx], dtype=np.float32)
    
    
    # set the surface pressure to wrf one
    tpreslev[0,:,:] = wrfpreslayer[:,:,0] 

    # relationship between pressure to avk
    tpreslev = tpreslev.values 
    scatw    = scatw.values
    wrfpreslayer = np.where((wrfpreslayer <=0.0), np.nan, wrfpreslayer)

    # shrink the satellite domain to WRF
    lb = np.where( (satlon >= np.nanmin(modlon)) & (satlon <= np.nanmax(modlon)) 
        & (satlat >= np.nanmin(modlat)) & (satlat <= np.nanmax(modlat)))

    vertical_pres = []
    vertical_scatw = []
    vertical_wrfp = []
    
    if len(lb[0]) == 0:
        print('Caution: There are no observations within the model domain')
    for llb in range(len(lb[0])):
        yy = lb[0][llb]
        xx = lb[1][llb]
        vertical_pres = tpreslev[:,yy,xx] # mli, update to new dimension
        vertical_scatw = scatw[yy,xx,:]
        vertical_wrfp = wrfpreslayer[yy,xx,:]
        f = interpolate.interp1d(np.log10(vertical_pres[:]),vertical_scatw[:], fill_value="extrapolate")# relationship between pressure to avk
        wrfavk[yy,xx,:] = f(np.log10(vertical_wrfp[:])) #wrf-chem averaging kernel

    for l in range(nz-1):  # noqa: E741
        # check if it's within tropopause
        preminus[:,:]         = wrfpreslayer[:,:,l] - troppres[:,:]

        # wrfpressure and wrfavk
        wrfpreslayer_slc[:,:] = wrfpreslayer[:,:,l]
        wrfavk_scl[:,:]       = wrfavk[:,:,l]

        ind_ak = np.where(np.isinf(wrfavk_scl) | (wrfavk_scl <= 0.0))
        # use the upper level ak 
        if (ind_ak[0].size >= 1):
            tmpvalue_sat[:,:]  = wrfavk[:,:,l+1]
            wrfavk_scl[ind_ak] = tmpvalue_sat[ind_ak]

        ind = np.where(preminus >= 0.0)
        # within tropopause
        if (ind[0].size >= 1):
            # select grids that this level is within tropopause
            tmpvalue_mod[:,:]  = wrfno2layer_molec[:,:,l]
            nume[ind] += wrfavk_scl[ind]*tmpvalue_mod[ind]
            deno[ind] += tmpvalue_mod[ind]
        else:
            break
            
    # tropospheric amf calculated based on model profile and TROPOMI averaging kernel
    amf_wrfchem = nume / deno * tamf_org

    # ratio
    ratio = tamf_org / amf_wrfchem 

    # exclude nan
    ratio = np.where(np.isnan(ratio), 1.0, ratio)

    print('Done with Averaging Kernel revision,', 'factor min:',np.nanmin(ratio), 'max:',np.nanmax(ratio)) 

    return ratio 

