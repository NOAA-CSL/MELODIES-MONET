# read all swath data for the time range
# developed for TROPOMI Level2 NO2

import xesmf as xe
import numpy as np
import xarray as xr
import pandas as pd
from datetime import datetime
from collections import OrderedDict

import logging
numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)

def trp_interp_swatogrd(obsobj, modobj):

    from datetime import datetime


    """
    interpolate sat swath to model grid
    
    Parameters
    ------
    obsobj  : satellite swath data
    modobj  : model data
    
    Output
    ------
    no2_modgrid_avg: Regridded satellite data at model grids for all datetime

    """
    
    # daily averaged sate data at model grids
    no2_modgrid_avg=xr.Dataset()

    # revised modobj with more variables
    modobj = model_columns(modobj)
    print(' mobel obj: ', modobj)
    print(' trpomi obj: ', obsobj)

    # model grids attributes
    nt, nz, ny, nx  = modobj['no2col'].shape # time, z, y, x, no2 columns at molec cm^-2
    modlat = modobj.coords['latitude']
    modlon = modobj.coords['longitude']

    time   = [ datetime.strptime(x,'%Y-%m-%d') for x in obsobj.keys()]
    ntime  = len(list(obsobj.keys()))

    no2_nt = np.zeros([ntime, ny, nx], dtype=np.float32)
    no2_mod = np.zeros([ntime, ny, nx], dtype=np.float32)

    for nd in range(ntime):
        datetime = list(obsobj.keys())[nd]
        # get model no2 trop. columns
        modobj_tm = modobj.sel(time=((modobj.time.dt.strftime("%Y-%m-%d") == datetime) & (modobj.time.dt.hour >= 13)
            & (modobj.time.dt.hour <= 14)))['no2col'].values
        no2col_satm = np.nanmean(modobj_tm, axis = 0)
        
        # sum up tropopause
        no2_mod[nd, :,:] = np.nansum(no2col_satm[0:49,:,:], axis=0)

    del(modobj)


    # loop over all datetime
    for nd in range(ntime):

        datetime = list(obsobj.keys())[nd]

        # number of swath
        nswath = len(obsobj[datetime])

        # array for all swaths
        no2_modgrid_all = np.zeros([ny, nx, nswath], dtype=np.float32)

        for ns in range(nswath):
            satlon = obsobj[datetime][ns]['lon']
            satlat = obsobj[datetime][ns]['lat']
            satno2 = obsobj[datetime][ns]['nitrogendioxide_tropospheric_column']

            # regridding from swath grid to model grids
            grid_in = {'lon':satlon.values, 'lat':satlat.values}
            grid_out= {'lon':modlon.values, 'lat':modlat.values}

            regridder = xe.Regridder(grid_in, grid_out,'bilinear',ignore_degenerate=True,reuse_weights=False)
            
            # regridded no2 trop. columns
            no2_modgrid = regridder(satno2) # , keep_attrs=True
            print('Done with TROPOMI regridding', datetime, ns)

            #regridder.destroy()
            del(regridder)
            regridder = None
 
            no2_modgrid_all[:,:,ns] = no2_modgrid[:,:]
            print(' no2 satellite:', np.nanmin(no2_modgrid), np.nanmax(no2_modgrid))

        # daily averaged no2 trop. columns at model grids
        no2_nt[nd,:,:] = np.nanmean(np.where(no2_modgrid_all > 0.0, no2_modgrid_all, np.nan), axis=2)

    del(obsobj)
 
    no2_modgrid_avg = xr.Dataset(
        data_vars = dict(
            no2sat=(["time", "x", "y"],no2_nt),
            latitude=(["x", "y"],modlat.values),
            longitude=(["x", "y"],modlon.values),
            no2mod = (["time", "x", "y"],no2_mod )
            ),
        coords = dict(
            time=time,
            lon=(["x", "y"], modlon.values),
            lat=(["x", "y"], modlat.values)),
        attrs=dict(description="daily tropomi data at model grids"),
        )

    return no2_modgrid_avg


def trp_interp_swatogrd_ak(obsobj, modobj):

    from datetime import datetime

    """
    interpolate sat swath to model grid applied with averaging kernel
    
    Parameters
    ------
    obsobj  : satellite swath data
    modobj  : model data
    
    Output
    ------
    no2_modgrid_avg: Regridded satellite data at model grids for all datetime

    """
    
    # daily averaged sate data at model grids
    no2_modgrid_avg=xr.Dataset()

    # revised modobj with more variables
    modobj = model_columns(modobj)

    # model grids attributes
    nt, nz, ny, nx  = modobj['no2'].shape
    modlat = modobj.coords['latitude']
    modlon = modobj.coords['longitude']

    tmpvalue = np.zeros([ny, nx], dtype = np.float)

    time   = [ datetime.strptime(x,'%Y-%m-%d') for x in obsobj.keys()]
    ntime  = len(list(obsobj.keys()))

    no2_nt = np.zeros([ntime, ny, nx], dtype=np.float32)
    no2_mod = np.zeros([ntime, ny, nx], dtype=np.float32)


    # loop over all datetime
    for nd in range(ntime):

        datetime = list(obsobj.keys())[nd]

        # --- model ---
        # get model no2 trop. columns
        modobj_tm = modobj.sel(time=((modobj.time.dt.strftime("%Y-%m-%d") == datetime) & (modobj.time.dt.hour >= 13)
            & (modobj.time.dt.hour <= 14)))
        no2col_satm = np.nanmean(modobj_tm['no2col'].values, axis = 0)
        
        # sum up tropopause, needs to be revised to tropopause
        no2_mod[nd, :,:] = np.nansum(no2col_satm[0:49,:,:], axis=0)

        # --- tropomi ---
        # number of swath
        nswath = len(obsobj[datetime])

        # array for all swaths
        no2_modgrid_all = np.zeros([ny, nx, nswath], dtype=np.float32)

        for ns in range(nswath):
            satlon = obsobj[datetime][ns]['lon']
            satlat = obsobj[datetime][ns]['lat']
            satno2 = obsobj[datetime][ns]['nitrogendioxide_tropospheric_column']            

            grid_sat = {'lon':satlon.values, 'lat':satlat.values}
            grid_mod= {'lon':modlon.values, 'lat':modlat.values}

            # -- applying averaging kernel ---
            # trop. amf in standard product
            tamf_org   = obsobj[datetime][ns]['air_mass_factor_troposphere']
            amf_total  = obsobj[datetime][ns]['air_mass_factor_total']
            troppres   = obsobj[datetime][ns]['troppres'] # TM5 tropopause pressure, Pa 
            tpreslev   = obsobj[datetime][ns]['preslev']
            scatwts    = obsobj[datetime][ns]['averaging_kernel']

            nysat, nxsat, nzsat = scatwts.shape

            # regridding from model grid to sat grid
            regridder_ms = xe.Regridder(grid_mod, grid_sat,'bilinear',ignore_degenerate=True,reuse_weights=False)

            # regridding for model pressure, and no2 vertical colums
            wrfpres        =  np.zeros([nysat, nxsat, nz], dtype = np.float32)
            wrfpres[:,:,:] =  np.nan 
            wrfno2         =  np.zeros([nysat, nxsat, nz], dtype = np.float32)
            wrfno2[:,:,:]  =  np.nan
            modvalue_pb2   =  np.nanmean(modobj_tm['PB2'].values, axis = 0)
            modvalue_no2   =  np.nanmean(modobj_tm['no2col'].values, axis = 0)

            for l in range(nz):
                tmpvalue[:,:]  = modvalue_pb2[l,:,:]         
                wrfpres[:,:,l] = regridder_ms(tmpvalue)
                tmpvalue[:,:]  = modvalue_no2[l,:,:]
                wrfno2[:,:,l]  = regridder_ms(tmpvalue)            


            print('checking wrfpres', np.nanmin(wrfpres),np.nanmax(wrfpres) )
            print('checking wrfpres', np.nanmin(modvalue_pb2),np.nanmax(modvalue_pb2) )

            # convert from aks to trop.aks
            for l in range(nzsat):
                scatwts[:,:,l] = scatwts[:,:,l] * amf_total[:,:] / tamf_org[:,:]

            # calculate the revised tamf_mod, and ratio = tamf_mod / tamf_org
            ratio = cal_amf_wrfchem(scatwts, wrfpres, tpreslev, troppres, wrfno2, tamf_org, satlon.values, satlat.values, modlon, modlat)

            # averaing kernel applied done
            satno2 = satno2 * ratio 

            # regridding from swath grid to model grids
            regridder = xe.Regridder(grid_sat, grid_mod,'bilinear',ignore_degenerate=True,reuse_weights=False)

            # regridded no2 trop. columns
            no2_modgrid = regridder(satno2, keep_attrs=True)
            no2_modgrid_all[:,:,ns] = no2_modgrid[:,:]

        # daily averaged no2 trop. columns at model grids
        no2_nt[nd,:,:] = np.nanmean(np.where(no2_modgrid_all > 0.0, no2_modgrid_all, np.nan), axis=2)


    no2_modgrid_avg = xr.Dataset(
        data_vars = dict(
            no2sat=(["time", "x", "y"],no2_nt),
            latitude=(["x", "y"],modlat.values),
            longitude=(["x", "y"],modlon.values),
            no2mod = (["time", "x", "y"],no2_mod)
            #ratio=(["time", "x", "y"],ratio_nt)
            ),
        coords = dict(
            time=time,
            lon=(["x", "y"], modlon.values),
            lat=(["x", "y"], modlat.values)),
        attrs=dict(description="daily tropomi data at model grids"),
        )
    return no2_modgrid_avg


def cal_amf_wrfchem(scatw, wrfpreslayer, tpreslev, troppres, wrfno2layer_molec, tamf_org, satlon, satlat, modlon, modlat):
    from scipy import interpolate

    nsaty, nsatx, nz    = wrfpreslayer.shape
    nsaty, nsatx, nsatz = tpreslev.shape

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
    
    
    # set the surface pressure to wrf one
    tpreslev[:,:,0] = wrfpreslayer[:,:,0] 

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

    for llb in range(len(lb[0])):
        yy = lb[0][llb]
        xx = lb[1][llb]
        vertical_pres = tpreslev[yy,xx,:]
        vertical_scatw = scatw[yy,xx,:]
        vertical_wrfp = wrfpreslayer[yy,xx,:]

        f = interpolate.interp1d(np.log10(vertical_pres[:]),vertical_scatw[:], fill_value="extrapolate")# relationship between pressure to avk
        wrfavk[yy,xx,:] = f(np.log10(vertical_wrfp[:])) #wrf-chem averaging kernel

        #if vertical_pres[0] != 0.0:
            #print(yy, xx, satlon[yy,xx], satlat[yy,xx])
            #print('pres', vertical_pres, len(vertical_pres))
            #print('scatw', vertical_scatw, len(vertical_scatw))
            #print('wrfp', vertical_wrfp, len(vertical_wrfp))
            #print('wrfavk', wrfavk[yy,xx,:])


    for l in range(nz-1):
        # check if it's within tropopause
        preminus[:,:]         = wrfpreslayer[:,:,l] - troppres[:,:]

        # wrfpressure and wrfavk
        wrfpreslayer_slc[:,:] = wrfpreslayer[:,:,l]
        wrfavk_scl[:,:]       = wrfavk[:,:,l]

        ind_ak = np.where((np.isinf(wrfavk_scl) == True) | (wrfavk_scl <= 0.0))
        # use the upper level ak 
        if (ind_ak[0].size >= 1):
            tmpvalue_sat[:,:]  = wrfavk[:,:,l+1]
            wrfavk_scl[ind_ak] = tmpvalue_sat[ind_ak]

        ind = np.where(preminus >= 0.0)
        # within tropopause
        if (ind[0].size >= 1):
            nume[:,:] += wrfavk_scl[:,:]*wrfno2layer_molec[:,:,l]
            deno[:,:] += wrfno2layer_molec[:,:,l]
        else:
            break
            
    # tropospheric amf calculated based on model profile and TROPOMI averaging kernel
    amf_wrfchem = nume / deno * tamf_org

    # ratio
    ratio = tamf_org / amf_wrfchem 

    # exclude nan
    ratio = np.where((np.isnan(ratio) == True), 1.0, ratio)

    print('Done with Averaging Kernel revision,', 'factor min:',np.nanmin(ratio), 'max:',np.nanmax(ratio)) 

    return ratio 

def model_columns(modobj):
    from datetime import tzinfo, timedelta

    # calculate the no2 tropospheric vertical columns and pressure from wrf-chem
    # local time averaged between 13:00 - 14:00
    no2    = modobj['no2']
    ph     = modobj['PH']
    phb    = modobj['PHB']
    pb     = modobj['PB']
    tdata  = modobj['T']
    pdata  = modobj['P']
    time   = modobj.coords['time']
    modlon = modobj.coords['longitude']

    # presure: base state + PB (KSMP)
    nt, nz, ny, nx = no2.shape
    pb2  = np.zeros([nt, nz, ny, nx],dtype=np.float)
    pb2  = pdata + pb

    # convert the perturbation potential temperature (from 300K reference) to temp
    tb = np.zeros([nt, nz, ny, nx],dtype=np.float)
    tb =(300.0+tdata)*((pb2/1.0e5)**0.286)

    modobj['PB2'] = xr.DataArray(pb2)
    modobj['TB']  = xr.DataArray(tb) 

    no2col = np.zeros([nt, nz, ny, nx], dtype = np.float32)
    value = np.zeros([nt, ny, nx], dtype = np.float32)

    # average between 13:00 and 14:00 localtime
    localtm   = np.zeros([nt,ny,nx], dtype='datetime64[s]')
    tdlist    = np.zeros([ny], dtype=np.int32)
    tdlt      = np.zeros([ny, nx], dtype='timedelta64[ms]')

    for xx in range(nx):
         tdlist[:] = np.array(modlon[:,xx]/15.0).astype(int)
         tdlt[:,xx] = pd.TimedeltaIndex(tdlist, 'h')

    for tt in range(nt):
        localtm[tt,:,:] = pd.Timestamp(time.values[tt]) + tdlt[:,:]

    modobj['localtime'] = xr.DataArray(localtm, dims=["time", "y", "x"])

    # convert to ppm
    no2 = no2 / 1000.0
    for vl in range(nz):
        ad = pb2[:,vl,:,:] * (28.97e-3)/(8.314*tb[:,vl,:,:])
        zh = ((ph[:,vl+1,:,:] + phb[:,vl+1,:,:]) - (ph[:,vl,:,:]+phb[:,vl,:,:]))/9.81
        value[:,:,:]= no2[:,vl,:,:]*zh[:,:,:]*6.022e23/(28.97e-3)*1e-10*ad[:,:,:] # timex y x x
        no2col[:,vl,:,:] = value[:,:,:]

    modobj['no2col'] = xr.DataArray(no2col,dims=["time", "z", "y","x"])

    return modobj
