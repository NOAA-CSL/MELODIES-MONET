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

logger = logging.getLogger(__name__)

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
        
        # sum up tropopause
        if 'pres_pa_trop' in list(modobj.keys()):
            no2_modgrid_avg[f'{no2varname}trpcol'][nd, :,:] = modobj_tm[f'{no2varname}_col'].where(modobj_tm['pres_pa_mid'] >= modobj_tm['pres_pa_trop']).sum(dim='z').values.squeeze()

        else:
            print('Caution: model tropospheric NO2 column was calculated assuming the model top is the tropopause')
            no2_modgrid_avg[f'{no2varname}trpcol'][nd, :,:] = modobj_tm[f'{no2varname}_col'].sum(dim='z').values.squeeze()
            
        # --- TROPOMI
        # number of swath
        nswath = len(obsobj[days])

        # intermediate array for all swaths
        no2_modgrid_all = np.zeros([ny, nx, nswath], dtype=np.float64)

        for ns in range(nswath):
            satlon = obsobj[days][ns]['lon']
            satlat = obsobj[days][ns]['lat']
            satno2 = obsobj[days][ns]['nitrogendioxide_tropospheric_column']

            # regridding from swath grid to model grids
            grid_in = {'lon':satlon.values, 'lat':satlat.values}

            regridder = xe.Regridder(grid_in, no2_modgrid_avg[['lat','lon']],'bilinear',ignore_degenerate=True,reuse_weights=False)
            
            # regridded no2 trop. columns
            no2_modgrid = regridder(satno2) # , keep_attrs=True
            print('Done with TROPOMI regridding', days, ns)

            #regridder.destroy()
            del regridder
 
            no2_modgrid_all[:,:,ns] = no2_modgrid
            print(' no2 satellite:', np.nanmin(no2_modgrid), np.nanmax(no2_modgrid))

        # daily averaged no2 trop. columns at model grids
        no2_modgrid_avg['nitrogendioxide_tropospheric_column'][nd,:,:] = np.nanmean(np.where(no2_modgrid_all > 0.0, no2_modgrid_all, np.nan), axis=2)

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
              
        # sum up tropopause, needs to be revised to tropopause
        if 'pres_pa_trop' in list(modobj.keys()):
            no2_modgrid_avg[f'{no2varname}trpcol'][nd, :,:] = modobj_tm[f'{no2varname}_col'].where(modobj_tm['pres_pa_mid'] >= modobj_tm['pres_pa_trop']).sum(dim='z').values.squeeze()

        else:
            print('Caution: model tropospheric NO2 column was calculated assuming the model top is the tropopause')
            no2_modgrid_avg[f'{no2varname}trpcol'][nd, :,:] = modobj_tm[f'{no2varname}_col'].sum(dim='z').values.squeeze()
        # --- tropomi ---
        # number of swath
        nswath = len(obsobj[days])

        # array for all swaths
        no2_modgrid_all = np.zeros([ny, nx, nswath], dtype=np.float32)

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

            # regridding from swath grid to model grids
            regridder = xe.Regridder(grid_sat, grid_mod,'bilinear',ignore_degenerate=True,reuse_weights=False)

            # regridded no2 trop. columns
            no2_modgrid = regridder(satno2, keep_attrs=True)
            no2_modgrid_all[:,:,ns] = no2_modgrid[:,:]

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

# Conservative / unstructured TROPOMI NO2 operator (model -> TROPOMI column)
# Reuses the model-agnostic regrid() (regrid_util) + uxarray helpers

def _to_molmol(da):
    """Return model mixing ratio in mol/mol, deciding by VALUE MAGNITUDE.

    The CAM-unstructured reader scales NO2 to ppbV (~0.01-300) but can leave
    the units attribute as 'mol/mol' (stale-attr bug). mol/mol NO2 is
    ~1e-9-1e-7, so a max above ~1e-3 means the values are really ppbV.
    """
    mx = float(np.nanmax(np.abs(da.values)))
    return da * 1e-9 if mx > 1e-3 else da

def interp_vertical_mod2tropomi(obsobj, modobj_swath, variables=("NO2",)):
    """Interpolate model mixing ratio (intensive) onto TROPOMI's layers.
    
    interpolate the model *concentration* in log-pressure onto those layers; the partial
    column is integrated afterward (in :func:`apply_weights_mod2tropomi_no2`)
    using TROPOMI's own layer thickness 
    
    Parameters
    ----------
    obsobj : xr.Dataset
        One TROPOMI granule (time squeezed) with ``pres_pa_mid`` (z, y, x).
    modobj_swath : xr.Dataset
        Model already regridded to the swath pixels, with ``pres_pa_mid``
        (z_mod, y, x) and the requested ``variables``.

    Returns
    -------
    xr.Dataset
        ``variables`` interpolated onto TROPOMI's z layers, dims (z, y, x).
    """
    
    p_mod = np.asarray(modobj_swath["pres_pa_mid"].transpose("z", "y", "x").values)
    p_trop = np.asarray(obsobj["pres_pa_mid"].transpose("z", "y", "x").values)
    nz_t, ny, nx = p_trop.shape
    out = xr.Dataset()
    for var in list(variables):
        src = np.asarray(modobj_swath[var].transpose("z", "y", "x").values)
        dest = np.full((nz_t, ny, nx), np.nan, dtype=float)
        for j in range(ny):
            for i in range(nx):
                mp = p_mod[:, j, i]; mc = src[:, j, i]; tp = p_trop[:, j, i]
                good = np.isfinite(mp) & np.isfinite(mc)
                if good.sum() < 2 or not np.isfinite(tp).any():
                    continue
                order = np.argsort(mp[good])
                
                # Edge behavior: CLAMP (np.interp default), matching TEMPO's
                # _interp_vert. Retrieval levels below the model's lowest
                # mid-level get the lowest-layer mixing ratio; levels above
                # the model top get the top layer (tropopause-masked for
                # tropospheric columns anyway). NaN edges here silently
                # discarded ~85% of pixels: the retrieval surface level
                # (~actual sfc pressure) usually sits below the model's
                # lowest MID-level, so z=0 went NaN and the AK step's
                # isfinite(vmr[0]) mask dropped the entire pixel.
                dest[:, j, i] = np.interp(
                    np.log10(tp), np.log10(mp[good][order]), mc[good][order],
                    #left=np.nan, right=np.nan,
                )
        out[var] = xr.DataArray(dest, dims=("z", "y", "x"))
    return out

# regrid to lat lon 
# In the future, this might need a seperate util file that enables it to be generalizeable to other sat products 
_DEG_PER_M = 1.0 / 111320.0   # degrees latitude per metre (mean Earth)

def _model_lonlat_extent(modobj, pad=0.0):
    """(lonmin, lonmax, latmin, latmax) of the model domain, optionally padded"""
    mlon = np.asarray(modobj["longitude"].values)
    mlat = np.asarray(modobj["latitude"].values)
    return (float(np.nanmin(mlon)) - pad, float(np.nanmax(mlon)) + pad,
            float(np.nanmin(mlat)) - pad, float(np.nanmax(mlat)) + pad)

def _crop_swath_to_extent(o, extent, pad=0.5):
    """Slice a swath granule to the (y,x) index window overlapping the extent box.
    """
    if extent is None:
        return o
    lonmin, lonmax, latmin, latmax = extent
    
    lon_n = "longitude" if "longitude" in o.variables else "lon"  # TROPOMI vs TEMPO
    lat_n = "latitude" if "latitude" in o.variables else "lat"
    lon = np.asarray(o[lon_n].values)
    lat = np.asarray(o[lat_n].values)
    
    inbox = ((lon >= lonmin - pad) & (lon <= lonmax + pad) &
             (lat >= latmin - pad) & (lat <= latmax + pad))
    if not inbox.any():
        return None
    hdims = o[lon_n].dims                      # (y, x) or (x, y)
    a0 = np.where(inbox.any(axis=1))[0]
    a1 = np.where(inbox.any(axis=0))[0]
    sl = {hdims[0]: slice(int(a0.min()), int(a0.max()) + 1),
          hdims[1]: slice(int(a1.min()), int(a1.max()) + 1)}
    return o.isel(sl)

def _swath2latlon(swath, data_vars, res, extent, units = "deg", method="radius_mean"):
    """Regrid swath-paired fields (y, x) onto a regular lat/lon grid 

    use radius mean averaging 
    
    returns a Dataset on dims (lat, lon).
    """
    from melodies_monet.util.regrid_util import regrid

    lonmin, lonmax, latmin, latmax = extent
    if units in ("km", "m"):
        res_m = float(res) * 1000.0 if units == "km" else float(res)
        dlat = res_m * _DEG_PER_M
        midlat = 0.5 * (latmin + latmax)
        dlon = res_m * _DEG_PER_M / max(np.cos(np.deg2rad(midlat)), 0.1)
    else:  # degrees
        dlat = dlon = float(res)

    tlon1 = np.arange(lonmin, lonmax + dlon, dlon)
    tlat1 = np.arange(latmin, latmax + dlat, dlat)
    tlon2, tlat2 = np.meshgrid(tlon1, tlat1)        # (lat, lon)

    # handle the lat lon naming between satellites 
    # tropomi carry long / lat on (y,x) where tempo does lon / lat on (x,y)
    # rename to longintude latitude 
    ren = {}
    if "longitude" not in swath.variables and "lon" in swath.variables:
        ren["lon"] = "longitude"
    if "latitude" not in swath.variables and "lat" in swath.variables:
        ren["lat"] = "latitude"
    if ren:
        swath = swath.rename(ren)
    _to_coord = [c for c in ("longitude", "latitude") if c in swath.data_vars]
    if _to_coord:
        swath = swath.set_coords(_to_coord)

    # insert a conservative regriding option 
    # Builds the swath source mesh from corner bounds and a rectilinear target
    # mesh, then mesh-to-mesh conservative regrid
    # Fall back to radius_mean if bounds are missing
    
    if method in _CONSERVATIVE and "longitude_bounds" in swath.variables:
        try:
            import uxarray as ux

            swath_grid, _ = _tropomi_swath_mesh(swath)
            _hd = tuple(swath["longitude"].dims)
            _flat = (
                swath[data_vars].stack(n_face=_hd).reset_index("n_face", drop=True)
            )
            _flat = _flat.drop_vars(
                [c for c in _flat.coords if "n_face" in _flat[c].dims], errors="ignore")
            _src = ux.UxDataset(_flat, uxgrid=swath_grid)
            _tgt = ux.Grid.from_structured(lon=tlon1, lat=tlat1)
            _out = regrid(_src, method=method, target_grid=_tgt)
            _out = _out.where(_out != 0)        # empty cells go to  NaN, not 0

            _nlat, _nlon = tlat1.size, tlon1.size
            _ncell = _nlat * _nlon
            _fd = next((d for d in _out.dims if _out.sizes.get(d) == _ncell), None)
            if _fd is None:
                raise ValueError("conservative target face dim not found after regrid")
            _res = xr.Dataset()
            for v in _out.data_vars:
                da = _out[v]
                if _fd not in da.dims:
                    continue
                t = da.transpose(..., _fd)
                arr = np.asarray(t.values).reshape(t.shape[:-1] + (_nlat, _nlon))
                _res[v] = xr.DataArray(arr, dims=t.dims[:-1] + ("y", "x"),
                                       attrs=dict(da.attrs))
            return _res.assign_coords(
                longitude=(("y", "x"), tlon2), latitude=(("y", "x"), tlat2),
                x=("x", tlon1), y=("y", tlat1))
        except Exception as e:
            print(f"_swath2latlon: conservative regrid failed ({e!r}); "
                  "falling back to radius_mean.", flush=True)

    hdims = list(swath["longitude"].dims)            # (y, x) or (x, y)

    # instead of radius mean use a box mean gridding (faster)
    
    # flat = (
    #     swath[data_vars]
    #     .stack(pixel=hdims)
    #     .reset_index("pixel", drop=True)
    # )

    plon = np.asarray(swath["longitude"].values).ravel()
    plat = np.asarray(swath["latitude"].values).ravel()
    nlat, nlon = tlat1.size, tlon1.size
    ncell = nlat * nlon
    
    # search radius ~ cell size, but at least ~5 km (sensor footprint) so a
    # finer-than-sensor grid fills from the nearest pixel instead of going empty.
    # radius = max(float(dlat), 0.05)
    # out = regrid(flat, target={"lon": tlon2, "lat": tlat2},
    #              method="radius_mean", radius=radius, target_dims=("lat", "lon"))
    # out = out.where(out != 0)        # empty cells -> NaN, not 0

    # out = out.assign_coords(lon=("lon", tlon1), lat=("lat", tlat1))
    ix = np.floor((plon - lonmin) / dlon).astype(np.intp)
    iy = np.floor((plat - latmin) / dlat).astype(np.intp)
    ingrid = ((ix >= 0) & (ix < nlon) & (iy >= 0) & (iy < nlat)
              & np.isfinite(plon) & np.isfinite(plat))
    cell = iy * nlon + ix                            # flat (lat-major) cell id
    
    out = xr.Dataset()
    for v in data_vars:
        da = swath[v]
        extra = [d for d in da.dims if d not in hdims]
        arr = np.asarray(da.transpose(*extra, *hdims).values, dtype=float).reshape(-1, plon.size)
        res = np.full((arr.shape[0], ncell), np.nan)
        for k in range(arr.shape[0]):
            good = ingrid & np.isfinite(arr[k])
            c = cell[good]
            ssum = np.bincount(c, weights=arr[k][good], minlength=ncell)
            scnt = np.bincount(c, minlength=ncell)
            with np.errstate(invalid="ignore", divide="ignore"):
                m = ssum / scnt
            m[scnt == 0] = np.nan
            res[k] = m
        out[v] = xr.DataArray(
            res.reshape(tuple(da.sizes[d] for d in extra) + (nlat, nlon)),
            dims=tuple(extra) + ("y", "x"), attrs=dict(da.attrs))

    # 2-D lon/lat + 1-D x/y so the result flows through structured-sat plotting
    return out.assign_coords(
        longitude=(("y", "x"), tlon2), latitude=(("y", "x"), tlat2),
        x=("x", tlon1), y=("y", tlat1))

    # # want to make sure these regrided lat lon pairs can just run through the existing plotting 
    # return out.rename({"lat": "y", "lon": "x"})

def _wmean(da, w):
    """Weighted spatial mean of da with weights w, ignoring NaN

    Returns a 0-d (scalar) DataArray; NaN if nothing valid.
    """
    w = w.where(np.isfinite(da) & np.isfinite(w) & (w > 0))
    den = w.sum(skipna=True)
    return xr.where(den > 0, (da * w).sum(skipna=True) / den, np.nan)

def _attach_obs_err(swath, o, obs_var):
    """Attach the per-pixel retrieval error to swath as _obs_err
    
    Used by the 'series' target to inverse-variance weight the observations
    
    Looks for '<obs_var>_precision' then '<obs_var>_uncertainty' in the
    granule. No-op if neither was read in (then 'series' falls back to area
    weighting for obs).
    
    """
    for _name in (obs_var + "_precision", obs_var + "_uncertainty"):
        if _name in o.variables:
            err = o[_name]
            if "time" in err.dims:
                err = err.squeeze("time", drop=True)
            swath["_obs_err"] = err
            return

def _swath2series(swath, data_vars, obs_var=None,
                  obs_weight="inverse_variance", model_weight="area"):
    
    """Collapse a swath to ONE weighted domain value per variable (a time point).

    The 'series' target saves a *time vector* instead of a map: every
    granule/overpass is reduced to a single representative number, and the
    orchestrator stamps it with the overpass time and concatenates to (time,).
    Obs and model are reduced over the *same* sampled pixels, so they stay
    directly comparable (the absolute level still reflects which footprint each
    overpass sampled).

    Default weighting:
      * model field        -> AREA-weighted by cos(latitude), so the domain
        mean is area-representative rather than biased toward where swath
        pixels happen to be dense.
      * obs field (obs_var) -> INVERSE-VARIANCE weighted (weight = 1/sigma^2)
        using the per-pixel retrieval error in swath['_obs_err'] (see
        _attach_obs_err), so noisy retrievals count for less. Falls back to
        area weighting if no error field is attached.

    obs_weight / model_weight may be overridden: 'inverse_variance' | 'area'
    | 'equal'. Returns a Dataset of scalars (one per variable in data_vars).
    """
    
    lat = swath["latitude"] if "latitude" in swath.variables else swath["lat"]
    w_area = np.cos(np.deg2rad(lat))

    out = xr.Dataset()
    for v in data_vars:
        da = swath[v]
        mode = obs_weight if v == obs_var else model_weight
        if v == obs_var and mode == "inverse_variance" and "_obs_err" in swath.variables:
            out[v] = _wmean(da, 1.0 / (swath["_obs_err"] ** 2))   # inverse-variance
            used = "inverse_variance"
        elif mode == "equal":
            out[v] = da.mean(skipna=True)                          # unweighted
            used = "equal"
        else:
            out[v] = _wmean(da, w_area)                            # area (cos lat)
            # requested inverse_variance but no _obs_err read in to area
            used = ("area (inverse_variance fallback: no _obs_err)"
                    if v == obs_var and mode == "inverse_variance" else "area")
        # record how this scalar was built so saved series files are auditable
        out[v].attrs = {**dict(da.attrs), "series_weighting": used}
        
    return out

def _swath_to_target(swath, modobj, method, data_vars, target, res, extent, units="deg" ):
    """Regrid the paired swath onto the requested target space

    if target = model, use model's native unstructured grid via tropomi_swath2mod

    if target = obs, use a regular lat lon grid via _swath2latlon 

    if target = series, single weighted domain value per overpass (time vector)

    """
    if target == "model":
        return _tropomi_swath2mod(swath, modobj, method, data_vars)
    if target == "obs":
        return _swath2latlon(swath, data_vars, res, extent, units=units, method=method)
    if target == "series":
        _ov = data_vars[1] if len(data_vars) > 1 else None
        return _swath2series(swath, data_vars, obs_var=_ov)
    raise ValueError(f"regrid_target {target!r} not understood; use 'model', 'obs', or 'series'.")

def apply_weights_mod2tropomi_no2(obsobj, modobj_on_tropomi_layers, species="NO2"):
    """Apply the TROPOMI averaging kernel to a model NO2 profile.

    Mirrors :func:`apply_weights_mod2tempo_no2` but uses TROPOMI's
    averaging kernel instead of scattering weights:

      AK_trop = averaging_kernel * (amf_total / amf_troposphere)
      subcol  = vmr * (dp/g) * (NA/M_air) / 1e4         # molec/cm2 per layer
      VCD     = sum over tropospheric layers (p >= tropopause) of AK_trop*subcol

    Parameters
    ----------
    obsobj : xr.Dataset
        One TROPOMI granule (time squeezed): averaging_kernel (z,y,x),
        air_mass_factor_troposphere/_total (y,x), pres_pa_mid (z,y,x),
        pres_pa_int (z_stagg,y,x), tm5_tropopause_pressure (y,x).
    modobj_on_tropomi_layers : xr.Dataset
        Model ``species`` mixing ratio on TROPOMI's z layers (z,y,x), from
        :func:`interp_vertical_mod2tropomi`. Units auto-detected (ppbV/mol/mol).

    Returns
    -------
    xr.DataArray
        Model NO2 tropospheric column with the AK applied, molec/cm2, (y, x).
    """
    g, M_air, NA = 9.80665, 0.0289644, 6.022e23

    vmr = _to_molmol(modobj_on_tropomi_layers[species]).transpose("z", "y", "x")
    pint = obsobj["pres_pa_int"].transpose("z_stagg", "y", "x")
    dp = np.abs(
        pint.isel(z_stagg=slice(0, -1)).values
        - pint.isel(z_stagg=slice(1, None)).values
    )
    dp = xr.DataArray(dp, dims=("z", "y", "x"))
    subcol = vmr * dp * (NA / (g * M_air) / 1e4)            # molec/cm2 per layer

    ak = obsobj["averaging_kernel"].transpose("z", "y", "x")
    ak_trop = ak * (obsobj["air_mass_factor_total"]
                    / obsobj["air_mass_factor_troposphere"])
    trop = (obsobj["pres_pa_mid"].transpose("z", "y", "x")
            >= obsobj["tm5_tropopause_pressure"])

    vcd = (ak_trop * subcol).where(trop).sum("z", skipna=True)
    vcd = vcd.where(np.isfinite(vmr.isel(z=0)))

    # AK sanity diagnostic (enable with logging DEBUG): AK/raw equals
    # AMF_model/AMF_retrieval. This one-liner exposed the xregrid doubling.
    if logger.isEnabledFor(logging.DEBUG):
        try:
            raw = subcol.where(trop).sum("z", skipna=True)
            ratio = (vcd / raw).where(raw != 0)
            logger.debug(
                "[AK] TROPOMI %s: raw_col=%.2e AK-applied=%.2e AK/raw=%.2f "
                "(=AMF_mod/AMF_ret)", species,
                float(np.nanmean(raw.values)), float(np.nanmean(vcd.values)),
                float(np.nanmean(ratio.values)))
        except Exception as e:  # noqa: BLE001
            logger.debug("[AK] TROPOMI %s ratio diag skipped: %r", species, e)
        
    vcd.attrs = {
        "units": "molecules/cm2",
        "description": "model NO2 tropospheric column after applying TROPOMI averaging kernel",
        "history": "Created by MELODIES-MONET, apply_weights_mod2tropomi_no2",
    }
    return vcd.where(np.isfinite(vcd))

# Orchestrator: conservative/unstructured TROPOMI NO2 pairing.
# Mirrors the TEMPO regrid_and_apply_weights flow but with the TROPOMI
# averaging-kernel operator. Reuses only the shared, product-agnostic regrid
# primitives (regrid_util.regrid, uxarray_util.open_uxgrid /
# uxgrid_from_corner_bounds) -- no dependency on the TEMPO utility.

_TROPOMI_NO2_VAR = "nitrogendioxide_tropospheric_column"
_MOL_M2_TO_MOLEC_CM2 = 6.02214e19
_CONSERVATIVE = ("conservative", "conservative_normed")

def _tropomi_swath_mesh(o):
    """Build a uxarray Grid from a TROPOMI granule's pixel corner bounds.

    Returns (grid, (ny, nx)). Flatten order is row-major over (y, x), matching
    how the n_face result is reshaped back.
    """
    from melodies_monet.util.uxarray_util import uxgrid_from_corner_bounds

    olon = np.asarray(o["longitude"].values)  # (y, x)
    ny, nx = olon.shape
    clon = np.asarray(o["longitude_bounds"].values).reshape(ny * nx, -1)
    clat = np.asarray(o["latitude_bounds"].values).reshape(ny * nx, -1)
    return uxgrid_from_corner_bounds(clon, clat), (ny, nx)


def _mod2tropomi_swath(modobj, o, method, mod_vars, grid_file):
    """Regrid model fields onto the TROPOMI swath pixels (y, x).

    conservative -> mesh-to-mesh via xregrid (model SCRIP/MPAS mesh -> swath
    cells from bounds); nearest/radius -> cKDTree at pixel centers. Returns a
    Dataset with each requested var on (..., y, x).
    """
    from melodies_monet.util.regrid_util import regrid

    msrc = modobj[[v for v in mod_vars if v in modobj.variables]]
    olon = np.asarray(o["longitude"].values)
    olat = np.asarray(o["latitude"].values)
    ny, nx = olon.shape

    if method in _CONSERVATIVE:
        from melodies_monet.util.uxarray_util import (
            faces_to_grid, subset_model_source)

        swath_grid, _ = _tropomi_swath_mesh(o)
        # Subset the model mesh to the swath bbox before building conservative
        # weights (far faces contribute zero). Mirrors the TEMPO forward path.
        src, src_grid = subset_model_source(
            msrc, grid_file,
            float(np.nanmin(olon)), float(np.nanmax(olon)),
            float(np.nanmin(olat)), float(np.nanmax(olat)),
            label="_mod2tropomi_swath",
        )
        out = regrid(src, method=method, src_grid=src_grid,
                     target_grid=swath_grid)

        res = faces_to_grid(out, (ny, nx), ("y", "x"))
        return res.assign_coords(longitude=(("y", "x"), olon),
                                 latitude=(("y", "x"), olat))

    # nearest family
    out = regrid(msrc, target={"lon": olon, "lat": olat},
                 method="nearest_s2d", target_dims=("y", "x"))
    return out


def _tropomi_swath2mod(swath, modobj, method, data_vars):
    """Regrid swath-paired fields (y, x) onto the unstructured model columns.

    conservative -> mesh-to-mesh via xregrid (swath cells -> model mesh);
    else -> cKDTree within-radius mean. Returns a Dataset on the model
    column dim with longitude/latitude attached.
    """
    import uxarray as ux
    from melodies_monet.util.regrid_util import regrid
    from melodies_monet.util.uxarray_util import open_uxgrid

    col_dim = modobj["longitude"].dims[0]
    mlon = np.asarray(modobj["longitude"].values).ravel()
    mlat = np.asarray(modobj["latitude"].values).ravel()
    grid_file = (modobj.attrs.get("mio_scrip_file")
                 or modobj.attrs.get("mio_mesh_file"))

    if method in _CONSERVATIVE and "longitude_bounds" in swath.variables:
        from melodies_monet.util.uxarray_util import (
            flatten_to_faces, subset_mesh_to_bbox)
        swath_grid, (ny, nx) = _tropomi_swath_mesh(swath)
        flat = flatten_to_faces(swath[data_vars], ("y", "x"))
        src = ux.UxDataset(flat, uxgrid=swath_grid)
        model_grid = open_uxgrid(grid_file)
        n_col = int(model_grid.n_face)

        # Subset the TARGET mesh to the swath bbox before building weights, then
        # scatter results back to the full mesh (faces outside the swath go NaN).
        # Mirrors the TEMPO backward path.
        out = None
        try:
            _slon = np.asarray(swath["longitude"].values)
            _slat = np.asarray(swath["latitude"].values)
            _keep, _subgrid = subset_mesh_to_bbox(
                grid_file,
                float(np.nanmin(_slon)), float(np.nanmax(_slon)),
                float(np.nanmin(_slat)), float(np.nanmax(_slat)))
            if _subgrid is not None:
                out_sub = regrid(src, method=method, target_grid=_subgrid)
                _sfd = next((d for d in out_sub.dims
                             if out_sub.sizes[d] == _keep.size), None)
                out = xr.Dataset(attrs=dict(out_sub.attrs))
                for v in out_sub.data_vars:
                    da = out_sub[v]
                    if _sfd is None or _sfd not in da.dims:
                        out[v] = da
                        continue
                    da = da.transpose(..., _sfd)
                    arr = np.asarray(da.values)
                    full = np.full(arr.shape[:-1] + (n_col,), np.nan, dtype="float64")
                    full[..., _keep] = arr
                    out[v] = xr.DataArray(full, dims=da.dims[:-1] + (_sfd,),
                                          attrs=dict(da.attrs))
                print(f"_tropomi_swath2mod: target mesh subset {n_col} -> "
                      f"{_keep.size} faces", flush=True)
        except Exception as e:  # noqa: BLE001
            print(f"_tropomi_swath2mod: target subset skipped ({e!r}); full mesh.",
                  flush=True)
            out = None
        if out is None:
            out = regrid(src, method=method, target_grid=model_grid)
            
        # Conservative regrid fills model cells with no swath overlap with exactly 0.
        out = out.where(out != 0)
        
        # n_col = int(model_grid.n_face)
        d = next((dd for dd in out.dims if out.sizes[dd] == n_col), None)
        if d is not None and d != col_dim:
            out = out.rename({d: col_dim})
        return out.assign_coords({"longitude": (col_dim, mlon),
                                  "latitude": (col_dim, mlat)})

    # nearest / radius_mean fallback
    flat = (
        swath[data_vars]
        .stack(pixel=("y", "x"))
        .reset_index("pixel", drop=True)
        .set_coords(["longitude", "latitude"])
    )
    out = regrid(flat, target={"lon": mlon, "lat": mlat},
                 method="radius_mean", radius=0.1, target_dims=(col_dim,))

    # make sure to fill 0s with nans
    return out.where(out != 0)

# 07062026 simplify this code and combine where redundant 

# time should be non specific to the satellite 
def _granule_time(o):
    """Granule overpass time for model matching.

    Prefer 'time_granule' (real per-measurement times, take the mean); fall
    back to the 'time' coord (start-of-day reference) if absent. None if
    neither is usable.
    """
    if "time_granule" in o.variables:
        tg = np.asarray(o["time_granule"].values).ravel().astype("datetime64[ns]")
        tg = tg[~np.isnat(tg)]
        return (np.array(tg.astype("int64").mean(), dtype="int64")
                .astype("datetime64[ns]") if tg.size else None)
    return o["time"].values if "time" in o.coords else None
    
# model time at the granule sat overpass time 
def _model_at_time(modobj, gtime):
    """Model state at a granule's overpass time.

    Linear time interpolation, clamped to the model's time span (an overpass
    just outside the model window snaps to the first/last model time);
    nearest-neighbor as a fallback; passthrough when the model has no time
    dim (or first step when the overpass time is unknown).
    """
    if "time" not in modobj.dims:
        return modobj
    if gtime is None:
        return modobj.isel(time=0)
    tsel = gtime if np.ndim(gtime) == 0 else np.asarray(gtime).ravel()[0]
    mtimes = modobj["time"].values
    tsel = min(max(tsel, mtimes.min()), mtimes.max())
    try:
        return modobj.interp(time=tsel)
    except Exception:
        return modobj.sel(time=tsel, method="nearest")


###################### Regrid and application of weights ###################### 
########################### instrument specific ############################### 

def regrid_and_apply_weights_tropomi(obsobj, modobj, species=["NO2"],
                                     method="conservative", qa_min=0.75,
                                     regrid_target="model", obs_grid_res=0.1, obs_grid_units="deg", obs_grid_extent=None, 
                                     save_vars=None):
    
    """Pair an unstructured model with TROPOMI L2 NO2 (AK applied).

    wrapper for _pair_tropomi. for each granule, forward regrid model to the swath, interpolate model profile 
    onto tropomi layers, apple the AK, then regrid the AK'd model column AND obs column to requested 
    regrid targets. Concatenated along time 

    Parameters
    ----------
    obsobj : dict[str, list[xr.Dataset]]
        Output of the generic TROPOMI reader (tropomi_l2.open_datasets):
        keyed by date, each value a list of orbit granules.
    modobj : xr.Dataset
        Unstructured model with longitude/latitude on its column dim,
        NO2 (ppbV or mol/mol), pres_pa_mid, and mio_scrip_file/mio_mesh_file.
    species : list[str]
        Model species name(s); species[0] is paired.
    method : str
        Regrid method (conservative recommended; nearest_s2d/radius_mean ok).

    Returns
    -------
    dict[str, xr.Dataset]
        Per regrid target: paired model + obs NO2 tropospheric columns
        (molec/cm2), stacked along ``time`` (one entry per granule).
    """
    
    return _pair_tropomi("no2", obsobj, modobj, species=species, method=method,
                         qa_min=qa_min, regrid_target=regrid_target,
                         obs_grid_res=obs_grid_res,
                         obs_grid_units=obs_grid_units,
                         obs_grid_extent=obs_grid_extent, save_vars=save_vars)
    
# tropomi HCHO 
# single tropo AMF and total columng averaging kernel 
# use generic tropomi_l2 reader + shared regrid / interp 

def apply_weights_mod2tropomi_hcho(obsobj, modobj_on_tropomi_layers, species="CH2O"):
    """Apply the TROPOMI HCHO averaging kernel to a model formaldehyde profile.

    Returns the model HCHO column with the AK applied (molec/cm2), dims (y, x).
    """
    g, M_air, NA = 9.80665, 0.0289644, 6.022e23

    vmr = _to_molmol(modobj_on_tropomi_layers[species]).transpose("z", "y", "x")

    pint = obsobj["pres_pa_int"].transpose("z_stagg", "y", "x")
    dp = np.abs(
        pint.isel(z_stagg=slice(0, -1)).values
        - pint.isel(z_stagg=slice(1, None)).values
    )
    dp = xr.DataArray(dp, dims=("z", "y", "x"))
    subcol = vmr * dp * (NA / (g * M_air) / 1e4)            # molec/cm2 per layer

    # AK * AMF 
    ak = obsobj["averaging_kernel"].transpose("z", "y", "x")
    ak = ak * obsobj["formaldehyde_tropospheric_air_mass_factor"] # hopefully this isnt hardcoded otherwise will need to pull this in from YAML

    vcd = (ak * subcol).sum("z", skipna=True)
    vcd = vcd.where(np.isfinite(vmr.isel(z=0)))
    
    vcd.attrs = {
        "units": "molecules/cm2",
        "description": "model HCHO column after applying TROPOMI averaging kernel",
        "history": "Created by MELODIES-MONET, apply_weights_mod2tropomi_hcho",
    }
    return vcd.where(np.isfinite(vcd))

def regrid_and_apply_weights_tropomi_hcho(obsobj, modobj, species=["CH2O"],
                                          method="conservative", qa_min=0.5,
                                          regrid_target="model", obs_grid_res=0.1, obs_grid_units="deg", 
                                          obs_grid_extent=None, save_vars = None):
    """
    Pair an unstructured model with TROPOMI L2 HCHO (AK applied).

    defaults to regridding to model space via "regrid_target"

    regrid_target="model" ; model space 
    regrid_target="obs" ; obs space
    regrid_target="series" ; one weighted domain value per overpass
    
    """
    return _pair_tropomi("hcho", obsobj, modobj, species=species, method=method,
                         qa_min=qa_min, regrid_target=regrid_target,
                         obs_grid_res=obs_grid_res,
                         obs_grid_units=obs_grid_units,
                         obs_grid_extent=obs_grid_extent, save_vars=save_vars)

# co 
_TROPOMI_CO_VAR = "carbonmonoxide_total_column"

def apply_weights_mod2tropomi_co(obsobj, modobj_on_tropomi_layers, species="CO"):
    """Apply the TROPOMI CO column averaging kernel to a model CO profile.

    """
    g, M_air, NA = 9.80665, 0.0289644, 6.022e23

    vmr = _to_molmol(modobj_on_tropomi_layers[species]).transpose("z", "y", "x")

    pint = obsobj["pres_pa_int"].transpose("z_stagg", "y", "x")
    dp = np.abs(
        pint.isel(z_stagg=slice(0, -1)).values
        - pint.isel(z_stagg=slice(1, None)).values
    )
    dp = xr.DataArray(dp, dims=("z", "y", "x"))
    subcol = vmr * dp * (NA / (g * M_air) / 1e4)            # molec/cm2 per layer

    ak = obsobj["column_averaging_kernel"].transpose("z", "y", "x")  # dimensionless 
    vcd = (ak * subcol).sum("z", skipna=True)
    vcd = vcd.where(np.isfinite(vmr.isel(z=0)))
    
    vcd.attrs = {
        "units": "molecules/cm2",
        "description": "model CO column after applying TROPOMI column averaging kernel",
        "history": "Created by MELODIES-MONET, apply_weights_mod2tropomi_co",
    }
    return vcd.where(np.isfinite(vcd))

def regrid_and_apply_weights_tropomi_co(obsobj, modobj, species=["CO"],
                                        method="conservative", qa_min=0.5,
                                        regrid_target="model", obs_grid_res=0.1, obs_grid_units="deg", obs_grid_extent=None, 
                                        save_vars=None):
    """Pair an unstructured model with TROPOMI L2 CO (column AK applied).

    Thin wrapper over _pair_tropomi('co', ...). Prefers the destriped
    'carbonmonoxide_total_column_corrected' when it was read in.
    """
    return _pair_tropomi("co", obsobj, modobj, species=species, method=method,
                         qa_min=qa_min, regrid_target=regrid_target,
                         obs_grid_res=obs_grid_res,
                         obs_grid_units=obs_grid_units,
                         obs_grid_extent=obs_grid_extent, save_vars=save_vars)



# generalize the ability to easily add more tropomi species 
# to add new satelite variables ncdump -h "$F" | grep -iE "group: (PRODUCT|DETAILED_RESULTS|INPUT_DATA)|carbonmonoxide_total_column|column_averaging_kernel|pressure_levels|qa_value|layer =|:units|:multiplication_factor"

_TROPOMI_PRODUCTS = {
    "no2": dict(obs_var="nitrogendioxide_tropospheric_column",
                ak_apply=apply_weights_mod2tropomi_no2,
                species="NO2", qa_min=0.75, prefer_corrected=False),
    "hcho": dict(obs_var="formaldehyde_tropospheric_vertical_column",
                 ak_apply=apply_weights_mod2tropomi_hcho,
                 species="CH2O", qa_min=0.5, prefer_corrected=False),
    "co": dict(obs_var="carbonmonoxide_total_column",
               ak_apply=apply_weights_mod2tropomi_co,
               species="CO", qa_min=0.5, prefer_corrected=True),
}

def _tropomi_swath_pixels(swath, sp, obs_var, gtime, save_vars=None):
    """Flatten one paired TROPOMI granule to a native-pixel vector.

    The 'swath' target keeps every footprint (no gridding): a 1-D ``obs``
    vector carrying longitude/latitude/time as per-pixel coords, plus the
    per-pixel retrieval uncertainty (renamed ``<obs_var>_uncertainty``).
    All-NaN pixels dropped.
    """
    sdims = tuple(swath[obs_var].dims)                       # (y, x)

    # can save out QA values 
    _extra = [v for v in (save_vars or []) if v in swath.variables]
    keep_vars = ([sp, obs_var]
                 + (["_obs_err"] if "_obs_err" in swath.variables else [])
                 + _extra)
    
    flat = swath[keep_vars].stack(obs=sdims).reset_index("obs", drop=True)
    lon = np.asarray(
        swath["longitude"].stack(obs=sdims).reset_index("obs", drop=True).values)
    lat = np.asarray(
        swath["latitude"].stack(obs=sdims).reset_index("obs", drop=True).values)
    n = flat.sizes["obs"]
    _t = (np.datetime64(gtime if np.ndim(gtime) == 0
                        else np.asarray(gtime).ravel()[0])
          if gtime is not None else np.datetime64("NaT"))
    flat = flat.assign_coords(longitude=("obs", lon), latitude=("obs", lat),
                              time=("obs", np.full(n, _t)))
    if "_obs_err" in flat.data_vars:
        flat = flat.rename({"_obs_err": obs_var + "_uncertainty"})
        
    # carry pixel CORNER bounds (obs, corner) for footprint-polygon rendering
    for _bn in ("longitude_bounds", "latitude_bounds"):
        if _bn in swath.variables:
            _bb = swath[_bn].stack(obs=sdims).reset_index("obs", drop=True)
            _cd = next((d for d in _bb.dims if d != "obs"), None)
            if _cd is not None:
                flat[_bn] = _bb.transpose("obs", _cd)
                
    keep = np.zeros(n, dtype=bool)
    for v in (sp, obs_var):
        if v in flat.data_vars:
            keep |= np.isfinite(np.asarray(flat[v].values))
    return flat.isel(obs=np.where(keep)[0])


def _pair_tropomi(product, obsobj, modobj, species=None, method="conservative",
                  qa_min=None, regrid_target="model", obs_grid_res=0.1,
                  obs_grid_units="deg", obs_grid_extent=None, save_vars = None):
    
    """Shared TROPOMI pairing loop; the per-product science comes from
    _TROPOMI_PRODUCTS[product].

    For each granule: crop to the domain extent (if specified), take the model at the
    overpass time, forward-regrid model to the swath, interpolate to the
    retrieval layers, apply the averaging kernel, QA-filter obs and model
    together, then regrid the paired swath to every requested target
    ('model' / 'obs' / 'series'), stamp it with the overpass time and
    concatenate along time.
    """
    
    cfg = _TROPOMI_PRODUCTS[product]
    obs_var = cfg["obs_var"]
    sp = (species or [cfg["species"]])[0]
    if qa_min is None:
        qa_min = cfg["qa_min"]
    save_vars = list(save_vars or [])

    grid_file = (modobj.attrs.get("mio_scrip_file")
                 or modobj.attrs.get("mio_mesh_file"))

    targets = [regrid_target] if isinstance(regrid_target, str) else list(regrid_target)
    extent = (tuple(obs_grid_extent) if obs_grid_extent
              else _model_lonlat_extent(modobj))

    # Flatten dict[date -> list[granule]] to a flat granule list.
    granules = []
    for v in obsobj.values():
        granules.extend(v if isinstance(v, list) else [v])

    # warn once when the product has a destriped column but only the raw
    # (stripey) one was read in
    # this particularly seems to affect the CO product. 
    if cfg["prefer_corrected"]:
        if not any(obs_var + "_corrected" in o.variables for o in granules):
            print(f"TROPOMI {product}: pairing the raw (stripey) column. The "
                  f"product ships a destriped variable -- add "
                  f"'{obs_var}_corrected: {{}}' to the obs variables in the "
                  "control YAML.", flush=True)

    out_by = {t: [] for t in targets}

    # warn about heavy 3d vars saving out 
    _warned_3d = False
    
    for o in granules:
        if "time" in o.dims:
            o = o.squeeze("time", drop=False)

        # crop the orbit to the region of interest before the expensive pairing
        o = _crop_swath_to_extent(o, extent)
        if o is None:
            continue

        gtime = _granule_time(o)
        mod_t = _model_at_time(modobj, gtime)

        # model to swath to retrieval layers to AK-applied column (molec/cm2)
        on_swath = _mod2tropomi_swath(mod_t, o, method, [sp, "pres_pa_mid"],
                                      grid_file)
        prof_t = interp_vertical_mod2tropomi(o, on_swath, [sp])

        model_col = cfg["ak_apply"](o, prof_t, sp)

        src = (obs_var + "_corrected"
               if cfg["prefer_corrected"] and obs_var + "_corrected" in o.variables
               else obs_var)
        obs_col = o[src] * _MOL_M2_TO_MOLEC_CM2               # molec/cm2
        obs_col.attrs["units"] = "molecules/cm2"   # granule attr says mol m-2

        # QA filter: obs and model masked together so pixel sets stay aligned
        if qa_min and "qa_value" in o.variables:
            qa = o["qa_value"]
            if float(np.nanmax(np.asarray(qa.values))) > 1.5:
                qa = qa / 100.0
            good = qa >= qa_min
            obs_col = obs_col.where(good)
            model_col = model_col.where(good)

        swath = xr.Dataset({sp: model_col, obs_var: obs_col,
             "latitude_bounds": o["latitude_bounds"],
             "longitude_bounds": o["longitude_bounds"]},
            coords={"longitude": o["longitude"], "latitude": o["latitude"]},
        )
        swath[obs_var].attrs["source_variable"] = src

        _save_present = [v for v in save_vars
                         if v in o.variables and v not in swath.variables]
        for v in _save_present:
            swath[v] = o[v]
            
        # per-pixel retrieval error for inverse-variance weighting in 'series'
        _attach_obs_err(swath, o, obs_var)
        
        for t in targets:
            if t == "swath":
                # native pixel vector (no gridding); time is a per-pixel coord,
                # so it concatenates along 'obs', not 'time'.
                out_by[t].append(_tropomi_swath_pixels(swath, sp, obs_var, gtime,
                                                        save_vars=_save_present))
                continue
                
            # For gridded targets only carry 2-D per-pixel save vars (AMF, sza,
            # cloud, qa). keep the full
            # AK on the 'swath' product instead.
            _grid_save = [v for v in _save_present if swath[v].ndim <= 2]
            if not _warned_3d and len(_grid_save) < len(_save_present):
                _dropped = [v for v in _save_present if swath[v].ndim > 2]
                print(f"TROPOMI {product}: not regridding layer-resolved save "
                      f"var(s) {_dropped} to '{t}' (preserved on 'swath' only).",
                      flush=True)
                _warned_3d = True
                
            on = _swath_to_target(swath, modobj, method, [sp, obs_var] + _grid_save,
                                  t, obs_grid_res, extent, units=obs_grid_units)
            
            if gtime is not None:
                tval = gtime if np.ndim(gtime) == 0 else np.asarray(gtime).ravel()[0]
                on = on.expand_dims(time=[np.datetime64(tval)])
            out_by[t].append(on)

    def _cat(t, lst):
        if not lst:
            return xr.Dataset()
        return xr.concat(lst, dim="obs" if t == "swath" else "time")

    return {t: _cat(t, lst) for t, lst in out_by.items()}
    


    