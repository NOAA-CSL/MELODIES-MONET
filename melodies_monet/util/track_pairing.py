
"""
track_pairing.py
================

Pairing utilities for **along-track model output** — CAM "1s_1pt" history
streams that sample the model at an aircraft's lat/lon/time directly

Needs to be paired via time-matching + vertical interpolation to the obs
pressure

this file is specific to cesm 

"""

import json

import numpy as np
import pandas as pd
import xarray as xr

# dry gas const & gravity
RD = 287.05
G0 = 9.80665

def _track_times(m):
    """Model 'time' as a pandas DatetimeIndex.

    CESM history uses a 'noleap'/'365_day' calendar, so xarray decodes time as
    cftime objects that ``pd.to_datetime`` can't parse
    """
    vals = m["time"].values
    if np.issubdtype(vals.dtype, np.datetime64):
        return pd.DatetimeIndex(vals)
    return pd.DatetimeIndex([pd.Timestamp(str(t)) for t in vals])
    
def pair_track_model(model_file, obs_file, mapping, resample="60s",
                     obs_scale=None, fill_below=-9999.0, gap_warn_hours=1.0):
    """Pair a CAM '1s_1pt' along-track model file with aircraft obs.

    Parameters
    ----------
    model_file : str
        Path to the along-track model file (
            dims ``ncol``, ``lev``
            per-sample ``time``/``lat``/``lon`` coords and ``PMID`` in Pa).
    obs_file : str
        Path to the aircraft obs (netCDF). Needs time, pressure_obs (Pa), altitude, lat/lon, and
        the mapped obs variables.
    mapping : dict
        model_var: obs_var
    resample : str
        Pandas resample rule applied to the obs before pairing
    obs_scale : dict, optional
        Per-obs-var unit scaling applied after fill-masking
    fill_below : float
        Obs values ``<= fill_below`` are set to NaN Default -9999.
    gap_warn_hours : float
        Warn if obs file doesn't match the model's flight day

    Returns
    -------
    pandas.DataFrame
        One row per obs point
        Model gases in mol/mol are auto-scaled to ppb.
    """
    
    obs_scale = obs_scale or {"NO2": 0.001}

    # model ncol to time, (time, lev)
    m = xr.open_dataset(model_file, decode_times=True).swap_dims({"ncol": "time"}).sortby("time")
    pmid = m["PMID"].values                          # (time, lev) Pa
    #mtime = pd.to_datetime(m["time"].values)
    mtime = _track_times(m)

    # obs mask fills before resample/scale
    odf = xr.open_dataset(obs_file).to_dataframe().reset_index()
    odf["time"] = pd.to_datetime(odf["time"])
    for c in odf.select_dtypes("number").columns:
        odf.loc[odf[c] <= fill_below, c] = np.nan
    odf = (odf.set_index("time").resample(resample).mean(numeric_only=True)
              .dropna(subset=["pressure_obs"]).reset_index())
    for v, s in obs_scale.items():
        if v in odf:
            odf[v] = odf[v] * s

    # nearest model track sample per obs time
    mi = mtime.values.astype("datetime64[ns]").astype("int64")
    oi = odf["time"].values.astype("datetime64[ns]").astype("int64")
    pos = np.clip(np.searchsorted(mi, oi), 1, len(mi) - 1)
    j = np.where((oi - mi[pos - 1]) <= (mi[pos] - oi), pos - 1, pos)

    gap = np.abs(mi[j] - oi).max() / 1e9
    if gap > gap_warn_hours * 3600:
        print(f"WARNING: max obs-model time gap = {gap / 3600:.1f} h "
              f"-- wrong obs file for this model day?")

    out = {"time": odf["time"].values,
           "latitude": odf.get("latitude", odf.get("lat")).values,
           "longitude": odf.get("longitude", odf.get("lon")).values,
           "altitude": odf["altitude"].values,
           "pressure_obs": odf["pressure_obs"].values}

    for mod_v, obs_v in mapping.items():
        if mod_v not in m:
            continue
        da = m[mod_v]
        scale = 1e9 if "mol/mol" in da.attrs.get("units", "").lower() else 1.0  # -> ppb
        prof = da.values * scale
        mvals = np.full(len(odf), np.nan)
        for i, (jj, p) in enumerate(zip(j, odf["pressure_obs"].values)):
            order = np.argsort(pmid[jj])
            mvals[i] = np.interp(p, pmid[jj][order], prof[jj][order],
                                 left=np.nan, right=np.nan)
        out[obs_v] = odf[obs_v].values if obs_v in odf else np.nan
        out[f"{obs_v}_mod"] = mvals
    return pd.DataFrame(out)

def save_paired_mm(df, mapping, obs_label, model_label, out_nc):
    """Write track-paired DataFrame to be compatible with MM's aircraft paired-file 

    Produces dims (time, x) 
    
    read_analysis can reconstruct the pair
    """
    
    obs_vars = list(mapping.values())
    model_vars = list(mapping.keys())
    obs_names = set(obs_vars)
    n = len(df)
    dims = ("time", "x")
    f64 = lambda a: (dims, np.asarray(a, float).reshape(n, 1))
    f32 = lambda a: (dims, np.asarray(a, np.float32).reshape(n, 1))

    ds = xr.Dataset()
    ds["latitude"] = f64(df["latitude"])
    ds["longitude"] = f64(df["longitude"])
    ds["pressure_obs"] = f32(df["pressure_obs"])
    ds["altitude"] = f64(df["altitude"])

    for mod_v, obs_v in mapping.items():
        if obs_v in df:
            ds[obs_v] = f64(df[obs_v])
        mcol = f"{mod_v}_new" if mod_v in obs_names else mod_v
        ds[mcol] = f32(df[f"{obs_v}_mod"])

    t = pd.to_datetime(df["time"].values)
    ds = ds.assign_coords(time=("time", t))
    t0 = t[0]
    ds["time"].encoding = {"units": f"minutes since {t0:%Y-%m-%d %H:%M:%S}",
                           "calendar": "proleptic_gregorian", "dtype": "int64"}

    meta = {"type": "aircraft", "radius_of_influence": None,
            "obs": obs_label, "model": model_label,
            "model_vars": model_vars, "obs_vars": obs_vars,
            "filename": f"{obs_label}_{model_label}.nc"}
    ds.attrs = {"title": "", "format": "NetCDF-4",
                "dict_json": json.dumps(meta, indent=4),
                "group_name": f"{obs_label}_{model_label}"}
    ds.to_netcdf(out_nc)
    print("wrote", out_nc)

# helper function to help with curtain plot 

def build_track_curtain(model_file, mod_var, times, num_levels=100, to_ppb=True,
                        vert_coord="pressure", temp_var="T", phis_var="PHIS",
                        ps_var="PS"):
    """Build (target_pressures, model_data_2d) for a curtain plot from a track file

    For each obs time, the nearest model track sample's column is vertically
    interpolated onto a regular pressure grid — giving the model's time x
    pressure curtain along the flight, the same array shape MM's
    ``make_curtain_plot`` expects (``(n_time, n_levels)``).

    Parameters
    ----------
    model_file : str
        Along-track '1s_1pt' model file.
    mod_var : str
        **Raw** model variable name in the file (e.g. 'O3', not 'O3_new').
    times : array-like
        Obs times (the paired-data times) to align the curtain columns to.
    num_levels : int
        Number of vertical interpolation levels (curtain y-resolution).
    to_ppb : bool
        Scale mol/mol gases to ppb (default True).
    vert_coord : 
        {'pressure', 'altitude'}
    temp_var : str ; temp K 
    phis_var : str ; sfc geopoetntial 
    ps_var : str ; sfc pressure 

    temp, phis, ps req when vert_coord = altitude 

    Returns
    -------
    target_pressures : numpy.ndarray
        (num_levels,) Pa, descending (surface high p ; top low p).
    model_data_2d : numpy.ndarray
        (len(times), num_levels) model field for the contourf.
    """
    m = xr.open_dataset(model_file, decode_times=True).swap_dims({"ncol": "time"}).sortby("time")
    pmid = m["PMID"].values                          # (time, lev) Pa
    #mt = pd.to_datetime(m["time"].values)
    mt = _track_times(m)
    da = m[mod_var]
    scale = 1e9 if (to_ppb and "mol/mol" in da.attrs.get("units", "").lower()) else 1.0
    prof = da.values * scale                         # (time, lev)

    mi = mt.values.astype("datetime64[ns]").astype("int64")
    oi = pd.to_datetime(times).values.astype("datetime64[ns]").astype("int64")
    pos = np.clip(np.searchsorted(mi, oi), 1, len(mi) - 1)
    j = np.where((oi - mi[pos - 1]) <= (mi[pos] - oi), pos - 1, pos)

    if vert_coord == "altitude":
        if temp_var not in m:
            raise ValueError(
                f"vert_coord='altitude' needs temperature var '{temp_var}' in {model_file}"
            )
        t_all = m[temp_var].values                   # (time, lev) K
        zsurf = (m[phis_var].values / G0) if phis_var in m else np.zeros(len(mt))
        psurf = m[ps_var].values if ps_var in m else None
        
        # geopotential height (m) for each matched column
        vert = np.empty((len(oi), pmid.shape[1]))
        for i, jj in enumerate(j):
            zs = float(np.ravel(zsurf)[jj]) if np.ndim(zsurf) else float(zsurf)
            ps = float(psurf[jj]) if psurf is not None else None
            vert[i] = _column_altitude(pmid[jj], t_all[jj], zs, ps)
        target_levels = np.linspace(np.nanmin(vert), np.nanmax(vert), num_levels)
        model_data_2d = np.full((len(oi), num_levels), np.nan)
        for i, jj in enumerate(j):
            order = np.argsort(vert[i])
            model_data_2d[i] = np.interp(target_levels, vert[i][order], prof[jj][order],
                                         left=np.nan, right=np.nan)
        return target_levels, model_data_2d
            
    target_pressures = np.linspace(np.nanmax(pmid), np.nanmin(pmid), num_levels)
    model_data_2d = np.full((len(oi), num_levels), np.nan)
    for i, jj in enumerate(j):
        order = np.argsort(pmid[jj])
        model_data_2d[i] = np.interp(target_pressures, pmid[jj][order], prof[jj][order],
                                     left=np.nan, right=np.nan)
    return target_pressures, model_data_2d

# altitude column using hypso metric eqn 
def _column_altitude(pmid_col, t_col, z_surf=0.0, p_surf=None):
    """Geopotential height (m) at each model midpoint 

    Integrates the hypsometric equation upward from the surface:
    dz = (Rd * Tv / g) * ln(p_below / p) 
    
    ASSUME: Tv ~ T 
    
    Parameters
    ----------
    pmid_col, t_col : array-like
        1D mid-layer pressure (Pa) and temperature (K) for one column,
        in any vertical ordering.
    z_surf : float
        Surface geopotential height (m), e.g. ``PHIS / g``. Default 0.
    p_surf : float, optional
        Surface pressure (Pa). If None, the lowest mid-layer is used as the
        anchor (so its height == ``z_surf``).

    Returns
    -------
    numpy.ndarray
        Height (m) at each level
    """
    
    pmid_col = np.asarray(pmid_col, float)
    t_col = np.asarray(t_col, float)
    order = np.argsort(pmid_col)[::-1]          # surface (high p) ; top (low p)
    p = pmid_col[order]
    tv = t_col[order]
    z = np.empty_like(p)
    p_below = p_surf if p_surf is not None else p[0]
    z_below = z_surf
    for k in range(len(p)):
        z[k] = z_below + (RD * tv[k] / G0) * np.log(p_below / p[k])
        z_below, p_below = z[k], p[k]
    out = np.empty_like(z)
    out[order] = z
    return out


