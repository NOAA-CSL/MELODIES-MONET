"""
Need a way to compare different obs simultaneously rather than one at a time. Because of how involved this may become, i wrote a septerate 
obs2obs util folder to keep any issues with it seperate (as much as possible) from the rest of the MM 

Goal: Use this file to put diverse obs in mm plotting compliant format using the saved paired dfs. As many plots as possible should route throug
normal satplots or sfcplots to reduce upkeep

obs2obs: compare multiple obs AND multiple models simultaneously.

MELODIES-MONET pairing is one-obs to N-models; this util works AFTER pairing,
on the saved pair objects (an.paired), treating each (pair label, variable) as
an independent *data series*. Because the satellite obs-grid pairs share one
lat/lon grid (and model-space pairs share the model mesh)

multi-obs comparison should be pure alignment -- no new pairing needed 

Series spec (each entry in a group's ``series`` list):
    - {label: <pair label>, var: <variable>, name: <legend>, color: <opt>}
      raw values ("value" mode). Only mix series with the SAME units.
    - {label: ..., obs_var: X, mod_var: Y, name: ..., mode: bias}
      model - obs from within one pair.
    - {label: ..., obs_var: X, mod_var: Y, name: ..., mode: relative_bias}
      100*(model-obs)/obs. Dimensionless units will allow to put surface ppb,
      satellite columns, and aircraft data on one axis.

Handles pair objects that are xarray Datasets (satellite: (time,y,x),
(time,n_face), or (time,) series) and pandas DataFrames (surface/aircraft
point obs). Time matching: 'none', 'daily', or {window_local: [h0, h1]}
(solar local hour ~ UTC + lon/15; the honest way to co-locate a
geostationary full-day product with a ~13:30 LT polar orbiter).

Driver file: 

    an.read_control(); an.read_analysis()
    from melodies_monet.util import obs2obs_util
    obs2obs_util.run(an.paired, an.control_dict["obs2obs"],
                     default_outdir=an.output_dir)
                     
"""

import os
import gc

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# use native plotting capabilities within MM 
from melodies_monet.plots import surfplots as splots
from melodies_monet.plots import xarray_plots as xrplots

from melodies_monet.plots import savefig as _mm_savefig

def _save(outname_png):
    """Native MM save (adds the logo), same call pattern as xarray_plots."""
    try:
        _mm_savefig(outname_png, logo_height=100, bbox_inches="tight", dpi=200)
    except Exception as e:  # never lose a figure over logo trouble
        print(f"obs2obs: mm savefig failed ({e!r}); plain matplotlib save.",
              flush=True)
        plt.savefig(outname_png, dpi=200, bbox_inches="tight")

def _fnum(x):
    """PyYAML parses '5.0e16' (no '+') as a STRING, which then compares
    silently-wrong; coerce every user-supplied numeric defensively."""
    return None if x is None else float(x)

def _clip2(clip):
    """clip: [lo, hi] from YAML (float, float) or None.
       Add an additioanl way for user to clip values they may think are extraneous 
    """
    return None if clip is None else (float(clip[0]), float(clip[1]))
    
def _pair_obj(paired, label):
    if label not in paired:
        raise KeyError(
            f"obs2obs: pair label {label!r} not in paired data. "
            f"Available: {sorted(paired)}. Check read.paired.filenames.")
    return paired[label].obj

def _as_dataframe(obj):
    """Surface/aircraft pairs may carry a (Multi)Index; normalize to columns."""
    return obj.reset_index() if not isinstance(obj.index, pd.RangeIndex) else obj

def _local_hour_mask_xr(da, dset, h0, h1):
    """(time, ...) mask where solar local hour (UTC + lon/15) is in [h0, h1]."""
    if "longitude" not in dset.coords and "longitude" not in dset.variables:
        print("obs2obs: window_local requested but no longitude on this pair; "
              "skipping the window for it.", flush=True)
        return None
    hours = da["time"].dt.hour + da["time"].dt.minute / 60.0
    local = (hours + dset["longitude"] / 15.0) % 24.0
    mask = (local >= h0) & (local <= h1)
    
    try:
        frac = float(mask.mean().values)
        mean_kept = float(local.where(mask).mean().values)
        print(f"obs2obs: window_local [{h0},{h1}] LST -> kept {frac * 100:.0f}% "
              f"of (time x cell) samples, mean kept local hour {mean_kept:.1f} "
              "(TROPOMI should read ~13.5)", flush=True)
    except Exception:  # diagnostics must never sink a plot
        pass
    return mask

def _subsample(flat, max_points):
    """Deterministic strided subsample so month-scale swaths (~1e8 points)
    don't blow up the boxplot frame; quartiles are unaffected."""
    if max_points and flat.size > max_points:
        step = int(np.ceil(flat.size / max_points))
        return flat[::step]
    return flat

def _extract(spec, paired, time_match, max_points=2_000_000):
    """One series: dict(name, color, ts: pandas Series over time,
    flat: 1-D valid values, grid: time-mean (y,x) DataArray or None)."""
    label = spec["label"]
    mode = spec.get("mode", "value")
    name = spec.get("name", f"{label}:{spec.get('var', mode)}")
    obj = _pair_obj(paired, label)

    if isinstance(obj, xr.Dataset):
        if mode == "value":
            da = obj[spec["var"]]
        else:
            da_o, da_m = obj[spec["obs_var"]], obj[spec["mod_var"]]
            da = (da_m - da_o) if mode == "bias" \
                else 100.0 * (da_m - da_o) / da_o.where(da_o > 0)

        # temporary filter string
        clip = _clip2(spec.get("clip"))   # [lo, hi] physical range -> NaN outside
        if clip is not None:
            da = da.where((da >= clip[0]) & (da <= clip[1]))
            
        bounds = spec.get("bounds")
        if bounds is not None and "longitude" in obj.variables:
            w, e, s, n = bounds
            da = da.where((obj["longitude"] >= w) & (obj["longitude"] <= e)
                          & (obj["latitude"] >= s) & (obj["latitude"] <= n),
                          drop=True)

        if isinstance(time_match, dict) and "window_local" in time_match:
            h0, h1 = time_match["window_local"]
            m = _local_hour_mask_xr(da, obj, h0, h1)
            if m is not None:
                da = da.where(m)

        spatial = [d for d in da.dims if d != "time"]
        ts = da.mean(dim=spatial, skipna=True).to_series() if "time" in da.dims \
            else pd.Series(dtype=float)
        if time_match == "daily" and not ts.empty:
            ts = ts.resample("1D").mean()

        vals = da.values.ravel()
        flat = _subsample(vals[np.isfinite(vals)], max_points)

        grid = None
        if {"y", "x"} <= set(da.dims):
            grid = da.mean(dim="time", skipna=True) if "time" in da.dims else da

        # da is the masked / windowed array itself for the plots that need 
        # point level alignmet between timeseries (e.g. the scatter plots) rather than just a pre-reduced ts/flat/grid view
        return dict(name=name, color=spec.get("color"), ts=ts, flat=flat,
                    grid=grid, da = da, src=obj)

    # pandas path (surface / aircraft point pairs)
    df = _as_dataframe(obj)
    if mode == "value":
        v = df[spec["var"]]
    else:
        o, m = df[spec["obs_var"]], df[spec["mod_var"]]
        v = (m - o) if mode == "bias" else 100.0 * (m - o) / o.where(o > 0)
    
    clip = _clip2(spec.get("clip"))   # [lo, hi] physical range -> NaN outside
    if clip is not None:
        v = v.where((v >= clip[0]) & (v <= clip[1]))
        
    work = pd.DataFrame({"time": pd.to_datetime(df["time"]), "v": v})
    if "longitude" in df.columns:
        work["longitude"] = df["longitude"].values
        work["latitude"] = df["latitude"].values

    bounds = spec.get("bounds")
    if bounds is not None and "longitude" in work.columns:
        w, e, s, n = bounds
        work = work[(work.longitude >= w) & (work.longitude <= e)
                    & (work.latitude >= s) & (work.latitude <= n)]

    if isinstance(time_match, dict) and "window_local" in time_match:
        if "longitude" in work.columns:
            h0, h1 = time_match["window_local"]
            lh = (work.time.dt.hour + work.time.dt.minute / 60.0
                  + work.longitude / 15.0) % 24.0
            work = work[(lh >= h0) & (lh <= h1)]
        else:
            print(f"obs2obs: window_local skipped for {label} (no longitude).",
                  flush=True)

    ts = work.set_index("time")["v"].resample(
        "1D" if time_match == "daily" else "1h").mean().dropna()
    flat = work["v"].values
    flat = _subsample(flat[np.isfinite(flat)], max_points)
    return dict(name=name, color=spec.get("color"), ts=ts, flat=flat,
                grid=None, da=None, src=obj)


# MM plotting 

def _series_dataset(s):
    """1-var Dataset ('v') with a time dim for make_timeseries: prefer the
    masked (time, spatial) DataArray from _extract; fall back to the reduced
    domain-mean timeseries for point pairs that produced no 'da'."""
    da = s.get("da")
    if da is not None and "time" in getattr(da, "dims", ()):
        return da.to_dataset(name="v")
    ts = s.get("ts")
    if ts is None or ts.empty:
        return None
    idx = pd.DatetimeIndex(ts.index, name="time")
    return xr.DataArray(ts.to_numpy(dtype=float), dims=("time",),
                        coords={"time": idx}, name="v").to_dataset()

def multi_timeseries(series, outname, ylabel=None, title=None, avg_window=None,
                     vmin=None, vmax=None, set_axis=False, fig_dict=None,
                     text_dict=None):
    """Overlay every series on one axis via xrplots.make_timeseries."""

    tk = {"fontsize": 16, **(text_dict or {})}
    fig, ax = plt.subplots(figsize=(fig_dict or {}).get("figsize", (12, 6)))
    drawn = 0
    for s in series:
        dset = _series_dataset(s)
        if dset is None:
            print(f"obs2obs: no timeseries for {s['name']}; skipping.",
                  flush=True)
            continue
        # always pass our ax (overlay branch) and a non-None plot_dict
        xrplots.make_timeseries(
            dset, varname="v", label=s["name"], ax=ax, avg_window=avg_window,
            ylabel=ylabel, vmin=vmin, vmax=vmax,
            plot_dict=({"color": s["color"]} if s["color"] else {}),
            text_dict=text_dict, debug=False)
        drawn += 1
    if not drawn:
        plt.close(fig)
        return
    if set_axis:                       # keep the pinned vmin/vmax limits
        ax.set_ylim(vmin, vmax)
    else:                              # span every series (obs + all models)
        ax.relim()
        ax.autoscale(enable=True, axis="y")
    if title:
        ax.set_title(title, fontweight="bold", fontsize=tk["fontsize"])
    plt.tight_layout()
    _save(outname + ".png")
    plt.close(fig)
    print(f"obs2obs: wrote {outname}.png ({drawn} series)", flush=True)


def multi_boxplot(series, outname, ylabel=None, title=None, vmin=None,
                  vmax=None, fig_dict=None, text_dict=None, showfliers=False,
                        hline=None): # add outliers as an option in the yaml
    """One box per series via surfplots.make_boxplot.

    The combined frame is built with pd.concat (NOT column assignment as in
    calculate_boxplot) because our series have different lengths.
    """
    cols, label_bx = {}, []
    for s in series:
        if s["flat"].size == 0:
            print(f"obs2obs: no valid data for {s['name']}; skipping box.",
                  flush=True)
            continue
        cols[s["name"]] = pd.Series(s["flat"])
        label_bx.append(dict(color=s["color"] or "gray", linestyle="-",
                             marker="x", linewidth=1.2, markersize=6.0,
                             column=s["name"], label=s["name"]))
    if not cols:
        return
    comb_bx = pd.concat(cols, axis=1)   # union index; shorter series NaN-pad
    splots.make_boxplot(comb_bx, label_bx, ylabel=ylabel, vmin=vmin, vmax=vmax,
                        outname=outname, domain_type="all",
                        domain_name=title or "", fig_dict=fig_dict,
                        text_dict=text_dict, showfliers=showfliers, hline = hline, debug=False)
    print(f"obs2obs: wrote {outname} (boxplot)", flush=True)


def diff_map(series, outname, ylabel=None, title=None, vdiff=None,
             domain_name=None, fig_dict=None, text_dict=None):
    """series[0] - series[1] mean map via xrplots.make_spatial_bias_gridded."""
    a, b = series[0], series[1]
    if a["grid"] is None or b["grid"] is None:
        print("obs2obs: diff_map needs two gridded (y,x) series "
              f"({a['name']}, {b['name']}); skipping.", flush=True)
        return
    if a["grid"].shape != b["grid"].shape:
        print(f"obs2obs: diff_map grids differ {a['grid'].shape} vs "
              f"{b['grid'].shape}; regrid to a common obs grid first.",
              flush=True)
        return
    ds = xr.Dataset(
        {"s_minuend": a["grid"], "s_subtrahend": b["grid"]},
        coords={"longitude": a["src"]["longitude"],
                "latitude": a["src"]["latitude"]},
    )

    if "time" in a["src"].variables or "time" in a["src"].coords:
        _t = np.asarray(a["src"]["time"].values)
        if _t.size:
            ds = ds.assign_coords(time=("time", [_t.min(), _t.max()]))
            
    # plots s_minuend - s_subtrahend; the auto-region branch zooms to the
    # finite footprint of the difference.
    xrplots.make_spatial_bias_gridded(
        ds, varname_o="s_subtrahend", label_o=b["name"],
        varname_m="s_minuend", label_m=a["name"],
        ylabel=ylabel or f"{a['name']} - {b['name']}", vdiff=vdiff,
        outname=outname, domain_type="auto-region:box",
        domain_name=domain_name or (title or "obs2obs"),
        fig_dict=fig_dict or {"states": True, "counties": False},
        text_dict=text_dict, debug=False)
    print(f"obs2obs: wrote {outname} (diff_map)", flush=True)

def scatter(series, outname, ylabel=None, xlabel=None, title=None,
            fig_dict=None, text_dict=None, density=False, xylim=None,
            time_avg="1D", max_plot_points=200_000):
    """Point-matched scatter of series[1] (y) against series[0] (x).

    Alignment before plotting

    - two GRIDDED series (obsgrid pairs): daily-mean per cell (``time_avg``),
      then inner-join on the shared (time, y, x); every point is the same
      cell on the same day seen by both series. Requires the same grid --
      e.g. obs vs model within one pair (always same grid), or two pairs on
      the same obsgrid resolution/extent.
    - otherwise: daily domain-mean timeseries join (one point per day).

    Draws the 1:1 line and an N/r/slope/MB/RMSE stats box (computed on ALL
    matched points), then renders via surfplots.make_scatter_density_plot
    (best-fit line included; ``density: true`` -> seaborn KDE fill).
    
    """
    a, b = series[0], series[1]

    def _daily(s):
        da = s.get("da")
        if isinstance(da, xr.DataArray) and "time" in da.dims and time_avg:
            return da.resample(time=time_avg).mean(skipna=True)
        return da

    xa, ya = _daily(a), _daily(b)
    if isinstance(xa, xr.DataArray) and isinstance(ya, xr.DataArray):
        if (xa.sizes.get("y"), xa.sizes.get("x")) != (ya.sizes.get("y"), ya.sizes.get("x")):
            print(f"obs2obs: scatter {a['name']!r} vs {b['name']!r}: grids "
                  f"differ ({dict(xa.sizes)} vs {dict(ya.sizes)}). Cell-matched "
                  "scatter needs the same obsgrid -- pair both products at the "
                  "same res/extent, or compare within one pair (obs vs model).",
                  flush=True)
            return
        xa, ya = xr.align(xa, ya, join="inner")
        x, y = xa.values.ravel(), ya.values.ravel()
    else:
        ja = a["ts"].resample("1D").mean() if not a["ts"].empty else a["ts"]
        jb = b["ts"].resample("1D").mean() if not b["ts"].empty else b["ts"]
        j = pd.concat([ja, jb], axis=1, join="inner").dropna()
        if j.empty:
            print(f"obs2obs: scatter {a['name']!r} vs {b['name']!r}: no "
                  "overlapping days.", flush=True)
            return
        x, y = j.iloc[:, 0].to_numpy(float), j.iloc[:, 1].to_numpy(float)

    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size < 3:
        print(f"obs2obs: scatter {a['name']!r} vs {b['name']!r}: only "
              f"{x.size} matched points; skipping.", flush=True)
        return

    # stats on the FULL matched set (plot may be subsampled below)
    n = int(x.size)
    mb = float(np.mean(y - x))
    rmse = float(np.sqrt(np.mean((y - x) ** 2)))
    r = float(np.corrcoef(x, y)[0, 1])
    slope, intercept = (float(v) for v in np.polyfit(x, y, 1))
    xp, yp = _subsample(x, max_plot_points), _subsample(y, max_plot_points)

    if xylim is not None:
        lo, hi = float(xylim[0]), float(xylim[1])
    else:
        lo, hi = float(min(x.min(), y.min())), float(max(x.max(), y.max()))
        pad = 0.05 * (hi - lo or 1.0)
        lo, hi = lo - pad, hi + pad

    tk = {"fontsize": 14, **(text_dict or {})}
    fig, ax = plt.subplots(figsize=(fig_dict or {}).get("figsize", (7.5, 7)))
    ax.plot([lo, hi], [lo, hi], color="grey", lw=1.2, label="1:1")
    ax.text(0.97, 0.97,
            (f"N = {n}\nr = {r:.2f}\nslope = {slope:.2f}\n"
             f"MB = {mb:.2e}\nRMSE = {rmse:.2e}"),
            transform=ax.transAxes, va="top", ha="right",
            fontsize=tk["fontsize"] * 0.8,
            bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.85))

    # make_scatter_density_plot's scatter branch parses "(units)" out of the
    # ylabel for its colorbar -- guarantee the parentheses exist.
    yl = ylabel or b["name"]
    units = yl[yl.find("(") + 1: yl.find(")")] if "(" in yl else "value"
    if "(" not in yl:
        yl = f"{b['name']} ({units})"
    xl = xlabel or (f"{a['name']} ({units})" if "(" not in (xlabel or a["name"])
                    else a["name"])

    df = pd.DataFrame({"x": xp, "y": yp})
    splots.make_scatter_density_plot(
        df, mod_var="x", obs_var="y", ax=ax, fill=bool(density),
        xlabel=xl, ylabel=yl, title=title,
        vmin_x=lo, vmax_x=hi, vmin_y=lo, vmax_y=hi,
        outname=outname + ".png")
    print(f"obs2obs: wrote {outname}.png (scatter, N={n})", flush=True)

def _series_points(spec, paired, tm, max_points=2_000_000):
    """Tidy per-record DataFrame [obs, mod, longitude, latitude, time,
    local_hour] for one obs_var/mod_var series -- works for both gridded
    (satellite obsgrid, flattened) and point (surface) pairs, with clip,
    bounds, and the local-hour window already applied and NaN rows dropped.

    The common currency for the platform-aware plot types (taylor, diurnal,
    xplatform_scatter): each needs obs and model side-by-side per record,
    not the pre-reduced ts/flat/grid of _extract.
    """
    obj = _pair_obj(paired, spec["label"])
    # value-mode series (single ``var``, as in the emissions boxplots) map
    # both obs and mod to that one field, so callers can use df["obs"].
    if "var" in spec and "mod_var" not in spec:
        ov = mv = spec["var"]
    else:
        ov, mv = spec["obs_var"], spec["mod_var"]
        
    if isinstance(obj, xr.Dataset):
        # build via numpy ravel 
        o = obj[ov]
        m = obj[mv].broadcast_like(o) if obj[mv].dims != o.dims else obj[mv]
        # cols = {"obs": np.asarray(o.values).ravel(),
        #         "mod": np.asarray(m.values).ravel()}

        # CROP to bounds in xarray FIRST (drop=True shrinks the array to the
        # box). build columns by boolean mask 
        bnds = spec.get("bounds")
        inbox = None
        if bnds is not None and ("longitude" in obj.coords
                                 or "longitude" in obj.variables):
            w, e, s, n = bnds
            inbox = ((obj["longitude"] >= w) & (obj["longitude"] <= e)
                     & (obj["latitude"] >= s) & (obj["latitude"] <= n))
            o = o.where(inbox, drop=True)
            m = m.where(inbox, drop=True)
        ov_vals, mv_vals = np.asarray(o.values), np.asarray(m.values)
        keep = np.isfinite(ov_vals) & np.isfinite(mv_vals)   # mask before cols
        cols = {"obs": ov_vals[keep], "mod": mv_vals[keep]}
        
        for c in ("longitude", "latitude"):
            if c in obj.coords or c in obj.variables:
                # cols[c] = np.broadcast_to(
                #     obj[c].broadcast_like(o).values, o.shape)[keep]
                # crop the coord with the SAME inbox so it aligns with the
                # cropped obj
                cc = obj[c].where(inbox, drop=True) if inbox is not None else obj[c]
                cols[c] = np.broadcast_to(cc.broadcast_like(o).values, o.shape)[keep]
                
        if "time" in o.coords or "time" in obj.coords:
            tt = o["time"] if "time" in o.coords else obj["time"]
            cols["time"] = np.broadcast_to(tt.broadcast_like(o).values, o.shape)[keep]
        df = pd.DataFrame(cols)
    else:
        df = _as_dataframe(obj).rename(columns={ov: "obs", mv: "mod"})
        keep = ["obs", "mod"] + [c for c in ("longitude", "latitude", "time")
                                 if c in df.columns]
        df = df[keep].copy()

    # drop non-finite obs/mod BEFORE the rest -> the frame is small from here on
    df = df[np.isfinite(df["obs"]) & np.isfinite(df["mod"])]
    if "time" in df.columns:
        df["time"] = pd.to_datetime(df["time"])
        
    clip = _clip2(spec.get("clip"))
    if clip is not None:
        df = df[df["obs"].between(*clip) & df["mod"].between(*clip)]

    bnds = spec.get("bounds")
    if bnds is not None and "longitude" in df.columns:
        w, e, s, n = bnds
        df = df[(df.longitude >= w) & (df.longitude <= e)
                & (df.latitude >= s) & (df.latitude <= n)]

    if "longitude" in df.columns and "time" in df.columns:
        df["local_hour"] = ((df.time.dt.hour + df.time.dt.minute / 60.0
                             + df.longitude / 15.0) % 24.0)
        if isinstance(tm, dict) and "window_local" in tm:
            h0, h1 = tm["window_local"]
            df = df[(df.local_hour >= h0) & (df.local_hour <= h1)]

    df = df.dropna(subset=["obs", "mod"])
    if max_points and len(df) > max_points:
        df = df.iloc[:: int(np.ceil(len(df) / max_points))]
    return df.reset_index(drop=True)

def taylor(specs, paired, tm, outname, title=None, ylabel=None, text_dict=None,
           fig_dict=None, ty_scale=1.6, max_points=2_000_000):
    """Normalized Taylor diagram: one point per platform (surface/TEMPO/
    TROPOMI). Each platform's std is normalized by its OWN obs std, so the
    reference std is 1 for all and different-unit platforms share one plot.
    """
    _markers = ["o", "s", "^", "D", "v", "P", "*"]
    tk = {"fontsize": 14, **(text_dict or {})}
    dia, drawn = None, 0
    
    for i, spec in enumerate(specs):
        try:
            df = _series_points(spec, paired, tm, max_points=max_points)
        except KeyError as e:
            print(f"obs2obs taylor: {e}; skipping.", flush=True)
            continue
        o, m = df["obs"].to_numpy(float), df["mod"].to_numpy(float)
        if o.size < 3 or o.std() == 0:
            print(f"obs2obs taylor: {spec.get('name', spec['label'])}: "
                  "too few points / zero variance; skipping.", flush=True)
            continue
            
        # 1-D record Dataset; make_taylor stacks all dims + dropna internally.
        dset = xr.Dataset({"obs": ("rec", o), "mod": ("rec", m)})
        dia = xrplots.make_taylor(
            dset, varname_o="obs", varname_m="mod",
            label_o="obs (ref)", label_m=spec.get("name", spec["label"]),
            dia=dia, normalize=True, ty_scale=ty_scale, ylabel=ylabel or "",
            plot_dict={"marker": _markers[i % len(_markers)], "ms": 11,
                       "ls": "", "color": spec.get("color")},
            fig_dict=fig_dict, text_dict=text_dict, debug=False)
        drawn += 1
    if not drawn or dia is None:
        plt.close("all")
        return

    if title:
        plt.title(title, fontweight="bold", fontsize=tk["fontsize"])
    _save(outname + ".png")
    plt.close(all)
    print(f"obs2obs: wrote {outname}.png (taylor, {drawn} platforms)", flush=True)

def diurnal(specs, paired, tm, outname, ylabel=None, title=None,
            normalize=False, range_shading="IQR", time_offset=0,
            vmin=None, vmax=None, text_dict=None, fig_dict=None, set_axis = False,
            max_points=2_000_000):
    """Mean value vs local solar hour, obs and model, one line pair per
    series. ``normalize: true`` divides each curve by its own 24-h mean so
    different-unit platforms (surface ppb vs column) can share the axis as
    diurnal SHAPES. The window in ``tm`` is ignored here , so pass time_match: none/daily at the group.
    """
    tk = {"fontsize": 14, **(text_dict or {})}
    fig, ax = plt.subplots(figsize=(fig_dict or {}).get("figsize", (10, 6)))
    drawn = 0
    for spec in specs:
        try:
            e = _extract(spec, paired, "none", max_points=max_points)  # no window
        except KeyError as ex:
            print(f"obs2obs diurnal: {ex}; skipping.", flush=True)
            continue
        da = e.get("da")
        if da is None or "time" not in getattr(da, "dims", ()):
            print(f"obs2obs diurnal: {e['name']} has no time dimension; "
                  "skipping.", flush=True)
            continue
        if normalize:
            mn = float(da.mean(skipna=True))
            if mn:
                da = da / mn
        var = spec.get("var", "value")
        xrplots.make_diurnal_cycle(
            da.to_dataset(name=var), var, ax=ax, time_offset=time_offset,
            range_shading=range_shading, label=e["name"],
            plot_dict=({"color": e["color"]} if e["color"] else {}),
            ylabel=ylabel, vmin=vmin, vmax=vmax,
            text_dict=(text_dict or {}), fig_dict=(fig_dict or {}), debug=False)
        drawn += 1
    if not drawn:
        plt.close(fig)
        return
    # Axis control (MM set_axis). make_diurnal_cycle pins ylim per call to the
    # sparse/noisy bands (e.g. the isolated 00Z hour), which stretches the axis.
    if set_axis:                       # pinned limits (vmin/vmax)
        ax.set_ylim(vmin, vmax)
    else:                              # autoscale to the CENTER LINES only:
        ax.relim()                     # relim ignores fill_between (Collections)
        ax.autoscale(enable=True, axis="y")   # so noisy std bands spill off-axis
    if title:
        ax.set_title(title, fontweight="bold", fontsize=tk["fontsize"])
    plt.tight_layout()
    _save(outname + ".png")
    plt.close(fig)
    print(f"obs2obs: wrote {outname}.png (diurnal, {drawn} series)", flush=True)
    
def _coupling_block(block, gbounds=None, gclip=None):
    
    """Normalize a surface:/column: config into (obs_spec, [model_specs]).

    Accepts three forms:
      - {label, obs_var, mod_var}             single run (obs + its one model)
      - {obs: {...}, models: [{...}, ...]}    explicit obs + N model runs
      - [obs_entry, model_entry, ...]         list; the FIRST entry is the obs
    Every returned spec is value-mode (carries 'var'), so _series_points yields
    the chosen field directly in df['obs']. Group bounds/clip are pushed on.
    Model runs carry 'name' so surface and column can be matched across blocks.
    """
    def _vspec(d, default_name=None):
        s = dict(d)
        if "var" not in s:                       # accept obs_var/mod_var too
            s["var"] = s.get("obs_var") or s.get("mod_var")
        s.pop("obs_var", None)                   # force value-mode so
        s.pop("mod_var", None)                   # _series_points uses df['obs']
        s.setdefault("name", default_name or s.get("label"))
        return _with_group_defaults(s, gbounds, gclip)

    if isinstance(block, dict) and "obs" in block:            # explicit form
        obs = _vspec(block["obs"], "obs")
        mods = [_vspec(m) for m in block.get("models", [])]
    elif isinstance(block, dict):                             # single-run dict
        obs = _vspec({"label": block["label"], "var": block["obs_var"]}, "obs")
        mods = [_vspec({"label": block["label"], "var": block["mod_var"],
                        "name": block.get("name", "model"),
                        "color": block.get("color", "purple")})]
    else:                                                     # list; first=obs
        obs = _vspec(block[0], "obs")
        mods = [_vspec(m) for m in block[1:]]
    return obs, mods

def _colocate_fields(sfc_spec, col_spec, paired, tm, max_points, max_sep):
    """Daily, nearest-cell colocation of a surface field with a column field.
    Returns (sfc_vals, col_vals) numpy arrays -- one pair per matched site-day
    within ``max_sep`` degrees. Uses each spec's value field (df['obs'] in
    value mode); cKDTree matches each surface site to its nearest column cell.
    """
    from scipy.spatial import cKDTree
    try:
        sdf = _series_points(sfc_spec, paired, tm, max_points=max_points)
        cdf = _series_points(col_spec, paired, tm, max_points=max_points)
    except KeyError as e:
        print(f"obs2obs coupling: {e}; skipping.", flush=True)
        return np.array([]), np.array([])
    if (sdf.empty or cdf.empty or "longitude" not in cdf
            or "time" not in cdf or "time" not in sdf):
        return np.array([]), np.array([])
    sdf = sdf.assign(day=sdf.time.dt.floor("D"))
    cdf = cdf.assign(day=cdf.time.dt.floor("D"))
    s_day = (sdf.groupby(["day", "longitude", "latitude"])["obs"]
             .mean().reset_index())
    c_day = (cdf.groupby(["day", "longitude", "latitude"])["obs"]
             .mean().reset_index())
    xs, ys = [], []
    for day, cg in c_day.groupby("day"):
        sg = s_day[s_day.day == day]
        if sg.empty or cg.empty:
            continue
        tree = cKDTree(cg[["longitude", "latitude"]].to_numpy())
        d, idx = tree.query(sg[["longitude", "latitude"]].to_numpy(), k=1)
        keep = d < max_sep
        xs.append(sg["obs"].to_numpy()[keep])
        ys.append(cg["obs"].to_numpy()[idx][keep])
    if not xs:
        return np.array([]), np.array([])
    x, y = np.concatenate(xs), np.concatenate(ys)
    m = np.isfinite(x) & np.isfinite(y)
    return x[m], y[m]

def xplatform_scatter(grp, paired, tm, outname, text_dict=None, fig_dict=None,
                      max_points=2_000_000):
    """Surface<->column coupling scatter. x = surface conc, y = colocated
    column. One black cloud for OBS (surface obs vs column obs) plus one cloud
    per emissions run (surface model vs column model), each with a best-fit
    slope. A model that preserves the vertical coupling (mixing height +
    profile) reproduces the observed surface->column slope.

    surface:/column: accept a single {label,obs_var,mod_var} dict (one run) or
    the multi-run {obs:{...}, models:[...]} / list forms (see _coupling_block);
    model runs are matched between surface and column by 'name'.
    """
    max_sep = float(grp.get("max_sep_deg", 0.5))
    gbounds, gclip = grp.get("bounds"), grp.get("clip")
    s_obs, s_mods = _coupling_block(grp["surface"], gbounds, gclip)
    c_obs, c_mods = _coupling_block(grp["column"], gbounds, gclip)

    tk = {"fontsize": 14, **(text_dict or {})}
    fig, ax = plt.subplots(figsize=(fig_dict or {}).get("figsize", (8, 7)))

    def _cloud(x, y, color, label):
        if x.size < 3:
            print(f"obs2obs xplatform_scatter: <3 pairs for {label}; skip.",
                  flush=True)
            return 0
        ax.scatter(x, y, s=10, alpha=0.35, color=color, edgecolors="none")
        sl, ic = np.polyfit(x, y, 1)
        r = np.corrcoef(x, y)[0, 1]
        xs = np.linspace(np.nanmin(x), np.nanmax(x), 50)
        ax.plot(xs, sl * xs + ic, color=color, lw=2,
                label=f"{label}: slope={sl:.2e}, r={r:.2f}")
        return 1

    drawn = 0
    xo, yo = _colocate_fields(s_obs, c_obs, paired, tm, max_points, max_sep)
    drawn += _cloud(xo, yo, s_obs.get("color", "black"),
                    grp.get("obs_label", "observed"))
    c_by = {m.get("name"): m for m in c_mods}          # match runs by name
    for sm in s_mods:
        cm = c_by.get(sm.get("name"))
        if cm is None:
            print(f"obs2obs xplatform_scatter: no column run named "
                  f"{sm.get('name')!r}; skip.", flush=True)
            continue
        x, y = _colocate_fields(sm, cm, paired, tm, max_points, max_sep)
        drawn += _cloud(x, y, sm.get("color", "purple"), sm.get("name", "model"))

    if not drawn:
        plt.close(fig)
        print("obs2obs xplatform_scatter: nothing drawn.", flush=True)
        return
    ax.set_xlabel(grp.get("xlabel", "surface"), fontweight="bold",
                  fontsize=tk["fontsize"])
    ax.set_ylabel(grp.get("ylabel", "column"), fontweight="bold",
                  fontsize=tk["fontsize"])
    ax.legend(fontsize=tk["fontsize"] * 0.7)
    ax.set_title(grp.get("title", "surface-column coupling"),
                 fontweight="bold", fontsize=tk["fontsize"])
    plt.tight_layout()
    _save(outname + ".png")
    plt.close(fig)
    print(f"obs2obs: wrote {outname}.png (xplatform_scatter, {drawn} clouds)",
          flush=True)

def _norm_point_df(obj, var):
    """[longitude, latitude, var] rows, whether the pair is a pandas DataFrame
    or an xarray Dataset (netcdf-backed pairs load as the latter)."""
    if isinstance(obj, xr.Dataset):
        ds = obj[[var]]
        for c in ("longitude", "latitude"):
            if c in obj.coords or c in obj.variables:
                ds = ds.assign_coords({c: obj[c]})
        d = ds.to_dataframe().reset_index()
    else:
        d = _as_dataframe(obj)
    return d[["longitude", "latitude", var]].dropna()
    
def norm_spatial(grp, paired, tm, outname, text_dict=None, fig_dict=None):
    """Normalized multiscale map: column time-mean field + surface site
    time-means, each divided by its own domain mean (dimensionless ratio) on a
    shared 0-2 scale. One panel for OBS, then one panel per emissions run
    (surface/column matched by 'name'). Shows whether surface and column
    hotspots co-locate and whether each run reproduces both.

    surface:/column: use the same forms as xplatform_scatter (single dict or
    multi-run {obs, models}/list); see _coupling_block.
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    bnds = grp.get("bounds")
    tk = {"fontsize": 14, **(text_dict or {})}
    s_obs, s_mods = _coupling_block(grp["surface"], bnds, grp.get("clip"))
    c_obs, c_mods = _coupling_block(grp["column"], bnds, grp.get("clip"))
    c_by = {m.get("name"): m for m in c_mods}

    # panel list: obs first, then each model run present in BOTH blocks
    panels = [("observed", s_obs, c_obs)]
    for sm in s_mods:
        cm = c_by.get(sm.get("name"))
        if cm is not None:
            panels.append((sm.get("name", "model"), sm, cm))

    def _field(cspec):
        obj = _pair_obj(paired, cspec["label"])
        da = obj[cspec["var"]]
        da = da.mean("time", skipna=True) if "time" in da.dims else da
        mean = float(np.nanmean(da.values))
        return obj, (da / mean if mean else da)

    def _sites(sspec):
        d = _norm_point_df(_pair_obj(paired, sspec["label"]), sspec["var"])
        if bnds is not None:
            w, e, s, n = bnds
            d = d[(d.longitude >= w) & (d.longitude <= e)
                  & (d.latitude >= s) & (d.latitude <= n)]
        g = d.groupby(["longitude", "latitude"])[sspec["var"]].mean().reset_index()
        mean = g[sspec["var"]].mean()
        g["r"] = g[sspec["var"]] / mean if mean else g[sspec["var"]]
        return g

    n = len(panels)
    fig, axes = plt.subplots(
        1, n, figsize=(fig_dict or {}).get("figsize", (8 * n, 6)),
        subplot_kw={"projection": ccrs.PlateCarree()})
    axes = np.atleast_1d(axes)
    pm = None
    for ax, (ttl, sspec, cspec) in zip(axes, panels):
        try:
            cobj, fld = _field(cspec)
            sit = _sites(sspec)
        except KeyError as e:
            print(f"obs2obs norm_spatial: {e}; blank panel {ttl}.", flush=True)
            continue
        pm = ax.pcolormesh(cobj["longitude"], cobj["latitude"], fld,
                           vmin=0, vmax=2, cmap="RdBu_r", shading="auto",
                           transform=ccrs.PlateCarree())
        ax.scatter(sit.longitude, sit.latitude, c=sit["r"], vmin=0, vmax=2,
                   cmap="RdBu_r", s=60, edgecolors="k", linewidths=0.6,
                   transform=ccrs.PlateCarree(), zorder=5)
        ax.coastlines(linewidth=0.5)
        ax.add_feature(cfeature.STATES, linewidth=0.3)
        if bnds is not None:
            ax.set_extent(bnds, crs=ccrs.PlateCarree())
        ax.set_title(ttl, fontweight="bold", fontsize=tk["fontsize"])
    if pm is None:
        plt.close(fig)
        print("obs2obs norm_spatial: nothing drawn.", flush=True)
        return
    cb = fig.colorbar(pm, ax=list(axes), orientation="horizontal", shrink=0.6,
                      pad=0.05, extend="both")
    cb.set_label("value / domain mean (column field, surface dots)",
                 fontweight="bold", fontsize=tk["fontsize"])
    if grp.get("title"):
        fig.suptitle(grp["title"], fontweight="bold", fontsize=tk["fontsize"])
    _save(outname + ".png")
    plt.close(fig)
    print(f"obs2obs: wrote {outname}.png (norm_spatial, {n} panels)", flush=True)

def operator_effect(grp, paired, tm, outname, ylabel=None, title=None,
                    text_dict=None, fig_dict=None, max_points=2_000_000):
    """Averaging-kernel operator effect, stacked two-panel.

    TOP: the SAME model column seen through each instrument's AK, grouped by
    emissions run (box color = run, hatch = instrument). Because the model is
    identical, any spread between the instrument boxes is the operator alone --
    this is the AK-driven TEMPO-vs-TROPOMI gap made explicit.
    BOTTOM: the observed column for each instrument (one box each), showing how
    much the instruments themselves disagree, independent of the model.
    Panels share the y-axis so operator spread and instrument disagreement are
    directly comparable.

    Config:
      instruments: {tempo: {obs_var, mod_var}, tropomi: {obs_var, mod_var}}
      runs: [{name, color, tempo: <pair label>, tropomi: <pair label>}, ...]
    """
    import matplotlib.patches as mpatches
    instruments = grp["instruments"]
    runs = grp["runs"]
    bnds, clip = grp.get("bounds"), grp.get("clip")
    inst_names = list(instruments)
    HATCH = ["", "//", "..", "xx", "\\\\"]
    hatch_of = {ins: HATCH[i % len(HATCH)] for i, ins in enumerate(inst_names)}
    plot_cap = min(int(max_points), 200_000)

    top, obs = [], {}          # top: (run_name, color, inst, values); obs: inst->vals
    for run in runs:
        for ins in inst_names:
            lbl = run.get(ins)
            if not lbl:
                continue
            spec = {"label": lbl, "obs_var": instruments[ins]["obs_var"],
                    "mod_var": instruments[ins]["mod_var"],
                    "bounds": bnds, "clip": clip}
            try:
                df = _series_points(spec, paired, tm, max_points=max_points)
            except KeyError as e:
                print(f"obs2obs operator_effect: {e}; skip.", flush=True)
                continue
            if df.empty:
                continue
            top.append((run.get("name", lbl), run.get("color", "gray"), ins,
                        _subsample(df["mod"].to_numpy(float), plot_cap)))
            if ins not in obs:
                obs[ins] = _subsample(df["obs"].to_numpy(float), plot_cap)
    if not top:
        print("obs2obs operator_effect: no data.", flush=True)
        return

    tk = {"fontsize": 13, **(text_dict or {})}
    fig, (axm, axo) = plt.subplots(
        2, 1, sharey=True, figsize=(fig_dict or {}).get("figsize", (10, 9)),
        gridspec_kw={"height_ratios": [3, 1]})

    K = len(inst_names)
    r_index = {run.get("name"): i for i, run in enumerate(runs)}
    positions, boxvals, facecolors, hatches = [], [], [], []
    for (rn, color, ins, vals) in top:
        r = r_index.get(rn, 0)
        positions.append(r * (K + 1) + inst_names.index(ins))
        v = vals[np.isfinite(vals)]
        boxvals.append(v if v.size else np.array([np.nan]))
        facecolors.append(color)
        hatches.append(hatch_of[ins])
    bp = axm.boxplot(boxvals, positions=positions, widths=0.8,
                     patch_artist=True, showfliers=False)
    for patch, fc, ha in zip(bp["boxes"], facecolors, hatches):
        patch.set_facecolor(fc); patch.set_alpha(0.55); patch.set_hatch(ha)
    for med in bp["medians"]:
        med.set_color("k")
    axm.set_xticks([r * (K + 1) + (K - 1) / 2 for r in range(len(runs))])
    axm.set_xticklabels([run.get("name", "") for run in runs])
    axm.set_ylabel(ylabel or "model column", fontweight="bold",
                   fontsize=tk["fontsize"])
    axm.set_title(title or grp.get("title",
                  "operator effect: model column by instrument AK"),
                  fontweight="bold", fontsize=tk["fontsize"])
    handles = [mpatches.Patch(facecolor="0.8", hatch=hatch_of[i], label=i)
               for i in inst_names]
    axm.legend(handles=handles, title="instrument AK",
               fontsize=tk["fontsize"] * 0.8)

    present = [i for i in inst_names if i in obs]
    obp = axo.boxplot([obs[i][np.isfinite(obs[i])] for i in present],
                      positions=list(range(len(present))), widths=0.6,
                      patch_artist=True, showfliers=False)
    for patch, i in zip(obp["boxes"], present):
        patch.set_facecolor("0.7"); patch.set_hatch(hatch_of[i])
    axo.set_xticks(range(len(present)))
    axo.set_xticklabels(present)
    axo.set_ylabel("obs column", fontweight="bold", fontsize=tk["fontsize"])
    axo.set_title("observed column by instrument",
                  fontweight="bold", fontsize=tk["fontsize"] * 0.9)

    plt.tight_layout()
    _save(outname + ".png")
    plt.close(fig)
    print(f"obs2obs: wrote {outname}.png (operator_effect, "
          f"{len(top)} model boxes)", flush=True)

# run

_PLOTTERS = {"multi_boxplot": multi_boxplot,
             "multi_timeseries": multi_timeseries,
             "diff_map": diff_map,
             "scatter": scatter}

# plot types that take the whole group because they need their own obs/model alignment
_GROUP_PLOTTERS = {"taylor": taylor, "diurnal": diurnal,
                   "xplatform_scatter": xplatform_scatter,
                   "norm_spatial": norm_spatial,
                   "operator_effect": operator_effect}
 
def _with_group_defaults(spec, gbounds, gclip):
    """Push group-level bounds/clip onto a series spec unless it sets its own.
    Used for BOTH the per-series plotters and the group plotters (taylor,
    diurnal) so a group's ``bounds:`` actually subsets those too."""
    if gbounds is not None and "bounds" not in spec:
        spec = {**spec, "bounds": gbounds}
    if gclip is not None and "clip" not in spec:
        spec = {**spec, "clip": gclip}
    return spec

def _apply_common_mask(series, max_points=2_000_000):
    """Co-mask the value-mode series so obs and every model compare on IDENTICAL
    (time, cell) samples: a cell is kept only where EVERY gridded series is
    finite (obs NaN -> model NaN and vice versa). Rewrites each series' flat/ts/
    da from the intersected mask. Non-destructive -- operates on the extracted
    arrays only, never the paired files 
    
    Series without a gridded time-varying da are left untouched.
    """
    idx_da = [(i, s["da"]) for i, s in enumerate(series)
              if isinstance(s.get("da"), xr.DataArray) and "time" in s["da"].dims]
    if len(idx_da) < 2:
        return series
    # a shared CELL mask only makes sense within one grid. If the series span
    # different grids (surface vs satellite, TEMPO vs TROPOMI), an inner align
    # intersects to nothing and would blank every series -> skip instead.
    sizes0 = dict(idx_da[0][1].sizes)
    if any(dict(d.sizes) != sizes0 for _, d in idx_da):
        print("obs2obs: common_mask skipped -- series are on different grids "
              "(surface/TEMPO/TROPOMI); co-masking applies within one grid only.",
              flush=True)
        return series
    try:
        aligned = xr.align(*[d for _, d in idx_da], join="inner")
    except Exception as e:  # noqa: BLE001 -- coords must match to co-mask
        print(f"obs2obs: common_mask align failed ({e!r}); left unmasked.",
              flush=True)
        return series
    valid = aligned[0].notnull()
    for a in aligned[1:]:
        valid = valid & a.notnull()
    for (i, _), a in zip(idx_da, aligned):
        m = a.where(valid)
        spatial = [d for d in m.dims if d != "time"]
        vals = m.values.ravel()
        series[i]["da"] = m
        series[i]["flat"] = _subsample(vals[np.isfinite(vals)], max_points)
        series[i]["ts"] = m.mean(dim=spatial, skipna=True).to_series()
    return series
    
def labels_used(config, only=None):
    """Every pair label referenced by the obs2obs groups (optionally only the
    groups whose name contains ``only``). Use this to prune
    read.paired.filenames BEFORE read_analysis so a run opens only the files
    it needs 
    """
    used = set()
    for gname, grp in config.get("groups", {}).items():
        if only and only not in gname:
            continue
        for s in grp.get("series", []):
            if s.get("label"):
                used.add(s["label"])
        for key in ("surface", "column"):        # xplatform/norm group configs
            blk = grp.get(key)
            if isinstance(blk, dict):
                if blk.get("label"):
                    used.add(blk["label"])
                if isinstance(blk.get("obs"), dict) and blk["obs"].get("label"):
                    used.add(blk["obs"]["label"])           # {obs, models} form
                for m in blk.get("models", []):
                    if isinstance(m, dict) and m.get("label"):
                        used.add(m["label"])
            elif isinstance(blk, list):          # list-of-series form
                for s in blk:
                    if isinstance(s, dict) and s.get("label"):
                        used.add(s["label"])
        for run in grp.get("runs", []):           # operator_effect config
            if isinstance(run, dict):
                for k, v in run.items():
                    if k not in ("name", "color") and isinstance(v, str):
                        used.add(v)
    return used
    
def run(paired, config, default_outdir=".", only = None):
    """Execute every group in the ``obs2obs:`` config against loaded pairs."""
    outdir = os.path.expandvars(config.get("output_dir", default_outdir))
    os.makedirs(outdir, exist_ok=True)
    default_tm = config.get("time_match", "none")
    max_points = int(config.get("max_points", 2_000_000))

    for gname, grp in config.get("groups", {}).items():
        if only and only not in gname:
            continue
        dp = grp.get("data_proc") or {}
        tm = grp.get("time_match", default_tm)
        gbounds = grp.get("bounds")
        gclip = grp.get("clip")

        # group-level plot types (taylor/diurnal/xplatform_scatter/norm_spatial)
        # handle their own extraction/colocation from the raw group config.
        types = grp.get("type", "multi_boxplot")
        types = [types] if isinstance(types, str) else list(types)
        _group_types = [t for t in types if t in _GROUP_PLOTTERS]

        # taylor/diurnal take a series list -> push group bounds/clip onto each
        # (this is the domain fix: without it a group's bounds: was ignored).
        gseries = [_with_group_defaults(s, gbounds, gclip)
                   for s in grp.get("series", [])]

        for t in _group_types:
            out = os.path.join(outdir, f"obs2obs_{gname}.{t}")
            kw = dict(text_dict=grp.get("text_kwargs"),
                      fig_dict=grp.get("fig_kwargs"), max_points=max_points)
            try:
                if t == "taylor":
                    _GROUP_PLOTTERS[t](gseries, paired, tm, out,
                                       title=grp.get("title", gname),
                                       ylabel=grp.get("ylabel"), **kw)
                elif t == "diurnal":
                    _GROUP_PLOTTERS[t](gseries, paired, tm, out,
                                       ylabel=grp.get("ylabel"),
                                       title=grp.get("title", gname),
                                       normalize=dp.get("normalize",
                                                         grp.get("normalize", False)),
                                       range_shading=grp.get("range_shading", "IQR"),
                                       time_offset=grp.get("time_offset", 0),
                                       vmin=_fnum(grp.get("vmin")),
                                       vmax=_fnum(grp.get("vmax")),
                                       set_axis=dp.get("set_axis",
                                                       grp.get("set_axis", False)),
                                       **kw)
                elif t == "operator_effect":
                    _GROUP_PLOTTERS[t](grp, paired, tm, out,
                                       ylabel=grp.get("ylabel"),
                                       title=grp.get("title", gname), **kw)
                else:                       # xplatform_scatter / norm_spatial
                    _kw = {k: v for k, v in kw.items() if k != "max_points"}
                    if t == "xplatform_scatter":
                        _kw["max_points"] = max_points
                    _GROUP_PLOTTERS[t](grp, paired, tm, out, **_kw)
            except Exception as e:  # noqa: BLE001
                print(f"obs2obs: {gname}.{t} failed: {e!r}", flush=True)
        types = [t for t in types if t not in _GROUP_PLOTTERS]
        if not types:
            continue
            
        series = []
        for spec in grp["series"]:
            if gbounds is not None and "bounds" not in spec:
                spec = {**spec, "bounds": gbounds}
            if gclip is not None and "clip" not in spec:
                spec = {**spec, "clip": gclip}
            try:
                series.append(_extract(spec, paired, tm, max_points=max_points))
            except KeyError as e:
                print(f"obs2obs: {gname}: {e}; skipping series.", flush=True)
        if not series:
            continue
        # co-mask obs + all models to identical samples (opt-in per group)
        if dp.get("common_mask", grp.get("common_mask", False)):
            series = _apply_common_mask(series, max_points=max_points)
        for t in types:
            fn = _PLOTTERS.get(t)
            if fn is None:
                print(f"obs2obs: unknown type {t!r} "
                      f"(use {sorted(_PLOTTERS)}).", flush=True)
                continue
            out = os.path.join(outdir, f"obs2obs_{gname}.{t}")
            kwargs = dict(ylabel=grp.get("ylabel"),
                          title=grp.get("title", gname),
                          fig_dict=grp.get("fig_kwargs"),
                          text_dict=grp.get("text_kwargs"))
            if t == "multi_boxplot":
                kwargs.update(
                    vmin=_fnum(grp.get("vmin")), vmax=_fnum(grp.get("vmax")),
                    showfliers=dp.get("showfliers",
                                      grp.get("showfliers", False)),
                    hline=_fnum(dp.get("hline", grp.get("hline"))))
            if t == "multi_timeseries":
                kwargs["avg_window"] = grp.get("avg_window")
                kwargs["vmin"] = _fnum(grp.get("vmin"))
                kwargs["vmax"] = _fnum(grp.get("vmax"))
            if t == "diff_map":
                kwargs.update(vdiff=_fnum(grp.get("vdiff")),
                              domain_name=grp.get("domain_name"))
            if t == "scatter":
                kwargs.update(xlabel=grp.get("xlabel"),
                              xylim=grp.get("xylim"),
                              density=dp.get("density", grp.get("density", False)),
                              time_avg=grp.get("time_avg", "1D"))
                
            fn(series, out, **kwargs)
        # free this group's extracted arrays before moving on
        del series
        gc.collect()
