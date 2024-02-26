# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#
"""
Check for consistency with the tutorial datasets and that options work.
"""
import subprocess

import numpy as np
import xarray as xr

from melodies_monet.tutorial import fetch_example

ds0_aeronet = xr.open_dataset(fetch_example("aeronet:2019-09"))
ds0_airnow = xr.open_dataset(fetch_example("airnow:2019-09"))


def test_get_aeronet(tmp_path):
    fn = "x.nc"
    cmd = [
        "melodies-monet", "get-aeronet",
        "-s", "2019-09-01", "-e", "2019-09-02",
        "--dst", tmp_path.as_posix(), "-o", fn,
        "--no-compress",
    ]
    subprocess.run(cmd, check=True)

    # To compare, need to have site ID as dim, not made-up integer dim 'x'
    # since positions may differ due to NaN-lat/lon dropping or such
    ds = xr.open_dataset(tmp_path / fn).squeeze().swap_dims(x="siteid")
    ds0 = ds0_aeronet.sel(time=ds.time).squeeze().swap_dims(x="siteid")
    # NOTE: -1 in ds0 indicates missing value, due to compress routine
    
    assert not ds.identical(ds0)
    assert ds.time.equals(ds0.time)
    ds0["aod_551nm"] = ds0["aod_551nm"].where(ds0["aod_551nm"] != -1)
    assert (np.abs(ds.aod_551nm - ds0.aod_551nm).to_series().dropna() < 1e-9).all()
    # - Many more site IDs in ds0 (400 vs 283), and one that is in ds but not ds0
    # - In the above, only two sites


def test_get_airnow(tmp_path):
    fn = "x.nc"
    cmd = [
        "melodies-monet", "get-airnow",
        "-s", "2019-09-01", "-e", "2019-09-02",
        "--dst", tmp_path.as_posix(), "-o", fn,
        "--no-compress",
    ]
    subprocess.run(cmd, check=True)

    ds = xr.open_dataset(tmp_path / fn).squeeze().swap_dims(x="siteid")
    ds0 = ds0_airnow.sel(time=ds.time).squeeze().swap_dims(x="siteid")

    assert ds.time.equals(ds0.time)

    for vn in ["NO2", "OZONE", "PM2.5"]:
        ds0[vn] = ds0[vn].where(ds0[vn] != -1)
        ds[vn] = ds[vn].where(~ ((ds[vn] == 0) & (ds0[vn] != 0)))
        assert (np.abs((ds[vn] - ds0[vn]) / ds0[vn]).to_series().dropna() < 2e-6).all()
        assert (np.abs(ds[vn] - ds0[vn]).to_series().dropna() < 3e-7).all()


def test_get_airnow_comp(tmp_path):
    fn = "x.nc"
    cmd = [
        "melodies-monet", "get-airnow",
        "-s", "2019-09-01", "-e", "2019-09-02",
        "--dst", tmp_path.as_posix(), "-o", fn,
        "--compress",
    ]
    subprocess.run(cmd, check=True)

    ds = xr.open_dataset(tmp_path / fn).squeeze().swap_dims(x="siteid")
    ds0 = ds0_airnow.sel(time=ds.time).squeeze().swap_dims(x="siteid")

    assert ds.time.equals(ds0.time)

    for vn in ["NO2", "OZONE", "PM2.5"]:
        ds0[vn] = ds0[vn].where(ds0[vn] != -1)
        ds[vn] = ds[vn].where(ds[vn] != -1)
        ds[vn] = ds[vn].where(~ ((ds[vn] == 0) & (ds0[vn] != 0)))
        # assert (np.abs((ds[vn] - ds0[vn]) / ds0[vn]).to_series().dropna() < 2e-6).all()
        assert (np.abs(ds[vn] - ds0[vn]).to_series().dropna() < 3e-7).all()


def test_get_ish_lite_box(tmp_path):
    fn = "x.nc"
    cmd = [
        "melodies-monet", "get-ish-lite",
        "-s", "2023-01-01", "-e", "2023-01-01 23:00",
        "--box", "39.5", "-105.75", "40.5", "-104.75",
        "--dst", tmp_path.as_posix(), "-o", fn,
    ]
    subprocess.run(cmd, check=True)

    ds = xr.open_dataset(tmp_path / fn)

    assert ds.time.size == 24
    assert np.unique(ds.state) == ["CO"]


def test_get_ish_box(tmp_path):
    fn = "x.nc"
    cmd = [
        "melodies-monet", "get-ish",
        "-s", "2023-01-01", "-e", "2023-01-01 23:00",
        "--box", "39.5", "-105.75", "40.5", "-104.75",
        "--dst", tmp_path.as_posix(), "-o", fn,
    ]
    subprocess.run(cmd, check=True)

    ds = xr.open_dataset(tmp_path / fn)

    assert ds.time.size == 24
    assert np.unique(ds.state) == ["CO"]
