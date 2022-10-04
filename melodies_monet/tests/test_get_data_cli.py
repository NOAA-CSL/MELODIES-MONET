"""
Check for consistency with the tutorial datasets and that options work.
"""
import subprocess

import numpy as np
import xarray as xr

from melodies_monet.tutorial import fetch_example

ds0_aeronet = xr.open_dataset(fetch_example("aeronet:2019-09"))


def test_get_aeronet(tmp_path):
    fn = "x.nc"
    cmd = [
        "melodies-monet", "get-aeronet",
        "-s", "2019-09-01", "-e", "2019-09-02",
        "--dst", tmp_path.as_posix(), "-o", fn,
        "--no-compress",
        # TODO: compress doesn't work with current xarray (starting with v2022.06)
    ]
    subprocess.run(cmd, check=True)

    # To compare, need to have site ID as dim, not made-up integer dim 'x'
    # since positions may differ due to NaN-lat/lon dropping or such
    ds = xr.open_dataset(tmp_path / fn).squeeze().swap_dims(x="siteid")
    ds0 = ds0_aeronet.sel(time=ds.time).squeeze().swap_dims(x="siteid")
    # TODO: seems original loading missing value as -1 (on purpose, due to compress routine)
    
    assert not ds.identical(ds0)
    assert ds.time.equals(ds0.time)
    # assert (np.abs(ds.aod_551nm - ds0.aod_551nm) < 1e-9).all()
    assert (np.abs(ds.aod_551nm - ds0.aod_551nm).to_series().dropna() < 1e-9).all()
    # - Many more site IDs in ds0 (400 vs 283), and one that is in ds but not ds0
    # - In the above, only two sites
