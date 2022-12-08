# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#

"""
file: sparse_grid.py

test grid_util.update_sparse_data_grid

requires obs datafiles generated with, e.g.
    python setup_obs.py --nfile 3 --control sparse_grid.yaml
"""

import numpy as np
import pandas as pd
import xarray as xr

from glob import glob
from melodies_monet import driver

from melodies_monet.util import grid_util

an = driver.analysis()

an.control = 'sparse_grid.yaml'
an.read_control()


# lines below emulate functionality that would be in the driver

start_time = pd.to_datetime(an.control_dict['analysis']['start_time'])
end_time = pd.to_datetime(an.control_dict['analysis']['end_time'])
start_timestamp = start_time.timestamp()
end_timestamp = end_time.timestamp()

ntime = an.control_dict['obs']['test_obs']['sparse_data_grid']['ntime']
nlat = an.control_dict['obs']['test_obs']['sparse_data_grid']['nlat']
nlon = an.control_dict['obs']['test_obs']['sparse_data_grid']['nlon']
lon0 = an.control_dict['obs']['test_obs']['sparse_data_grid']['lon0']

# generate uniform grid
time_edges = np.linspace(start_timestamp, end_timestamp, ntime+1, endpoint=True, dtype=float)
time_grid = 0.5 * (time_edges[0:ntime] + time_edges[1:ntime+1])

lat_edges = np.linspace(-90, 90, nlat+1, endpoint=True, dtype=float)
lat_grid = 0.5 * (lat_edges[0:nlat] + lat_edges[1:nlat+1])
lat_min, lat_max = lat_edges[0:nlat], lat_edges[1:nlat+1]

lon_edges = np.linspace(lon0, lon0 + 360, nlon+1, endpoint=True, dtype=float)
lon_grid = 0.5 * (lon_edges[0:nlon] + lon_edges[1:nlon+1])

# instantiate count and data dictionaries
count_grid_dict = dict()
data_grid_dict = dict()

# initialize count and data arrays
count_grid = np.zeros((ntime, nlat, nlon), dtype=np.int32)
data_grid = np.zeros((ntime, nlat, nlon), dtype=np.float32)

files = sorted(glob(an.control_dict['obs']['test_obs']['files']))
obs_var = an.control_dict['test_setup']['obs_var']

# read obs
for filename in files:
    print('reading ' + filename)
    obs_ds = xr.open_dataset(filename)
    print(obs_ds.info())

    grid_util.update_sparse_data_grid(time_edges, lat_edges, lon_edges,
        obs_ds['timestamps'], obs_ds['lat'], obs_ds['lon'], obs_ds[obs_var],
        count_grid_dict, data_grid_dict)

    grid_util.update_data_grid(time_edges, lat_edges, lon_edges,
        obs_ds['timestamps'], obs_ds['lat'], obs_ds['lon'], obs_ds[obs_var],
        count_grid, data_grid)

grid_util.normalize_sparse_data_grid(count_grid_dict, data_grid_dict)
count_np, data_np = grid_util.to_np_array(time_edges, lat_edges, lon_edges,
    count_grid_dict, data_grid_dict)

grid_util.normalize_data_grid(count_grid, data_grid)

count_diff = count_grid - count_np
data_diff = data_grid[count_grid > 0] - data_np[count_grid > 0]
print('count diff min, max = %d, %d' % (count_diff.min(), count_diff.max()))
print('data diff min, max = %f, %f' % (data_diff.min(), data_diff.max()))

time_da = xr.DataArray(time_grid,
    attrs={'longname': 'time', 'units': 'seconds since 1970 Jan 01 00:00:00'})
lat_da = xr.DataArray(lat_grid,
    attrs={'longname': 'latitude', 'units': 'degree North'})
lon_da = xr.DataArray(lon_grid,
    attrs={'longname': 'longitude', 'units': 'degree East'})

count_da = xr.DataArray(count_np, dims=['time', 'lat', 'lon'],
    coords=[time_da, lat_da, lon_da])
data_da = xr.DataArray(data_np, dims=['time', 'lat', 'lon'],
    coords=[time_da, lat_da, lon_da])
"""
count_da = xr.DataArray(count_grid, dims=['time', 'lat', 'lon'],
    coords=[time_da, lat_da, lon_da])
data_da = xr.DataArray(data_grid, dims=['time', 'lat', 'lon'],
    coords=[time_da, lat_da, lon_da])
"""
ds = xr.Dataset({'count': count_da, 'data': data_da})
ds.to_netcdf('sparse_grid.nc')
