import numpy as np
import pandas as pd
import xarray as xr

from glob import glob
from melodies_monet import driver

from melodies_monet.util import grid_util

an = driver.analysis()

an.control = 'test_grid.yaml'
an.read_control()

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
count_grid = dict()
data_grid = dict()

files = glob(an.control_dict['obs']['test_obs']['files'])
obs_var = an.control_dict['test_setup']['obs_var']

# read obs
for filename in files:
    print('reading ' + filename)
    obs_ds = xr.open_dataset(filename)
    print(obs_ds.info())

    grid_util.update_sparse_data_grid(time_edges, lat_edges, lon_edges,
        obs_ds['timestamps'], obs_ds['lat'], obs_ds['lon'], obs_ds[obs_var],
        count_grid, data_grid)

# print(count_grid)
# print(data_grid)
