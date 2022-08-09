import numpy as np
import xarray as xr

from glob import glob
from melodies_monet import driver

from melodies_monet.util import grid_util

an = driver.analysis()

an.control = 'test_grid.yaml'
an.read_control()

# generate uniform grid
nlon = an.control_dict['obs']['test_obs']['sparse_data_grid']['nlon']
nlat = an.control_dict['obs']['test_obs']['sparse_data_grid']['nlat']
lon0 = an.control_dict['obs']['test_obs']['sparse_data_grid']['lon0']
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

    grid_util.update_sparse_data_grid(lat_edges, lon_edges,
        obs_ds['lat'], obs_ds['lon'], obs_ds[obs_var],
        count_grid, data_grid)

print(count_grid)
print(data_grid)
