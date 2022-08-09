"""
file: grid_util.py
"""

import math
import numpy as np

def update_data_grid(lat_edges, lon_edges, lat_obs, lon_obs,
                     data_obs, count_grid, data_grid):
    """
    Accumulate obs data on a uniform grid

    Parameters
        lat_edges (np.array): grid latitude edges
        lon_edges (np.array): grid longitude edges
        lat_obs (np.array): obs latitudes
        lon_obs (np.array): obs longitudes
        data_obs (np.array): obs data values
        count_grid (np.array): number of obs points in grid cell
        data_grid (np.array): sum of data values in grid cell

    Returns
        None
    """
    logging.debug((len(lat_obs), len(lon_obs), len(data_obs)))
    lat_del = lat_edges[1] - lat_edges[0]
    lon_del = lon_edges[1] - lon_edges[0]
    nlat, nlon = data_grid.shape
    for i in range(len(data_obs)):
        if not np.isnan(data_obs[i]):
            i_lat = math.floor((lat_obs[i] - lat_edges[0]) / lat_del)
            i_lon = math.floor((lon_obs[i] - lon_edges[0]) / lon_del)
            i_lat = min(nlat - 1, i_lat)
            i_lon = min(nlon - 1, i_lon)
            count_grid[i_lat, i_lon] += 1
            data_grid[i_lat, i_lon] += data_obs[i]


def update_sparse_data_grid(lat_edges, lon_edges, lat_obs, lon_obs,
                            data_obs, count_grid, data_grid):
    """
    Accumulate obs data on a uniform grid

    Parameters
        lat_edges (np.array): grid latitude edges
        lon_edges (np.array): grid longitude edges
        lat_obs (np.array): obs latitudes
        lon_obs (np.array): obs longitudes
        data_obs (np.array): obs data values
        count_grid (dict): number of obs points in grid cell
        data_grid (dict): sum of data values in grid cell

    Returns
        None
    """
    logging.debug((len(lat_obs), len(lon_obs), len(data_obs)))
    lat_del = lat_edges[1] - lat_edges[0]
    lon_del = lon_edges[1] - lon_edges[0]
    nlat, nlon = data_grid.shape
    for i in range(len(data_obs)):
        if not np.isnan(data_obs[i]):
            i_lat = math.floor((lat_obs[i] - lat_edges[0]) / lat_del)
            i_lon = math.floor((lon_obs[i] - lon_edges[0]) / lon_del)
            i_lat = min(nlat - 1, i_lat)
            i_lon = min(nlon - 1, i_lon)
            if (i_lon, ilat) in count_grid.keys():
                count_grid[(i_lon, i_lat)] += 1
                data_grid[(i_lat, i_lon)] += data_obs[i]
            else:
                count_grid[(i_lon, i_lat)] = 1
                data_grid[(i_lat, i_lon)] = data_obs[i]

