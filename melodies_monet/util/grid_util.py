# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#

"""
file: grid_util.py
"""

import math
import numpy as np


def update_sparse_data_grid(time_edges, lat_edges, lon_edges,
                            time_obs, lat_obs, lon_obs, data_obs,
                            count_grid, data_grid):
    """
    Accumulate obs data on a uniform grid with dimensions time x lat x lon
    Store running counts and sums
        in dictionaries keyed by grid index tuples (i_time, i_lat, i_lon)

    Parameters
        time_edges (np.array): grid time edges
        lat_edges (np.array): grid latitude edges
        lon_edges (np.array): grid longitude edges
        time_obs (np.array): obs times
        lat_obs (np.array): obs latitudes
        lon_obs (np.array): obs longitudes
        data_obs (np.array): obs data values
        count_grid (dict): number of obs points in grid cell
        data_grid (dict): sum of data values in grid cell

    Returns
        None
    """
    time_del = time_edges[1] - time_edges[0]
    lat_del = lat_edges[1] - lat_edges[0]
    lon_del = lon_edges[1] - lon_edges[0]
    ntime, nlat, nlon = len(time_edges) - 1, len(lat_edges) - 1, len(lon_edges) - 1
    for i in range(len(data_obs)):
        if not np.isnan(data_obs[i]):
            i_time = math.floor((time_obs[i] - time_edges[0]) / time_del)
            i_lat = math.floor((lat_obs[i] - lat_edges[0]) / lat_del)
            i_lon = math.floor((lon_obs[i] - lon_edges[0]) / lon_del)
            i_time = min(ntime - 1, i_time)
            i_lat = min(nlat - 1, i_lat)
            i_lon = min(nlon - 1, i_lon)
            if (i_time, i_lon, i_lat) in count_grid.keys():
                count_grid[(i_time, i_lat, i_lon)] += 1
                data_grid[(i_time, i_lat, i_lon)] += data_obs[i].values
            else:
                count_grid[(i_time, i_lat, i_lon)] = 1
                data_grid[(i_time, i_lat, i_lon)] = data_obs[i].values


def update_data_grid_2d(lat_edges, lon_edges,
                        lat_obs, lon_obs, data_obs,
                        count_grid, data_grid):
    """
    Accumulate obs data on a uniform grid with dimensions lat x lon
    Store running counts and sums
        in 2d numpy arrays
    This method may become unused,
        it is retained for comparison to the sparse data grid

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

