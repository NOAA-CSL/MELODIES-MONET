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
    Store running counts and sums in dictionaries keyed by grid index tuples (i_time, i_lat, i_lon)

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
            i_time = np.clip(i_time, 0, ntime - 1)
            i_lat = np.clip(i_lat, 0, nlat - 1)
            i_lon = np.clip(i_lon, 0, nlon - 1)
            if (i_time, i_lat, i_lon) in count_grid.keys():
                count_grid[(i_time, i_lat, i_lon)] += 1
                data_grid[(i_time, i_lat, i_lon)] += data_obs[i].values
            else:
                count_grid[(i_time, i_lat, i_lon)] = 1
                data_grid[(i_time, i_lat, i_lon)] = data_obs[i].values


def normalize_sparse_data_grid(count_grid, data_grid):
    """
    Normalize accumulated data on a uniform grid

    Parameters
        count_grid (dict): number of obs points in grid cell
        data_grid (dict): sum of data values in grid cell

    Returns
        None
    """
    for index_tuple in count_grid.keys():
        data_grid[index_tuple] /= count_grid[index_tuple]


def sparse_data_to_array(time_edges, lat_edges, lon_edges,
                         count_grid, data_grid):
    """
    Convert sparse grid data to numpy arrays

    Parameters
        time_edges (np.array): grid time edges
        lat_edges (np.array): grid latitude edges
        lon_edges (np.array): grid longitude edges
        count_grid (dict): number of obs points in grid cell
        data_grid (dict): sum of data values in grid cell

    Returns
        count_grid_array (np.array): number of obs points in grid cell
        data_grid_array (np.array): sum of data values in grid cell
    """
    ntime, nlat, nlon = len(time_edges) - 1, len(lat_edges) - 1, len(lon_edges) - 1
    count_grid_array = np.zeros((ntime, nlat, nlon), dtype=np.int32)
    data_grid_array = np.full((ntime, nlat, nlon), np.nan, dtype=np.float32)
    for index_tuple in count_grid.keys():
        count_grid_array[index_tuple[0], index_tuple[1], index_tuple[2]] = count_grid[index_tuple]
        data_grid_array[index_tuple[0], index_tuple[1], index_tuple[2]] = data_grid[index_tuple]

    return count_grid_array, data_grid_array


def update_data_grid(time_edges, lat_edges, lon_edges,
                     time_obs, lat_obs, lon_obs, data_obs,
                     count_grid, data_grid):
    """
    Accumulate obs data on a uniform grid with dimensions time x lat x lon
    Store running counts and sums in numpy arrays

    Parameters
        time_edges (np.array): grid time edges
        lat_edges (np.array): grid latitude edges
        lon_edges (np.array): grid longitude edges
        time_obs (np.array): obs times
        lat_obs (np.array): obs latitudes
        lon_obs (np.array): obs longitudes
        data_obs (np.array): obs data values
        count_grid (np.array): number of obs points in grid cell
        data_grid (np.array): sum of data values in grid cell

    Returns
        None
    """
    time_del = time_edges[1] - time_edges[0]
    lat_del = lat_edges[1] - lat_edges[0]
    lon_del = lon_edges[1] - lon_edges[0]
    ntime, nlat, nlon = data_grid.shape
    for i in range(len(data_obs)):
        if not np.isnan(data_obs[i]):
            i_time = math.floor((time_obs[i] - time_edges[0]) / time_del)
            i_lat = math.floor((lat_obs[i] - lat_edges[0]) / lat_del)
            i_lon = math.floor((lon_obs[i] - lon_edges[0]) / lon_del)
            i_time = np.clip(i_time, 0, ntime - 1)
            i_lat = np.clip(i_lat, 0, nlat - 1)
            i_lon = np.clip(i_lon, 0, nlon - 1)
            count_grid[i_time, i_lat, i_lon] += 1
            data_grid[i_time, i_lat, i_lon] += data_obs[i]


def normalize_data_grid(count_grid, data_grid):
    """
    Normalize accumulated data on a uniform grid

    Parameters
        count_grid (np.array): number of obs points in grid cell
        data_grid (np.array): sum of data values in grid cell

    Returns
        None
    """
    data_grid[count_grid == 0] = np.nan
    data_grid[count_grid > 0] /= count_grid[count_grid > 0]

