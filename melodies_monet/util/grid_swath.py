"""
file:  grid_swath.py
"""

import math
import numpy as np

def update_grid_data(lat_edges, lon_edges, lat_swath, lon_swath,
                     data_swath, count_grid, data_grid):
    """
    Accumulate swath data on a uniform grid

    Parameters
        lat_edges (np.array): grid latitude edges
        lon_edges (np.array): grid longitude edges
        lat_swath (np.array): swath latitudes
        lon_swath (np.array): swath longitudes
        data_swath (np.array): swath data values
        count_grid (np.array): number of swath points in grid cell
        data_grid (np.array): sum of data values in grid cell

    Returns
        None
    """
    logging.debug((len(lat_swath), len(lon_swath), len(data_swath)))
    lat_del = lat_edges[1] - lat_edges[0]
    lon_del = lon_edges[1] - lon_edges[0]
    nlat, nlon = data_grid.shape
    for i in range(len(data_swath)):
        if data_swath[i] >= 0:
            i_lat = math.floor((lat_swath[i] - lat_edges[0]) / lat_del)
            i_lon = math.floor((lon_swath[i] - lon_edges[0]) / lon_del)
            i_lat = min(nlat - 1, i_lat)
            i_lon = min(nlon - 1, i_lon)
            count_grid[i_lat, i_lon] += 1
            data_grid[i_lat, i_lon] += data_swath[i]

