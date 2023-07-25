import os
import logging
import xarray as xr
from monetio.sat._gridded_eos_mm import read_gridded_eos

from .analysis_util import fill_date_template, find_file


def read_grid_models(config, date_str):
    """
    Read grid data models

    Parameters
        config (dict): configuration dictionary
        date_str (str yyyy-mm-m_abbr-dd-ddd): date string

    Returns
        model_datasets (dict of xr.Dataset): dictionary of model datasets
    """
    model_datasets = dict()

    for model_name in config['model']:

        datadir = config['model'][model_name]['datadir']
        filestr = fill_date_template(
            config['model'][model_name]['files'], date_str)
        filename = find_file(datadir, filestr)

        model_datasets[model_name] = xr.open_dataset(filename)

    return filename, model_datasets


def read_grid_obs(config, obs_vars, date_str):
    """
    Read grid data obs

    Parameters
        config (dict): configuration dictionary
        obs_vars (dict of dict):
            nested dictionary keyed by obs set name and obs variable name
        date_str (str yyyy-mm-m_abbr-dd-ddd): date string

    Returns
        obs_datasets (dict of xr.Dataset): dictionary of obs datasets
    """
    obs_datasets = dict()

    yyyy_str, mm_str, m_abbr_str, dd_str, ddd_str \
        = tuple(date_str.split('-'))

    for obs_name in obs_vars:

        data_format = config['obs'][obs_name]['data_format']
        datadir = config['obs'][obs_name]['datadir']
        filestr = fill_date_template(
            config['obs'][obs_name]['filename'], date_str)
        filename = find_file(datadir, filestr)

        file_extension = os.path.splitext(filename)[1]

        if data_format == 'gridded_eos':
            if file_extension == '.hdf':
                ds_obs = read_gridded_eos(
                    filename, obs_vars[obs_name])
                filename_nc = filename.replace('.hdf', '.nc')
                logging.info('writing ' + filename_nc)
                ds_obs.to_netcdf(filename_nc)
            else:
                ds_obs = xr.open_dataset(filename)
        else:
            ds_obs = xr.open_dataset(filename)

        obs_datasets[obs_name] = ds_obs

    return filename, obs_datasets

