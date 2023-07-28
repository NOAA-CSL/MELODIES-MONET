import os
import logging
import xarray as xr
from monetio.sat._gridded_eos_mm import read_gridded_eos

from .analysis_util import fill_date_template, find_file


def read_grid_models(config, date_str, model=None):
    """
    Read grid data models

    Parameters
        config (dict): configuration dictionary
        date_str (str yyyy-mm-m_abbr-dd-ddd): date string
        model: specific model to read
            optional, if not specified all model in config['models'] will be read

    Returns
        model_datasets (dict of xr.Dataset): dictionary of model datasets
        filenames (dict of str): dictionary of filenames
    """
    model_datasets = dict()
    filenames = dict()

    if model is not None:
        model_list = [model]
    else:
        model_list = config['model']

    for model_name in model_list:

        datadir = config['model'][model_name]['datadir']
        filestr = fill_date_template(
            config['model'][model_name]['files'], date_str)
        filename = find_file(datadir, filestr)

        model_datasets[model_name] = xr.open_dataset(filename)
        filenames[model_name] = filename

    return model_datasets, filenames


def read_grid_obs(config, obs_vars, date_str, obs=None):
    """
    Read grid data obs

    Parameters
        config (dict): configuration dictionary
        obs_vars (dict of dict):
            nested dictionary keyed by obs set name and obs variable name
        date_str (str yyyy-mm-m_abbr-dd-ddd): date string
        obs: specific observation to read
            optional, if not specified all obs in obs_vars will be read

    Returns
        obs_datasets (dict of xr.Dataset): dictionary of obs datasets
        filenames (dict of str): dictionary of filenames
    """
    obs_datasets = dict()
    filenames = dict()

    if obs is not None:
        obs_list = [obs]
    else:
        obs_list = obs_vars.keys()

    yyyy_str, mm_str, m_abbr_str, dd_str, ddd_str \
        = tuple(date_str.split('-'))

    for obs_name in obs_list:

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
        filenames[obs_name] = filename

    return obs_datasets, filenames

