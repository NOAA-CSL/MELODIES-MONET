import logging
import xarray as xr
from monetio.sat._gridded_eos_mm import read_gridded_eos

from analysis_util import fill_date_template, find_file


def get_obs_vars(config):
    """
    Get subset of obs variables from model to obs variable mapping

    Parameters
        config (dict): configuration dictionary

    Returns
        obs_vars_subset (dict of dicts):
            nested dictionary keyed by obs set name and obs variable name
    """
    obs_vars_subset = dict()

    for model_name in config['model']:

        mapping = config['model'][model_name]['mapping']

        for obs_name in mapping:
            obs_vars = config['obs'][obs_name]['variables']
            obs_vars_subset[obs_name] = dict()

            for model_var in mapping[obs_name]:
                obs_var = mapping[obs_name][model_var]
                obs_vars_subset[obs_name][obs_var] = obs_vars[obs_var]

    return obs_vars_subset


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
        filestr = config['model'][model_name]['filestr']
        filestr = fill_date_template(
            config['model'][model_name]['filestr'], date_str)
        filename = find_file(datadir, filestr)

        model_datasets[model_name] = xr.open_dataset(filename)

    return model_datasets
