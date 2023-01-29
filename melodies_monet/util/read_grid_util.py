import logging
import xarray as xr
from monetio.sat._gridded_eos_mm import read_gridded_eos

from analysis_util import fill_date_template, find_file


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
