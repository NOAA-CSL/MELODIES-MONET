import sys
import argparse
import logging
import yaml

import numpy as np
import pandas as pd
import xarray as xr

parser = argparse.ArgumentParser()
parser.add_argument('--control', type=str,
    default='control.yaml',
    help='yaml control file')
parser.add_argument('--logfile', type=str,
    default=sys.stdout,
    help='log file (default stdout)')
parser.add_argument('--debug', action='store_true',
    help='set logging level to debug')
parser.add_argument('--nfile', type=int,
    default=1,
    help='number of obs files')
args = parser.parse_args()

"""
Setup logging
"""
logging_level = logging.DEBUG if args.debug else logging.INFO
logging.basicConfig(stream=args.logfile, level=logging_level)

"""
Read YAML control
"""
with open(args.control, 'r') as f:
    control = yaml.safe_load(f)

start_time = pd.to_datetime(control['analysis']['start_time'])
end_time = pd.to_datetime(control['analysis']['end_time'])
datetime_indices = pd.date_range(start_time, end_time,
    freq=control['test_setup']['freq'])
ntime = len(datetime_indices)

ntime_file = ntime // args.nfile
logging.debug('ntime_file:%d' % ntime_file)

"""
Generate random test observations
"""
np.random.seed(control['test_setup']['random_seed'])

var_names = control['obs']['test_obs']['variables'].keys()

for ifile in range(args.nfile):

    df_dict = dict()

    for var_name in var_names:
 
        if 'range_min' in control['obs']['test_obs']['variables'][var_name]:
            range_min = control['obs']['test_obs']['variables'][var_name]['range_min']
        else:
            range_min = 0

        if 'range_max' in control['obs']['test_obs']['variables'][var_name]:
            range_max = control['obs']['test_obs']['variables'][var_name]['range_max']
        else:
            range_max = 1

        df_dict[var_name] = (range_max - range_min) * np.random.rand(ntime_file) + range_min

    istart = ifile * ntime_file
    iend = istart + ntime_file
    df = pd.DataFrame(df_dict, index=datetime_indices[istart:iend]).to_xarray()
    ds = xr.Dataset(df)
    logging.debug(ds)
    suffix = '_%d.nc' % ifile
    ds.to_netcdf(control['obs']['test_obs']['filename'].replace('.nc', suffix))

