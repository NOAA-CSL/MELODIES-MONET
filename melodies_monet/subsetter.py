import os
import sys
import argparse
import logging
import yaml
from glob import glob

parser = argparse.ArgumentParser()
parser.add_argument('--control', type=str,
    default='control.yaml',
    help='yaml control file')
parser.add_argument('--logfile', type=str,
    default=sys.stdout,
    help='log file (default stdout)')
parser.add_argument('--debug', action='store_true',
    help='set logging level to debug')
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

logging.debug(control)

for model in control['model']:
    logging.info('processing:' + model)

    var_str = '-v '
    for dataset in control['model'][model]['mapping']:
        for var in control['model'][model]['mapping'][dataset]:
            var_str += var + ','
    var_str = var_str.strip(',')
    logging.info(var_str)

    files = sorted(glob(control['model'][model]['files']))
    for file_in in files:
        file_out = file_in + '_subset'
        command = 'ncks -O ' + var_str + ' ' + file_in + ' ' + file_out
        logging.info(command)
        os.system(command)

