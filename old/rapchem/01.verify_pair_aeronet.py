#!/usr/bin/env python

__author__ = 'Patrick Campbell'
__email__ = 'patrick.c.campbell@noaa.gov'
__license__ = 'GPL'

#Simple MONET utility to command line pair model vs. observations

import os
import subprocess
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from distutils.spawn import find_executable
from glob import glob

import pandas as pd

import monet
from monet.util.tools import long_to_wide


def pair_point(da, df, sub_map, interp):
    dfpair = da[sub_map].monet.combine_point(
        df,col='aod_550nm',method=interp, reuse_weights=True)
    return dfpair


def get_aeronet(start, end):
    dates = pd.date_range(start=start, end=end, freq='H')
    df = monet.obs.aeronet.add_data(dates, freq='H')
    return df.dropna(subset=['latitude', 'longitude'])


def open_cmaq(finput, verbose=False):
    dset = monet.models.fv3chem.open_mfdataset(finput)
    return dset


if __name__ == '__main__':

    parser = ArgumentParser(
        description='pairs cmaq model data to aqs observations',
        formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '-f',
        '--files',
        help='string input model file directory/names',
        type=str,
        nargs='+',
        required=True)
    #    parser.add_argument('-s', '--startdates',  help='string input start date for pairing YYYY-MM-DD HH:MM:SS', type=str, required=True)
    #    parser.add_argument('-e', '--enddates',    help='string input end date for pairing YYYY-MM-DD HH:MM:SS', type=str, required=True)
    parser.add_argument(
        '-x',
        '--species',
        help='string input for obs species-variables to pair',
        type=str,
        nargs='+',
        required=False,
        default=['aod_550nm'])
    parser.add_argument(
        '-o',
        '--output',
        help='string output path for paired dataframe, stats, plots',
        type=str,
        required=False,
        default='paired_output.hdf')
    parser.add_argument(
        '-p',
        '--path',
        help='string path to director of network observations',
        type=str,
        required=False,
        default='/data/aqf2/barryb/5xpm/AQS_DATA/')
    parser.add_argument(
        '-n',
        '--network',
        help='string input data network name: airnow, aqs',
        type=str,
        required=False,
        default='airnow')
    parser.add_argument(
        '-l',
        '--label',
        help='Add this string to all data variable names',
        type=str,
        required=False,
        default='')
    parser.add_argument(
        '-m',
        '--model',
        help='input model: cmaq, fv3, hysplit (not-ready), or camx (not-ready)',
        type=str,
        required=False,
        default='cmaq')
    parser.add_argument(
        '-i',
        '--interp',
        help=
        'xesmf interpolation scheme, bilinear, conservative, nearest_s2d, nearest_d2s, patch',
        type=str,
        required=False,
        default='bilinear')
    parser.add_argument(
        '-v',
        '--verbose',
        help='print debugging information',
        action='store_true',
        required=False)
    args = parser.parse_args()

    finput = args.files
    #    start   = args.startdates
    #    end     = args.enddates
    species = args.species
    output = args.output
    datapath = args.path
    network = args.network
    model = args.model
    interp = args.interp
    verbose = args.verbose
    label= args.label
    #reads model output (cmaq default)

    if model == 'cmaq':
        da = open_cmaq(finput, verbose=verbose)
    else:
        print('Must enter cmaq model right now')
        raise RuntimeError

#retrieves data observations and formats pandas dataframe (airnow default)
    start = da.time.to_index().min() #- pd.Timedelta(1, unit='h')
    end = da.time.to_index().max() #+ pd.Timedelta(1, unit='h')
    # print(start.strftime('%Y%m%d%H'), end.strftime('%Y%m%d%H'))
    # print(start)
    # print(end)
    if os.path.isfile(output):
        df = pd.read_hdf(output)
    else:
        df = get_aeronet(start, end)
#    df = df.rename({'latitude':'lat','longitude':'lon'},axis=1) #pairs surface point-type observations with 2D model parameters

    mapping_table = {'aod_550nm': 'pm25aod550',
#        'dust25aod550': 'aod_550nm',
#        'salt25aod550': 'aod_550nm',
#        'sulf25aod550': 'aod_550nm',
#        'oc25aod550': 'aod_550nm',
#        'bc25aod550': 'aod_550nm'
    } 
    sub_map = {i: mapping_table[i] for i in species if i in mapping_table}
    use_these = [sub_map[i] for i in sub_map.keys()]
    dfpair = pair_point(da, df, use_these, interp)
    # dfpair = dfpair.dropna(subset=['aod_550nm', 'pm25aod550'])
    if label != '':
        dfpair = dfpair.rename({'pm25aod550':'pm25aod550_{}'.format(label),
#                            'dust25aod550':'dust25aod550_{}'.format(label),
#                            'salt25aod550':'salt25aod550_{}'.format(label),
#                            'sulf25aod550':'sulf25aod550_{}'.format(label),
#                            'oc25aod550':'oc25aod550_{}'.format(label),
#                            'bc25aod550':'bc25aod550_{}'.format(label)}, axis=1)
                             }, axis=1)
#    dfpair.to_hdf(output, 'df', format='table')

    dfpair.to_hdf(
            'AERONET_FV3CHEM_' + start.strftime('%Y-%m-%d-%H') + '_' +
            end.strftime('%Y-%m-%d-%H') + '_pair.hdf',
            'dfpair',
            format='table',
            mode='w')


    sys.exit(0)
