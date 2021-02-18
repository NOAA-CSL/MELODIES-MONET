#!/usr/bin/env python                                                                                                                                                                                                               

__author__ = 'Patrick Campbell'
__email__ = 'patrick.c.campbell@noaa.gov'
__license__ = 'GPL'

#Simple MONET utility to command line pair model vs. observations                                                                                                                                                                   

import os
from glob import glob
import sys

import subprocess
from distutils.spawn import find_executable
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import monet
from monet.util.tools import long_to_wide
import pandas as pd

if __name__ == '__main__':
    parser = ArgumentParser(
        description='combines paired data files for multiple model runs',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-p',
        '--paired_files',
        help='string input paired data files (>=2)',
        nargs='+',
        type=str,
        required=True)
    parser.add_argument(
        '-s',
        '--model_species',
        help='string/list input for obs species-variables to pair',
        type=str,
        nargs='+',
        required=False,
        default=['OZONE'])
    parser.add_argument(
        '-mdf',
        '--mergedf',
        help='boolean operator to merge entire dataframes (slow)',
        type=bool,
        required=False,
        default=False)
    parser.add_argument(
        '-o',
        '--output',
        help='string output path/filename for combined paired dataframe',
        type=str,
        required=False,
        default='./AIRNOW_CMAQ_merged_pair')
    parser.add_argument(
        '-avg',
        '--averagef',
        help='boolean operator to average entire dataframes (slow)',
        type=bool,
        required=False,
        default=False)
    parser.add_argument(
        '-v',
        '--verbose',
        help='print debugging information',
        action='store_true',
        required=False,
        default=False)
    args = parser.parse_args()

    paired_files = args.paired_files
    species = args.model_species
    output = args.output
    mergedf = args.mergedf
    avg = args.averagef
    
    mapping_table = {'OZONE':'O3', 'PM2.5':'PM25_TOT', 'PM10':'PM10_new', 'CO':'CO_new', 'NO':'NO_new', 'NO2':'NO2_new', 'SO2':'SO2_new','NOX':'NOX_new','NOY':'NOY_new','TEMP':'TEMP2','WS':'WSPD10','WD':'WDIR10','SRAD':'GSW','BARPR':'PRSFC','PRECIP':'RT','RHUM':'Q2'}
    sub_map  = {i: mapping_table[i] for i in species if i in mapping_table}
    sub_maps = pd.Series(sub_map)
    print('species...')
    print(sub_maps)
    print('combining paired files...')
    print(paired_files)

    count=0
    for i in paired_files:
        df=pd.read_hdf(i) 
        if count == 0:
                df_merge=df
        else: 
                if mergedf is False:
                	df_species = df[sub_maps]
                	df_species = df_species.add_suffix('_'+str(count+1))
                	if avg is True:
                        	print('averaging dataframes, both obs and mod columns (e.g., for multiple years)...')
                        	df_merge['date']=df_merge.time.dt.strftime('%m-%d-%H-%M')
                       		df['date']=df.time.dt.strftime('%m-%d-%H-%M')
                        	df_merge = pd.concat([df_merge, df]).groupby(['date','siteid','epa_region','state_name'], as_index=False).mean()
                        	df_merge['time']=df.time
                        	df_merge['time_local']=df.time_local
                	else:
               			print('joining dataframes, just add mod columns (e.g., for multiple model runs)...')
                		df_merge=df_merge.join(df_species)
                else:
                	df[species] = df[sub_maps].add_suffix('_'+str(count+1))
                	print('merging entire dataframes...slow...')
                	df_merge=df_merge.merge(df,on='time')
        count=count+1
    print('final merged data frame...')
    print(df_merge.keys())
    print(df_merge)

    df_merge.to_hdf(output+'.hdf',
            'df_merge',
            format='table',
            mode='w')
    df_merge.to_csv(output+'.csv')
    sys.exit(0)
