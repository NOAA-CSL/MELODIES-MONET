# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#
def subset_model_filelist(all_files,timeformat,timestep,timeinterval):
    '''Subset model filelist to within a given time interval. 
    Filename requirements:
    - individual files for each timestep
    - time must be in filename
    '''
    import pandas as pd
    subset_interval = pd.date_range(start=timeinterval[0],end=timeinterval[-1],freq=timestep)
    interval_files = []
    for i in subset_interval:
        flst = [fs for fs in all_files if i.strftime(timeformat) in fs]
        if len(flst) == 1:
            interval_files.append(flst[0])
        elif len(flst) >1:
            print('More than 1 file for {} in listing'.format(i.strftime(timeformat)))
    return interval_files

def subset_OMPS_l2(file_path,timeinterval):
    '''Dependent on filenaming convention
    '''
    import pandas as pd
    from glob import glob
    import fnmatch
    all_files = glob(file_path)
    interval_files = []
    subset_interval = pd.date_range(start=timeinterval[0],end=timeinterval[-1],freq='D',inclusive='left')
    
    for i in subset_interval:
        fst = fnmatch.filter(all_files,'*OMPS-NPP_NMTO3-L2_v*_{}*_o*'.format(i.strftime('%Ym%m%d')))
        fst.sort()
        for j in fst:
            interval_files.append(j)
    return interval_files

def subset_MODIS_l2(file_path,timeinterval):
    '''Dependent on filenaming convention
       MOD04_L2.AYYYYDDD.HHMM.collection.timestamp.hdf
    '''
    import pandas as pd
    from glob import glob
    import fnmatch

