#!/usr/bin/env python
# coding: utf-8

# ## MELODIES-MONET-Analysis Airnow prep 


import pandas as pd
import xarray as xr
import monetio as mio
from melodies_monet.util import write_util

# Here is an example of how you would use/import a variable from the master script
# ... python reformat_airnow.py $START_TIME $END_TIME ...
import sys
start_time_reformat=sys.argv[1]
end_time_reformat=sys.argv[2] 
print(sys.argv[1])
print(sys.argv[2])
dates = pd.date_range(start=start_time_reformat,end=end_time_reformat,freq='H')
#dates = pd.date_range(start='2021-05-28',end='2021-05-29',freq='H')


# helper function for local time.  Could be important for EPA statistics\n"
def get_local_time(ds):
    from numpy import zeros
    if 'utcoffset' in ds.data_vars:
        tim = t.time.copy()
        o = tim.expand_dims({'x':t.x.values}).transpose('time','x')
        on = xr.Dataset({'time_local':o,'utcoffset':t.utcoffset})
        y = on.to_dataframe()
        y['time_local'] = y.time_local + pd.to_timedelta(y.utcoffset, unit='H')
        time_local = y[['time_local']].to_xarray()
        ds = xr.merge([ds,time_local])
    return ds

df = mio.airnow.add_data(dates,wide_fmt=False,n_procs=12)

df = df.dropna(subset=['latitude','longitude']) # drop all values without an assigned latitude and longitude
dfp = df.rename({'siteid':'x'},axis=1).pivot_table(values='obs',index=['time','x'], columns=['variable']) # convert to wide format
dfx = dfp.to_xarray() # convert to xarray
# df.head()


# When converting to wide format we have to remerge the site data back into the file.
dfpsite = df.rename({'siteid':'x'},axis=1).drop_duplicates(subset=['x']) # droping duplicates and renaming
# convert sites to xarray
test = dfpsite.drop(['time','time_local','variable','obs'],axis=1).set_index('x').dropna(subset=['latitude','longitude']).to_xarray()
# merge sites back into the data
t = xr.merge([dfx,test])
# get local time
tt = get_local_time(t)


# add siteid back as a variable and create x as an array of integers
tt['siteid'] = (('x'),tt.x.values)
tt['x'] = range(len(tt.x))
# expand dimensions so that it is (time,y,x)
t = tt.expand_dims('y').set_coords(['siteid','latitude','longitude']).transpose('time','y','x')


#wite out to filename set in cell 2
write_util.write_ncf(t,'test5.nc')
