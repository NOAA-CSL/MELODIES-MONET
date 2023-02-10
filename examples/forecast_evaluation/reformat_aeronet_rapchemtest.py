import monetio as mio
import numpy as np
import pandas as pd
import xarray as xr
from melodies_monet.util import write_util
import os

# set the dates
### for some reason the 7th of august 2019 has an error reading in the parallel version.
# Not sure why but some of the dtypes are different when processing in parallel vs serial.
# Here we will use the serial version just to make it easy for everyone to understand.
### For faster processing in other periods add the `n_procs=n` kwarg where n is the number of processes
import sys
start_time_reformat=sys.argv[1]
end_time_reformat=sys.argv[2]
print(sys.argv[1])
print(sys.argv[2])
dates = pd.date_range(start=start_time_reformat,end=end_time_reformat,freq='H')

#dates = pd.date_range(start='2021-08-01',end='2021-09-01',freq='H')
# set the output filename
#outname = 'AERONET_L15_20190901_20190930.nc'
outname = 'test5.nc'
# set standard wavelengths
standard_wavelengths = np.array([0.34, 0.44, 0.55, 0.66, 0.86, 1.63, 11.1])* 1000. # convert from micron to nm
# get the data
df = mio.aeronet.add_data(dates, interp_to_aod_values=standard_wavelengths, freq='H') # ,n_procs=12)


# dfp = df.rename({'siteid':'x'},axis=1).set_index(['time','x']).drop_duplicates()
# df = df3
verbose=False
dfp = df.rename({'siteid':'x'},axis=1).set_index(['time','x'])
columns = dfp.columns.to_list()
columns2 = []
remove_columns = []
for i in np.arange(len(columns)):
    columns2.append(columns[i])
    try:
        dfp[columns2].to_xarray()
        if verbose:
            print('COLUMN SUCCESS:',columns[i])
    except:
        if verbose:
            print('COLUMN FAILURE:',columns[i])
        remove_columns.append(columns[i])
        columns2.remove(columns[i])
if verbose:
    print(columns2)
dft = df.drop(remove_columns,axis=1)
dfp = dfp.drop(remove_columns,axis=1).dropna(subset=['latitude','longitude'])
dfx = dfp.to_xarray()

dsets = []
for s in df.siteid.unique():
    dsets.append(dft.loc[df.siteid == s].set_index(['time']).to_xarray())


# now site_variable are the single element attributes of the individual site so lets simplify this

from numpy import unique
site_variable = ['siteid','latitude','longitude','aeronet_instrument_number','elevation']

def expand_dims(ds, index=0):
    # first set a new index for the siteid
    ds['x'] = index
    ds = ds.expand_dims(['x'])
    ds = ds.set_coords(['x'])
    # now reduce the site variables to single element variables
    for sv in site_variable:
        tmp = [unique(ds[sv])[0]]
        ds[sv] = (('x',), tmp)
    return ds

# now site_variable are the single element attributes of the individual site so lets simplify this
for index,d in enumerate(dsets):
    dsets[index] = expand_dims(d,index=index)

    # now combine all the datasets for each site into a single dataset
ds = xr.concat(dsets,dim='x').set_coords(site_variable)
# ds
os.path.join

# write the file
t = ds.expand_dims('y').transpose('time','y','x')
write_util.write_ncf(t,outname)
