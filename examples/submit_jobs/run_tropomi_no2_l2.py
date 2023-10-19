import os
import sys
sys.path.append('../../')
from melodies_monet import driver

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
import numpy as np

plt.set_loglevel (level = 'warning')

an = driver.analysis()
an.control = '../yaml/control_tropomi_l2_no2.yaml'
an.read_control()

# --- satobs ---
an.open_obs()
an.control_dict

# --- model ---
an.open_models()
an.models
#an.models['wrfchem_v4.2'].obj

# --- paring ---
an.pair_data()
paired_obs = an.paired['tropomi_l2_no2_wrfchem_v4.2'].obj
model = an.models['wrfchem_v4.2'].obj

# --- plotting ---
#lat     = paired_obs['latitude']
#lon     = paired_obs['longitude']
#no2grid = paired_obs['no2'][0,:,:]


paired_obs_stack = paired_obs.set_index(x=("ll", "y")).unstack("x")

lat = an.models['wrfchem_v4.2'].obj.coords['latitude']
lon = an.models['wrfchem_v4.2'].obj.coords['longitude']

print('!!! paired obs stack', paired_obs_stack)

no2grid = paired_obs_stack['nitrogendioxide_tropospheric_column']
no2grid = no2grid[0,:,:]
ind_x = paired_obs.coords['x']
ind_y = paired_obs.coords['ll']
print(no2grid, np.nanmin(no2grid), np.nanmax(no2grid))

fig = plt.figure(figsize=(20,10))
ax = plt.axes(projection=ccrs.PlateCarree())
clev = np.arange(0, 1e16, 1*1e15)
plt.contourf(lon, lat, no2grid, clev, cmap='Spectral_r', extend='both')
cbar=plt.colorbar(shrink=0.6)
plt.show()
fig.savefig('/Users/mengli/Work/melodies-monet/outdata/paried_trp_no2_20190715.png')
#print('paired_obs', paired_obs)

#print(paired_obs)
#lat = an.models['wrfchem_v4.2'].obj.coords['latitude']
#lon = an.models['wrfchem_v4.2'].obj.coords['longitude']
print('an models', an.models['wrfchem_v4.2'].obj)
no2grid = an.models['wrfchem_v4.2'].obj['no2trpcol'][1,0,:,:]
print(no2grid, np.nanmin(no2grid), np.nanmax(no2grid))

fig = plt.figure(figsize=(20,10))
ax = plt.axes(projection=ccrs.PlateCarree())
clev = np.arange(0, 1e16, 1*1e15)
plt.contourf(lon, lat, no2grid, clev, cmap='Spectral_r', extend='both')
cbar=plt.colorbar(shrink=0.6)
plt.show()
fig.savefig('/Users/mengli/Work/melodies-monet/outdata/paried_wrfchem_no2_20190715.png')
print('paired_obs', paired_obs)

# --- save paired data ---
an.save_analysis()

# --- read saved paired data ---
#an.read_analysis()
#an.paired['tropomi_l2_no2_wrfchem_v4.2'].obj
