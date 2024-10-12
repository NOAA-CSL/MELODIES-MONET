import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from melodies_monet import driver

import warnings
warnings.filterwarnings('ignore')

an = driver.analysis()
an.control = 'control_modis_l2.yaml'
an.read_control()

an.open_models()

an.setup_obs_grid()
# print(an.obs_grid)

an.setup_regridders()

for time_interval in an.time_intervals:

    print(time_interval)

    an.open_obs(time_interval=time_interval)
    an.update_obs_gridded_data()

an.normalize_obs_gridded_data()
# print(an.obs_gridded_dataset)

an.obs_gridded_dataset.to_netcdf('MODIS.nc')


varname_terra_aod = 'Terra_MODIS_AOD_550_Dark_Target_Deep_Blue_Combined_data'
varname_aqua_aod = 'Aqua_MODIS_AOD_550_Dark_Target_Deep_Blue_Combined_data'
varname_merra2_aod = 'TOTEXTTAU'

ax = plt.subplot(projection=ccrs.PlateCarree())
an.obs_gridded_dataset[varname_terra_aod].isel(time=19).plot.pcolormesh(
    cmap=plt.cm.turbo,
    cbar_kwargs={'location': 'bottom', 'label': 'AOD'},
    x='lon', y='lat', vmin=0, vmax=2, ax=ax)
ax.set_title('Terra MODIS')
ax.coastlines()
plt.savefig('Terra_MODIS_AOD.png', dpi=300)
plt.clf()

for model in an.models:
    # print(an.models[model].obj)
    regridder = an.model_regridders[model]
    # print(regridder)
    ds_model_regrid = regridder(an.models[model].obj)
    print(ds_model_regrid)

    ax = plt.subplot(projection=ccrs.PlateCarree())
    ds_model_regrid[varname_merra2_aod].isel(time=19).plot.pcolormesh(
        cmap=plt.cm.turbo,
        cbar_kwargs={'location': 'bottom', 'label': 'AOD'},
        x='lon', y='lat', vmin=0, vmax=2, ax=ax)
    ax.set_title(model)
    ax.coastlines()
    plt.savefig(model + '_AOD.png', dpi=300)

