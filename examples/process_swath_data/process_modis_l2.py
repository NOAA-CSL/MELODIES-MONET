from melodies_monet import driver

an = driver.analysis()
an.control = 'control_modis_l2.yaml'
an.read_control()

an.setup_obs_grid()
# print(an.obs_grid)

# an.setup_regridders()

for time_interval in an.time_intervals:

    print(time_interval)

    an.open_obs(time_interval=time_interval)
    an.update_obs_gridded_data()

an.normalize_obs_gridded_data()
print(an.obs_gridded_dataset)

for param in ['Terra_MODIS_AOD_550_Dark_Target_Deep_Blue_Combined',
              'Aqua_MODIS_AOD_550_Dark_Target_Deep_Blue_Combined']:
    param_data = an.obs_gridded_dataset[param + '_data'].values
    param_count = an.obs_gridded_dataset[param + '_count'].values
    mask = (param_count > 0)
    param_data[mask] = param_data[mask] / param_count[mask]
    an.obs_gridded_dataset[param + '_data'].values = param_data

an.obs_gridded_dataset.to_netcdf('MODIS.nc')

