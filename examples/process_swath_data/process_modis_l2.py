from melodies_monet import driver

an = driver.analysis()
an.control = 'control_modis_l2.yaml'
an.read_control()

an.setup_obs_grid()
# print(an.obs_grid)

"""
for time_interval in an.time_intervals:

    print(time_interval)

    an.open_obs(time_interval=time_interval)
    an.update_obs_gridded_data()

an.normalize_obs_gridded_data()
"""
print(an.obs_gridded_dataset)
