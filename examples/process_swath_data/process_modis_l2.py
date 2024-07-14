from melodies_monet import driver

an = driver.analysis()
an.control = 'control_modis_l2.yaml'
an.read_control()

an.setup_obs_grid()
# print(an.obs_grid)

for time_interval in an.time_intervals:

    print(time_interval)
    an.open_obs(time_interval=time_interval)

    print(time_interval)
    an.update_obs_grid()

    """
    for obs in an.obs:
        for obs_time in an.obs[obs].obj:
            print('obs_time:', obs_time)
            print(an.obs[obs].obj[obs_time])
    """
