from melodies_monet import driver

an = driver.analysis()
an.control = 'control_modis_l2.yaml'
an.read_control()


for time_interval in an.time_intervals:

    print(time_interval)

    an.open_obs(time_interval=time_interval)

    for obs in an.obs:
        print(an.obs[obs].obj)
