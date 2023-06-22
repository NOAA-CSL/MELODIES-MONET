from melodies_monet import driver

an = driver.analysis()
an.control = 'control_time_chunking_with_gridded_data.yaml'
an.read_control()
an.setup_regridders()

for time_interval in an.time_intervals:

    print(time_interval)

    an.open_obs(time_interval=time_interval)
    an.open_models(time_interval=time_interval)

    print(an.obs)
    for obs in an.obs:
        print(an.obs[obs].obj)
    print(an.models)
    for mod in an.models:
        print(an.models[mod].obj)
