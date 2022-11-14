from melodies_monet import driver

an = driver.analysis()
an.control = 'control_time_chunk.yaml'
an.read_control()

for time_interval in an.control.time_intervals:

    an.open_models(time_interval=time_interval)
    an.open_obs(time_interval=time_interval)
    an.pair_data()

an.concat_pairs()

an.plotting()
an.stats()
