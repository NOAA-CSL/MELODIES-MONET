from melodies_monet import driver

an = driver.analysis()
an.control = 'control_time_chunk.yaml'
an.read_control()

for n in range(len(an.time_intervals)-1):
    time_interval = [an.time_intervals[n], an.time_intervals[n+1]]

    an.open_models(time_interval=time_interval)
    an.open_obs(time_interval=time_interval)
    an.pair_data()

an.concat_pairs()

an.plotting()
an.stats()
