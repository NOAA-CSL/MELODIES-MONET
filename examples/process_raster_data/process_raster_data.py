from melodies_monet import driver

an = driver.analysis()
an.control = 'control_raster_data.yaml'
an.read_control()

for time_interval in an.time_intervals:

    print(time_interval)
