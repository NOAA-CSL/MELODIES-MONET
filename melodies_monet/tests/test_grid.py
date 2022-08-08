from melodies_monet import driver

an = driver.analysis()

an.control = 'test_grid.yaml'
an.read_control()

print(an.control_dict)
