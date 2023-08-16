from melodies_monet import driver

an = driver.analysis()

an.control = '../examples/yaml/control_wrfchem_aircraft_Latestfor_develop_aircraft.yaml'

an.read_control() 

an.control_dict   

an.open_models()

an.open_obs()

#This just pairs the data
an.pair_data()

an.plotting()
