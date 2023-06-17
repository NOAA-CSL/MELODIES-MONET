from melodies_monet import driver

an = driver.analysis()

an.control = '../yaml/control_wrfchem_aircraft_fromcharkincontrol_wrfchem_aircraft_vertprofile_Test_QZR.yaml'

an.read_control() 

an.control_dict   

an.open_models()

an.open_obs()

#This just pairs the data
an.pair_data()

an.plotting()
