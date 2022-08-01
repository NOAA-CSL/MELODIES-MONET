import sys
sys.path.append('..')
from melodies_monet import driver

an = driver.analysis()
an.control = 'mm_overview_wrfchem.yaml'
an.read_control()
an.control_dict

an.open_models()
an.models

an.open_obs()

an.pair_data()

an.plotting()
