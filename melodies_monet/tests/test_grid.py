import xarray as xr

from glob import glob
from melodies_monet import driver

an = driver.analysis()

an.control = 'test_grid.yaml'
an.read_control()

files = glob(an.control_dict['obs']['test_obs']['files'])

for filename in files:
    print('reading ' + filename)
    obs_ds = xr.open_dataset(filename)
    print(obs_ds.info())
