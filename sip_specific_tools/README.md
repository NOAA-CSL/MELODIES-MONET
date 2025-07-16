# Tools specifically designed for the SIP of Colorado

## Preprocessors
This preprocessors are designed for preprocessing non public datasets that are to be used by
the SIP of Colorado
### BoulderAir
BoulderAir data should be pre-processed to bo completely consistent with MELODIES-MONET's format.
This can be done by using the BoulderAir preprocessor, as follows

```
./boulderair.py -c coordinates_file.csv -p path_to_data_XYZ_etc.csv -v variable1,variable2,...,variablen -r resample_freq -m method -o output_name.nc
```

options: 

* `-c, --coordinates`

* `-p, --path`

* `-v, --variables`

* `-r, --resample_freq`

* `-m, --method`

* `-o, --output`



`coordinates_file.csv`: CSV containing the following fields: site_abbreviation,lat,lon,m_asl

`path_to_data_XYZ_etc.csv`: BoulderAir data. XYZ will be replaced automatically by the code with every site_abbreviation on the coordinates_file.csv. Wildcards (like `"*"` can be used).

`variablex`: Variables of interest.

`resample`(optional): If the data should be resampled. Options: `h`, `d`. If None is provided, the time in the files will be kept without resampling.

`method`: How to perform the resample. Default: `mean` -> The data will be averaged hourly. Other options: `median`, `max`, `min`, `inst`.

`output_name.nc`: Name of the output file. It will be netCDF


## Calculating ratios
Right now, proper ratio calculation is not implemented in the base code of MELODIES-MONET.
However, it can be manually done in the code, by defining objects (i.e., pair objects, although
specific observations/models objects can be defined) in the code.

An example of how to do this for each data type. This can possibly be done before the pairing, in
which case it would be paired automatically. However, the recommended procedure for satellite data
is to apply this directly to the paired object, to guarantee that the averaging kernels are properly
applied to the data before any operations are performed.

### Pair objects
```
ratio_hcho_no2 = xr.Dataset()
ratio_hcho_no2['model_hcho_no2'] = an.paired['tempo_l2_hcho_CAMX'].obj['FORM'] / an.paired['tempo_l2_no2_CAMX'].obj['NO2']
ratio_hcho_no2['tempo_hcho_no2'] = an.paired['tempo_l2_hcho_CAMX'].obj['vertical_column'] / an.paired['tempo_l2_no2_CAMX'].obj['vertical_column_troposphere']

p = driver.pair()
p.type = 'sat_swath_clm'
p.obs = 'tempo_l2_hcho_no2'
p.model = 'CAMX'
p.model_vars = ['model_hcho_no2']
p.obs_vars = ['tempo_hcho_no2']
p.obj = ratio_hcho_no2
p.filename = 'tempo_l2_hcho_no2_CAMX.nc'
an.paired['ratio_hcho_no2'] = p
```
