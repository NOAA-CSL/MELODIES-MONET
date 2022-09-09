#!/usr/bin/env python
# coding: utf-8

## Site specific aeronet analysis
## get an array of sites with a certain amount of data available

import xarray as xr

#import variables from master script
import sys
species_name = sys.argv[1]
#print(species_name)
#species_name = 'AOD_550'

#have to change the PM variable name
if species_name == 'AOD_550':
    species_name = 'aod_550nm'

#open the obs data using xarray
data = xr.open_dataset('test5.nc') #this script needs to be after the test5 file has been created for the day but before the control yaml file has been created

#define variables
site = data['siteid']
species = data[species_name]
latitude = data['latitude']
longitude = data['longitude']

for i in range(len(site.values)):

    #use the boolean to select the lat and lon corresponding to each site
    site_lat = latitude.sel(x=i).values
    site_lon = longitude.sel(x=i).values
    
    #skip the sites out of model domain
    if site_lat < 27 or site_lat > 49 or site_lon < -130 or site_lon > -66:
        continue

    else:
        
        #define species values for a single site
        species_values = species.sel(x = i).values

        #create an array of just the species values that aren't nan (-1 is nan in this dataset)
        species_notnan = species_values[species_values != -1]
        
        #finds the length of the notnull species data
        len_notnull = len(species_notnan)
        
        #finds the total length of data
        len_total = len(species_values)

        #finds the percent of the data avaiable
        if len_total > 0:
            percent_data = (len_notnull/len_total) * 100
        else:
            percent_data = 0

        #if more than 75% of the data is available print the site name
        if percent_data > 37.5:
#          siteswdata.append(site.sel(x=i).values[0])
            siteswdata = '"' + str(site.sel(x=i).values) + '"'
            print(siteswdata)
        else:
            continue


