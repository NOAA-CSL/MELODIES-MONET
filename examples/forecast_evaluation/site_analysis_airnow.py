#!/usr/bin/env python
# coding: utf-8

## Site specific airnow analysis
## get an array of sites with a certain amount of data available

import xarray as xr
import datetime as dt
import numpy as np
import re

#import variables from master script
import sys
species_name = sys.argv[1]
#print(species_name)
#species_name = 'OZONE'

#have to change the PM variable name
if species_name == 'PM25':
    species_name = 'PM2.5'

#define the region name and type to find the sites for
##could incorporate just finding the sites in specific regions into the bash script
region_name = ["R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10"] #R10, CA, etc
region_type = ["epa_region", "epa_region", "epa_region", "epa_region", "epa_region", "epa_region", "epa_region", "epa_region", "epa_region", "epa_region"] #epa_region, state

#open the obs data using xarray
data = xr.open_dataset('test5.nc') #this script needs to be after the test5 file has been created for the day but before the control yaml file has been created

#seperate out the site data 
site = data['site']

#define an array to read the sites into
sitesinregion = []
xregion = []

#loop through the specified regions
for r in range(len(region_name)):

    #define data for the specified region type
    region = data[region_type[r]]
    
    #loop through all the sites
    for i in range(len(site.values[0])):

        #select the region type data for one site and select only the values
        region_values = region.sel(x = i).values

        #If the site is in the specified region, append it to an array
        if region_values == region_name[r]:
            sitesinregion.append(site.sel(x=i).values[0])

            #append the x values to use to check if data is available for a specified species
            xregion.append(i)

#uncomment if you want an array of all sites in the region
#print(sitesinregion)

##selects the sites with data available for a specified species in a specifiec region
#selects only the data on the hour

#seperate out the species and time data
species = data[species_name]
time_local = data['time_local']
latitude = data['latitude']
longitude = data['longitude']

#define an array to read the sites into
#siteswdata = []

#for i in range(len(site.values[0])):
for i in xregion:

    #use the boolean to select the lat and lon corresponding to each site
    site_lat = latitude.sel(x=i).values[0]
    site_lon = longitude.sel(x=i).values[0]

    #skip sites that cause a parse error
    NON_PRINTABLE = re.compile('[^\x09\x0A\x0D\x20-\x7E\x85\xA0-\uD7FF\uE000-\uFFFD\U00010000-\U0010ffff]')
    match = NON_PRINTABLE.findall(site.sel(x=i).values[0])
    if match:
        continue

    #skip sites with a / in the name
    elif "/" in site.sel(x=i).values[0]:
        continue
   
    #skip the sites out of model domain
    elif site_lat < 27 or site_lat > 50 or site_lon < -130 or site_lon > -66:
        continue

    else:

        #selecting only the species values for one site on the hour
        ##create an empty array to read species values into
        species_values = []

        ##create an array of time values for one site
        time_array = time_local.sel(x = i).values

        ##loop through the times
        for j in range(len(time_array)):

            ##select an individual time from the array
            time = time_array[j][0]

            ##select the minute value
            t = dt.datetime.utcfromtimestamp(time.tolist() / 1e9)
            minute = t.minute

            #if the minute value is 0 append the species value to a list
            if minute == 0:
                species_values.append(species.sel(x = i).values[j])

        #convert the species values to an array
        species_values = np.array(species_values)

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
        if percent_data > 75:
#            siteswdata.append(site.sel(x=i).values[0])
            siteswdata = '"' + str(site.sel(x=i).values[0]) + '"'
            print(siteswdata)
        else:
            continue

#print(len(site_values))
#print(len(siteswdata))
#print(siteswdata)

