#!/usr/bin/env python
# coding: utf-8

#Create a interactive HTML file displaying the plots for individual sites

import folium
import base64
from folium import IFrame
from PIL import Image
import numpy as np
import xarray as xr
import sys

#import variables from the bash file
species=sys.argv[1]
start_date=sys.argv[2]
end_date=sys.argv[3]
todays_date=sys.argv[4]

#species="PM25"
#start_date= "2022-08-15"
#end_date="2022-08-16"
#todays_date="20220816"

#have to change the PM variable name
if species == 'PM25':
    species = 'PM2.5'

#creating a folium map centered over the US
##could change the basemap here
fmap = folium.Map(location=[40, -95],zoom_start=4,min_zoom=3)

#defining the resolution, width, and height of the popup
resolution, width, height = 42, 12, 6

#path to plots
path = "/wrk/csd4/rahmadov/RAP-Chem/MELODIES-MONET-analysis/run_sites/plots_out/" + todays_date + "/" 

#defining plot type
plot_type = "grp1.timeseries"

#get the site names from the txt file
fileObj = open("sitefile.txt", "r") #opens the file in read mode
site_name_list = fileObj.read().splitlines() #puts the file into an array
fileObj.close()
##remove the extra characters from either end
site_name_list = [site[1:-1] for site in site_name_list]
##replace the spaces in the site names with underscores
site_name = [site.replace(' ','_') for site in site_name_list]

#defining arrays for the lat and lon
##reading in the airnow file
data = xr.open_dataset('test5.nc')

##creating arrays of lat, lon, and site names
latitude = data['latitude'].values[0]
longitude = data['longitude'].values[0]

site = data['site'].values[0]

##loop through the sites with plots
for i in range(len(site_name)):

    #boolean to identify the index of each site
    sitename_bool = site == site_name_list[i]

    #use the boolean to select the lat and lon corresponding to each site
    site_lat = latitude[sitename_bool]
    site_lon = longitude[sitename_bool]
    
    #define the plot name
    plot_name="plot_" + plot_type + "." + species + "." + start_date + "_00." + end_date + "_00.site." + site_name[i]

    #resizing the image to fit the popup window size and saving as another file
    image = Image.open(path + plot_name + ".png") #path to original image file
    image = image.resize(((width*resolution),(height*resolution)),Image.ANTIALIAS)
    image.save(fp = path + plot_name + "_resize.png") #save back to the original directory

    #opening the resized plot and encoding it
    png = path + plot_name + "_resize.png"
    encoded = base64.b64encode(open(png, 'rb').read())

    #formating the encoded plot as html
    html = '<img src="data:image/png;base64,{}">'.format

    #decoding the plot into an iframe
    iframe = IFrame(html(encoded.decode('UTF-8')), width=(width*resolution)+20, height=(height*resolution)+20)

    #creating a popup of the iframe
    popup = folium.Popup(iframe, max_width=2650)

    #creating the marker and adding it to the map
    ##sometimes complications with latitude and longitude for sites with the same name
    ##need to update so the correct latitude and longitude are always selected
    marker = folium.CircleMarker(location=[site_lat[-1], site_lon[-1]], popup=popup, radius=7,color = 'blue', fill_color="blue", fill_opacity=0.5)
    marker.add_to(fmap)

#saving the map to an html file
output_name = todays_date + "_" + species + "_fullhtml_sfc_f000.html"
fmap.save(output_name)


