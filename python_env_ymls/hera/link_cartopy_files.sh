#!/bin/bash

####
# This file will link all the necessary shapefiles for cartopy plotting in matplotlib

#On hera and other computers with download restrictions, cartopy will error trying to
#download shapefiles for physical and cultural boundaries. This script solves this problem
#by manually linking the pre-downloaded cartopy shapefiles to your home directory
#(/home/user_name/.local/share/cartopy). Cartopy automatically checks this directory
#before trying to download the data.

#To run this script type: ./link_cartopy_files.sh
####

# local cartopy directory in home directory

mkdir -p ~/.local
mkdir -p ~/.local/share
mkdir -p ~/.local/share/cartopy/
mkdir -p ~/.local/share/cartopy/shapefiles
mkdir -p ~/.local/share/cartopy/natural_earth/natural_earth
mkdir -p ~/.local/share/cartopy/natural_earth/natural_earth/cultural
mkdir -p ~/.local/share/cartopy/natural_earth/natural_earth/physical


pushd ~/.local/share/cartopy/natural_earth/natural_earth/cultural

# Files are stored on rstprod for simplicity

orig_dir=/scratch1/NCEPDEV/rstprod/nexus_emissions/cartopy_shapefiles


# do cultural files first 
pushd ~/.local/share/cartopy/natural_earth/natural_earth/cultural
for res in 50 10 110; do
    ln -sf /scratch1/NCEPDEV/rstprod/nexus_emissions/cartopy_shapefiles/${res}m_cultural/* .
done 
popd
# now do the physical files 
pushd ~/.local/share/cartopy/natural_earth/natural_earth/physical
for res in 50 10 110; do
    ln -sf /scratch1/NCEPDEV/rstprod/nexus_emissions/cartopy_shapefiles/${res}m_physical/* .
done
popd 

mkdir -p ~/.local
mkdir -p ~/.local/share
mkdir -p ~/.local/share/cartopy/
mkdir -p ~/.local/share/cartopy/shapefiles
mkdir -p ~/.local/share/cartopy/shapefiles/natural_earth
mkdir -p ~/.local/share/cartopy/shapefiles/natural_earth/cultural
mkdir -p ~/.local/share/cartopy/shapefiles/natural_earth/physical

pushd ~/.local/share/cartopy/shapefiles/natural_earth/cultural
for res in 50 10 110; do
    ln -sf /scratch1/NCEPDEV/rstprod/nexus_emissions/cartopy_shapefiles/${res}m_cultural/* .
done
popd
# now do the physical files
pushd ~/.local/share/cartopy/shapefiles/natural_earth/physical
for res in 50 10 110; do
    ln -sf /scratch1/NCEPDEV/rstprod/nexus_emissions/cartopy_shapefiles/${res}m_physical/* .
done
popd
