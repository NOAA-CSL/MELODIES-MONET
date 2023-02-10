#!/bin/bash
#
# This file will link all the necessary shapefiles for cartopy plotting in matplotlib
#
# On hera and other computers with download restrictions, cartopy will error trying to
# download shapefiles for physical and cultural boundaries. This script solves this problem
# by manually linking the pre-downloaded cartopy shapefiles to your home directory
# (/home/<user name>/.local/share/cartopy). Cartopy automatically checks this directory
# before trying to download the data.
#
# To run this script type: ./link_cartopy_files.sh

# Location of Natural Earth shapefiles
# Directory structure: {res}_{cultural,physical} e.g. '110m_physical'
src_base=/scratch1/RDARCH/rda-arl-gpu/Barry.Baker/emissions/nexus/cartopy_shapefiles

# Cartopy user data dir, where cartopy checks for shapefiles by default
cartopy_data_dir=~/.local/share/cartopy

# The files are organized by resolution and type in the source dir
ress=( 10m 50m 110m )
types=( cultural physical )

# Link
for res in ${ress[*]}; do
  for type_ in ${types[*]}; do
    src=$src_base/${res}_${type_}
    dst=$cartopy_data_dir/shapefiles/natural_earth/$type_
    mkdir -p $dst
    files=$(find $src -maxdepth 1 -type f -not -path '*/.*' -not -name 'Archive.zip' 2> /dev/null)
    if [ "$files" != '' ]; then
      echo "Linking $(echo $files | wc -w) files (res=$res, type=$type_)"
      echo "  from $src"
      echo "  to   $dst"
      ln -sf $files $dst/
    else
      echo "error: no files found in $src"
      echo check that that dir exists and you have read permissions
      exit 1
    fi
  done
done
