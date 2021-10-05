#!/bin/bash -x
shopt -s extglob
#Script to link model files to a common directory for processing with MELODIES-MONET

#Set new combined directory to link files into
dir_combined=/scratch2/BMC/rcm1/rhs/fv3/regional/processed/surf_oper_rrfs_r131_v1/08_pm25/data/rrfs
#Set directory with raw model output data
dir_data=/scratch2/BMC/rcm1/rhs/fv3/regional/data/rrfs-cmaq/r131_v1
#Set year (needs to be four digits e.g., 2019)
year=2019
#Set month (needs to be two digits e.g., 08 or 10)
month=08
# go to local data directory
cd $dir_combined
# Update this to include days you want to include with 0 if needed e.g., {01..10}
# As you link, you need to rename the files to add year, month, day, so do not override. 
for dir in {01..10}
do
  for days in {01..24}
  do
    ln -sf $dir_data/${month}${dir}/dynf0${days}.nc ${year}${month}${dir}_dynf0${days}.nc
    ln -sf $dir_data/${month}${dir}/pm25/rcmaq.t12z.pm25tot_f0${days}.nc ${year}${month}${dir}_rcmaq.t12z.pm25tot_f0${days}.nc
  done
done
