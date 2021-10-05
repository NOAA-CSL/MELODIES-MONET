#!/bin/bash -x
shopt -s extglob
#Script to link model files to a common directory for processing with MELODIES-MONET

#Set new combined directory to link files into
dir_combined=/scratch2/BMC/rcm1/rhs/fv3/regional/processed/surf_oper_rrfs_r131_v1/08_pm25/data/cmaq_expt/
#Set directory with raw model output data
dir_data=/scratch2/BMC/rcm1/rhs/fv3/regional/data/cmaq_expt/fv3-cmaq-gbbpex1
#Set year (needs to be four digits e.g., 2019)
year=2019
#Set month (needs to be two digits e.g., 08 or 10)
month=08
# go to local data directory
cd $dir_combined
# Update this to include days you want to include with 0 if needed e.g., {01..10}
for dir in {01..10} 
do
  ln -sf $dir_data/${month}${dir}/aqm.${year}${month}${dir}.t12z.aconc-pm25_24.ncf .
done
#Note if CMAQ data includes more than first 24 hours, I run this command for each forecast
#ncks -d TSTEP,0,23 aqm.20190801.t12z.aconc-pm25.ncf aqm.20190801.t12z.aconc-pm25_24.ncf
