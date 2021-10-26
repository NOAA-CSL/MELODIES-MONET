#!/bin/bash -x
shopt -s extglob
#Script to link model files to a common directory for processing with MELODIES-MONET

#Set new combined directory to link files into
dir_combined=/home/rschwantes/MONET/processed_plots/CMAQ_RRFS_SURF/08/data/wrf
#Set directory with raw model output data
dir_data=/wrk/d2/rschwantes/MONET_test/WRF-Chem/run_CONUS_fv19_BEIS_1.0xISO_RACM_v4.2.2_racm_vcp_noI
#Set year (needs to be four digits e.g., 2019)
year=2019
#Set month (needs to be two digits e.g., 08 or 10)
month=08
# go to local data directory
cd $dir_combined
# Update this to include days you want to include. For WRF-Chem, have to seperated 
# leading zeros because need to iterate. Also last day of month do seperately.
# To compare with CMAQ choose 12 to 12 but for others choose 06 to 00 
for dir in {1..2} 
do
  #ln -sf $dir_data/${month}0${dir}/wrfout_d01_${year}-${month}-0${dir}_06:00:00 .
  ln -sf $dir_data/${month}0${dir}/wrfout_d01_${year}-${month}-0${dir}_12:00:00 .
  ln -sf $dir_data/${month}0${dir}/wrfout_d01_${year}-${month}-0${dir}_18:00:00 .
  if [ $dir -eq 9 ]
  then
    ln -sf $dir_data/${month}${dir}/wrfout_d01_${year}-${month}-$(($dir+1))_00:00:00 .
    ln -sf $dir_data/${month}${dir}/wrfout_d01_${year}-${month}-$(($dir+1))_06:00:00 .
  else
    ln -sf $dir_data/${month}0${dir}/wrfout_d01_${year}-${month}-0$(($dir+1))_00:00:00 .
    ln -sf $dir_data/${month}0${dir}/wrfout_d01_${year}-${month}-0$(($dir+1))_06:00:00 .
  fi
done
#for dir in {10..30}
#do
#done
#Note for the last day of the month, you will need to manually select these to deal with
#change of month and possibly year if month is Dec.
#last_day=31
##ln -sf $dir_data/${month}${last_day}/wrfout_d01_${year}-${month}-${last_day}_06:00:00 .
#ln -sf $dir_data/${month}${last_day}/wrfout_d01_${year}-${month}-${last_day}_12:00:00 .
#ln -sf $dir_data/${month}${last_day}/wrfout_d01_${year}-${month}-${last_day}_18:00:00 .
#ln -sf $dir_data/${month}${last_day}/wrfout_d01_${year}-09-01_00:00:00 ..
#ln -sf $dir_data/${month}${last_day}/wrfout_d01_${year}-09-01_06:00:00 .
