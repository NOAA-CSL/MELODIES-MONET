#!/bin/bash -x
shopt -s extglob
#added this so I can skip the last file which is a different size
#Example script to invoke all MONETv2 AIRNOW Evaluation package scripts

source /home/rschwantes/anaconda3/bin/activate py36

# get the 48-hour date to current
#date
#yyyymmdd=$(date -d "-2 days" +%Y-%m-%d)
#yyyymmdd48hr=$(date -d "-0 days" +%Y-%m-%d)
#yyyymmddform=$(date -d "-2 days" +%Y%m%d)

#Or, Set explicit 48 hr dates
yyyymmdd=2020-07-21
yyyymmdd48hr=2020-08-01
yyyymmddform=20200721

# go to local data directory
data_dir=/home/rschwantes/RAPchem/analysis/WRF_plots/AGU/conus_bau_base
mkdir $data_dir/$yyyymmddform
cd $data_dir
# link model output files
directory_output=/wrk/d2/rschwantes/wrf/covidaqs/run_CONUS_fv19_BEIS_0.5xISO 
for dir in {21..30}
do
ln -sf $directory_output/07${dir}/wrfout_d01_2020-07-${dir}_06:00:00 $yyyymmddform/.
ln -sf $directory_output/07${dir}/wrfout_d01_2020-07-${dir}_12:00:00 $yyyymmddform/.
ln -sf $directory_output/07${dir}/wrfout_d01_2020-07-${dir}_18:00:00 $yyyymmddform/.
ln -sf $directory_output/07${dir}/wrfout_d01_2020-07-$(($dir+1))_00:00:00 $yyyymmddform/.
done
ln -sf $directory_output/0731/wrfout_d01_2020-07-31_06:00:00 $yyyymmddform/.
ln -sf $directory_output/0731/wrfout_d01_2020-07-31_12:00:00 $yyyymmddform/.
ln -sf $directory_output/0731/wrfout_d01_2020-07-31_18:00:00 $yyyymmddform/.
ln -sf $directory_output/0731/wrfout_d01_2020-08-01_00:00:00 $yyyymmddform/.

# netcdf files (can be single or multiple files)
files=$yyyymmddform/* 

#naqfc AIRNOW verification scripts 
pair=/home/rschwantes/RAPchem/analysis/unchanged/naqfc_verify_scripts_wrfchem/01.verify_pair.py
stats=/home/rschwantes/RAPchem/analysis/unchanged/naqfc_verify_scripts_wrfchem/02.verify_stats.py
taylor=/home/rschwantes/RAPchem/analysis/unchanged/naqfc_verify_scripts_wrfchem/03.verify_taylor_plots.py
bias=/home/rschwantes/RAPchem/analysis/unchanged/naqfc_verify_scripts_wrfchem/04.verify_spatial_bias.py
spatial=/home/rschwantes/RAPchem/analysis/unchanged/naqfc_verify_scripts_wrfchem/05.verify_spatial_plots.py
box=/home/rschwantes/RAPchem/analysis/unchanged/naqfc_verify_scripts_wrfchem/06.verify_box_plots.py

## pair the data
${pair} -f ${files} -s 'OZONE' -p '/home/rschwantes/RAPchem/analysis/CMAQ_test/obs/' 

## period statistics
${stats} -p AIRNOW_RAP_${yyyymmdd}-06_${yyyymmdd48hr}-05_pair.hdf -s 'OZONE' -r True -sd "${yyyymmdd} 06:00:00" -ed "${yyyymmdd48hr} 05:00:00"
${stats} -p AIRNOW_RAP_${yyyymmdd}-06_${yyyymmdd48hr}-05_pair.hdf -s 'OZONE' -sd "${yyyymmdd} 06:00:00" -ed "${yyyymmdd48hr} 05:00:00"
${stats} -p AIRNOW_RAP_${yyyymmdd}-06_${yyyymmdd48hr}-05_pair.hdf -s 'OZONE' -r True -sd "${yyyymmdd} 06:00:00" -ed "${yyyymmdd48hr} 05:00:00" -b True -sn 'epa_region' -e 'R9' 
${stats} -p AIRNOW_RAP_${yyyymmdd}-06_${yyyymmdd48hr}-05_pair.hdf -s 'OZONE' -sd "${yyyymmdd} 06:00:00" -ed "${yyyymmdd48hr} 05:00:00" -b True -sn 'epa_region' -e 'R9'
${stats} -p AIRNOW_RAP_${yyyymmdd}-06_${yyyymmdd48hr}-05_pair.hdf -s 'OZONE' -r True -sd "${yyyymmdd} 06:00:00" -ed "${yyyymmdd48hr} 05:00:00" -b True -sn 'epa_region' -e 'R4'
${stats} -p AIRNOW_RAP_${yyyymmdd}-06_${yyyymmdd48hr}-05_pair.hdf -s 'OZONE' -sd "${yyyymmdd} 06:00:00" -ed "${yyyymmdd48hr} 05:00:00" -b True -sn 'epa_region' -e 'R4'

## period spatial overlay plots
${spatial} -f ${files} -p AIRNOW_RAP_${yyyymmdd}-06_${yyyymmdd48hr}-05_pair.hdf -s 'OZONE' -sd "${yyyymmdd} 06:00:00" -ed "${yyyymmdd48hr} 05:00:00"

## period spatial bias plots
${bias} -p AIRNOW_RAP_${yyyymmdd}-06_${yyyymmdd48hr}-05_pair.hdf -s 'OZONE' -r True -sd "${yyyymmdd} 06:00:00" -ed "${yyyymmdd48hr} 05:00:00" -miny -15.0 -maxy 15.0 

## period box plots
${box} -p AIRNOW_RAP_${yyyymmdd}-06_${yyyymmdd48hr}-05_pair.hdf -s 'OZONE' -r True -sd "${yyyymmdd} 06:00:00" -ed "${yyyymmdd48hr} 05:00:00"

