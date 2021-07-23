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
yyyymmdd=2020-09-26
yyyymmdd48hr=2020-09-28
yyyymmddform=2020092606

# go to local data directory
data_dir=/home/rschwantes/RAPchem/analysis/plots
mkdir $data_dir/$yyyymmddform
cd $data_dir
# link model output files
ln -sf /wrk/csd4/rahmadov/RAP-Chem/rapchem_covid/$yyyymmddform/* $yyyymmddform/.
 
# netcdf files (can be single or multiple files)
files=$yyyymmddform/!(*28_06*) 

#naqfc AIRNOW verification scripts 
pair=/home/rschwantes/RAPchem/analysis/unchanged/naqfc_verify_scripts_rapchem/01.verify_pair.py
stats=/home/rschwantes/RAPchem/analysis/unchanged/naqfc_verify_scripts_rapchem/02.verify_stats.py
taylor=/home/rschwantes/RAPchem/analysis/unchanged/naqfc_verify_scripts_rapchem/03.verify_taylor_plots.py
bias=/home/rschwantes/RAPchem/analysis/unchanged/naqfc_verify_scripts_rapchem/04.verify_spatial_bias.py
spatial=/home/rschwantes/RAPchem/analysis/unchanged/naqfc_verify_scripts_rapchem/05.verify_spatial_plots.py
box=/home/rschwantes/RAPchem/analysis/unchanged/naqfc_verify_scripts_rapchem/06.verify_box_plots.py

# pair the data
${pair} -f ${files} -s 'OZONE' 'PM2.5' -p '/home/rschwantes/RAPchem/analysis/CMAQ_test/obs/' 

# period statistics
${stats} -p AIRNOW_RAP_${yyyymmdd}-06_${yyyymmdd48hr}-05_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 06:00:00" -ed "${yyyymmdd48hr} 05:00:00"

# hourly spatial overlay plots
${spatial} -f ${files} -p AIRNOW_RAP_${yyyymmdd}-06_${yyyymmdd48hr}-05_pair.hdf -s 'OZONE' 'PM2.5'  -n ${yyyymmddform}
# period spatial overlay plots
${spatial} -f ${files} -p AIRNOW_RAP_${yyyymmdd}-06_${yyyymmdd48hr}-05_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 06:00:00" -ed "${yyyymmdd48hr} 05:00:00"


# hourly spatial bias plots
${bias} -p AIRNOW_RAP_${yyyymmdd}-06_${yyyymmdd48hr}-05_pair.hdf -s 'OZONE' 'PM2.5' -n ${yyyymmddform}
# period spatial bias plots
${bias} -p AIRNOW_RAP_${yyyymmdd}-06_${yyyymmdd48hr}-05_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 06:00:00" -ed "${yyyymmdd48hr} 05:00:00" 

# period taylor plots
${taylor} -p AIRNOW_RAP_${yyyymmdd}-06_${yyyymmdd48hr}-05_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 06:00:00" -ed "${yyyymmdd48hr} 05:00:00"
# period box plots
${box} -p AIRNOW_RAP_${yyyymmdd}-06_${yyyymmdd48hr}-05_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 06:00:00" -ed "${yyyymmdd48hr} 05:00:00"
