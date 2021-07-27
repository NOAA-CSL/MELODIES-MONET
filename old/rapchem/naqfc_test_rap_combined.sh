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
yyyymmdd=2020-07-17
yyyymmdd48hr=2020-07-30
yyyymmddform=bau_all

# go to local data directory
data_dir=/home/rschwantes/RAPchem/analysis/RAP_plots/AGU/reanalysis/$yyyymmddform
mkdir $data_dir/$yyyymmddform
cd $data_dir
# link model output files
directory_output=/wrk/csd4/rahmadov/RAP-Chem/retro_agu_jul2020_v1/bau
for dir in {17..29}
do
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_01_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_02_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_03_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_04_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_05_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_06_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_07_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_08_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_09_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_10_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_11_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_12_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_13_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_14_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_15_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_16_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_17_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_18_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_19_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_20_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_21_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_22_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-${dir}_23_00_00_surface $yyyymmddform/.
ln -sf $directory_output/202007${dir}00/wrfout_d01_2020-07-$((dir+1))_00_00_00_surface $yyyymmddform/.
done 
# netcdf files (can be single or multiple files)
files=$yyyymmddform/* 
#echo $files
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
${stats} -p AIRNOW_RAP_${yyyymmdd}-01_${yyyymmdd48hr}-00_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 01:00:00" -ed "${yyyymmdd48hr} 00:00:00"
${stats} -p AIRNOW_RAP_${yyyymmdd}-01_${yyyymmdd48hr}-00_pair.hdf -s 'OZONE' 'PM2.5' -sd "${yyyymmdd} 01:00:00" -ed "${yyyymmdd48hr} 00:00:00"
${stats} -p AIRNOW_RAP_${yyyymmdd}-01_${yyyymmdd48hr}-00_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 01:00:00" -ed "${yyyymmdd48hr} 00:00:00" -b True -sn 'epa_region' -e 'R9'
${stats} -p AIRNOW_RAP_${yyyymmdd}-01_${yyyymmdd48hr}-00_pair.hdf -s 'OZONE' 'PM2.5' -sd "${yyyymmdd} 01:00:00" -ed "${yyyymmdd48hr} 00:00:00" -b True -sn 'epa_region' -e 'R9'
${stats} -p AIRNOW_RAP_${yyyymmdd}-01_${yyyymmdd48hr}-00_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 01:00:00" -ed "${yyyymmdd48hr} 00:00:00" -b True -sn 'epa_region' -e 'R4'
${stats} -p AIRNOW_RAP_${yyyymmdd}-01_${yyyymmdd48hr}-00_pair.hdf -s 'OZONE' 'PM2.5' -sd "${yyyymmdd} 01:00:00" -ed "${yyyymmdd48hr} 00:00:00" -b True -sn 'epa_region' -e 'R4'

# period spatial overlay plots
${spatial} -f ${files} -p AIRNOW_RAP_${yyyymmdd}-01_${yyyymmdd48hr}-00_pair.hdf -s 'OZONE' 'PM2.5' -sd "${yyyymmdd} 01:00:00" -ed "${yyyymmdd48hr} 00:00:00"

# period spatial bias plots
${bias} -p AIRNOW_RAP_${yyyymmdd}-01_${yyyymmdd48hr}-00_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 01:00:00" -ed "${yyyymmdd48hr} 00:00:00" 
${bias} -p AIRNOW_RAP_${yyyymmdd}-01_${yyyymmdd48hr}-00_pair.hdf -s 'OZONE' 'PM2.5' -sd "${yyyymmdd} 01:00:00" -ed "${yyyymmdd48hr} 00:00:00"

# period taylor plots
${taylor} -p AIRNOW_RAP_${yyyymmdd}-01_${yyyymmdd48hr}-00_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 01:00:00" -ed "${yyyymmdd48hr} 00:00:00"
# period box plots
${box} -p AIRNOW_RAP_${yyyymmdd}-01_${yyyymmdd48hr}-00_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 01:00:00" -ed "${yyyymmdd48hr} 00:00:00"

