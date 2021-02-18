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
yyyymmdd=2019-07-28
yyyymmdd48hr=2019-07-30
yyyymmddform=0728

# go to local data directory
data_dir=/home/rschwantes/RAPchem/analysis/FV3/CMAQ_OPER/oper/$yyyymmddform
mkdir $data_dir/$yyyymmddform
cd $data_dir
# link model output files
ln -sf /wrk/d2/rschwantes/cmaq/surface/CMAQ_OPER/gpfs_prod_aqm/$yyyymmddform/* $yyyymmddform/.

# netcdf files (can be single or multiple files)
#files=$yyyymmddform/"aqm.t12z.aconc-sfc.ncf"
files=$yyyymmddform/"aqm.t12z.aconc_sfc.ncf"

#naqfc AIRNOW verification scripts 
pair=/home/rschwantes/RAPchem/analysis/unchanged/naqfc_verify_scripts_cmaq/01.verify_pair.py
stats=/home/rschwantes/RAPchem/analysis/unchanged/naqfc_verify_scripts_cmaq/02.verify_stats.py
taylor=/home/rschwantes/RAPchem/analysis/unchanged/naqfc_verify_scripts_cmaq/03.verify_taylor_plots.py
bias=/home/rschwantes/RAPchem/analysis/unchanged/naqfc_verify_scripts_cmaq/04.verify_spatial_bias.py
spatial=/home/rschwantes/RAPchem/analysis/unchanged/naqfc_verify_scripts_cmaq/05.verify_spatial_plots.py
box=/home/rschwantes/RAPchem/analysis/unchanged/naqfc_verify_scripts_cmaq/06.verify_box_plots.py

# pair the data
${pair} -f ${files} -s 'OZONE' 'PM2.5' -p '/home/rschwantes/RAPchem/analysis/CMAQ_test/obs/'

# period statistics
${stats} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00"
${stats} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00"
${stats} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00" -b True -sn 'epa_region' -e 'R1'
${stats} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00" -b True -sn 'epa_region' -e 'R2'
${stats} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00" -b True -sn 'epa_region' -e 'R3'
${stats} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00" -b True -sn 'epa_region' -e 'R4'
${stats} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00" -b True -sn 'epa_region' -e 'R5'
${stats} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00" -b True -sn 'epa_region' -e 'R6'
${stats} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00" -b True -sn 'epa_region' -e 'R7'
${stats} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00" -b True -sn 'epa_region' -e 'R8'
${stats} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00" -b True -sn 'epa_region' -e 'R9'
${stats} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00" -b True -sn 'epa_region' -e 'R10'
# period spatial overlay plots
${spatial} -f ${files} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00" -miny 20.0 -maxy 60.0
${spatial} -f ${files} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'PM2.5' -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00" -miny 0.0 -maxy 20.0

# period spatial bias plots
${bias} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00" -miny -20.0 -maxy 20.0
${bias} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00" -miny -20.0 -maxy 20.0
${bias} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'PM2.5' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00" -miny -15.0 -maxy 15.0
${bias} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'PM2.5' -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00" -miny -15.0 -maxy 15.0

# period taylor plots
${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00"
#period box plots
${box} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00"
