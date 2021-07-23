#!/bin/bash -x

#Example script to invoke all MONETv2 AIRNOW Evaluation package scripts

source /data/aqf2/barryb/anaconda2/bin/activate monetdev

# get the 48-hour date to current
#date
#yyyymmdd=$(date -d "-2 days" +%Y-%m-%d)
#yyyymmdd48hr=$(date -d "-0 days" +%Y-%m-%d)
#yyyymmddform=$(date -d "-2 days" +%Y%m%d)

#Or, Set explicit 48 hr dates
yyyymmdd=2019-04-09
yyyymmdd48hr=2019-04-11
yyyymmddform=20190409

# go to local data directory
data_dir=/data/aqf3/patrickc/5xpm
mkdir $data_dir/$yyyymmddform
cd $data_dir
# link model output files
ln -sf /data/aqf2/barryb/5xpm/$yyyymmddform/* $yyyymmddform/.
 
# netcdf files (can be single or multiple files)
files=$yyyymmddform/"aqm.t12z.aconc-sfc.ncf"

#naqfc AIRNOW verification scripts 
pair=/data/aqf/patrickc/naqfc_verify_scripts/01.verify_pair.py
stats=/data/aqf/patrickc/naqfc_verify_scripts/02.verify_stats.py
taylor=/data/aqf/patrickc/naqfc_verify_scripts/03.verify_taylor_plots.py
bias=/data/aqf/patrickc/naqfc_verify_scripts/04.verify_spatial_bias.py
spatial=/data/aqf/patrickc/naqfc_verify_scripts/05.verify_spatial_plots.py
box=/data/aqf/patrickc/naqfc_verify_scripts/06.verify_box_plots.py

# pair the data
${pair} -f ${files} -s 'OZONE' 'PM2.5' 

# period statistics
${stats} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00"

# hourly spatial overlay plots
#${spatial} -f ${files} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5'  -n ${yyyymmddform}
# period spatial overlay plots
${spatial} -f ${files} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00"

# hourly spatial bias plots
#${bias} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -n ${yyyymmddform}
# period spatial bias plots
${bias} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00" 

# hourly taylor plots
#${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -n ${yyyymmddform}
# period taylor plots
${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00"
