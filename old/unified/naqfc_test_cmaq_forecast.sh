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

#Specify Parameters
radius_of_influence=12000. #In meters
model='CMAQ'
mapping='{"OZONE":"O3","PM2.5":"PM25","PM10":"PM10","CO":"CO","NO":"NO","NO2":"NO2","SO2":"SO2","NOX":"NOX","NOY":"NOY","TEMP":"TEMP2","WS":"WSPD10","WD":"WDIR10","SRAD":"GSW","BARPR":"PRSFC","PRECIP":"RT","RHUM":"Q2"}'

# go to local data directory
data_dir=/home/rschwantes/RAPchem/analysis/FV3/CMAQ_OPER/oper_test/$yyyymmddform
mkdir $data_dir/$yyyymmddform
cd $data_dir
# link model output files
ln -sf /wrk/d2/rschwantes/cmaq/surface/CMAQ_OPER/gpfs_prod_aqm/$yyyymmddform/* $yyyymmddform/.

# netcdf files (can be single or multiple files)
#files=$yyyymmddform/"aqm.t12z.aconc-sfc.ncf"
files=$yyyymmddform/"aqm.t12z.aconc_sfc.ncf"

#naqfc AIRNOW verification scripts 
pair=/home/rschwantes/code/MONET-analysis/unified/01.verify_pair.py
stats=/home/rschwantes/code/MONET-analysis/unified/02.verify_stats.py
bias=/home/rschwantes/code/MONET-analysis/unified/04.verify_spatial_bias.py

# pair the data
${pair} -f ${files} -s 'OZONE' 'PM2.5' -p '/home/rschwantes/RAPchem/analysis/CMAQ_test/obs/' -r ${radius_of_influence} -m ${model} -map ${mapping}

# period statistics
${stats} -p AIRNOW_${model}_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' 'PM2.5' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00" -map ${mapping}
# period spatial bias plots
${bias} -p AIRNOW_${model}_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf -s 'OZONE' -r True -sd "${yyyymmdd} 13:00:00" -ed "${yyyymmdd48hr} 12:00:00" -miny -20.0 -maxy 20.0 -m ${model} -map ${mapping} -n "${model}_AIRNOW"
