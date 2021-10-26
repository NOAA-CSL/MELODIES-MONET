#!/bin/bash -l

#SBATCH --job-name=mlitrop
#SBATCH --partition=hera
#SBATCH --time=08:00:00
# -- Request 16 cores
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
# -- Specify under which account a job should run
#SBATCH --account=rcm2


echo '=== Processing TROPOMI and WRF-Chem NO2 Columns ==='

export Basedir_tropomi='/scratch1/BMC/rcm2/mli/tropomi_data/NO2_global/'
export Baseoutdir='/scratch1/BMC/rcm2/mli/outdir_12km_noPM_baseline_bmc_cams/'
export Basedir_wrfoutput='/scratch1/BMC/rcm2/mli/nyc18_cams/run_12km_five18_bmcdVCP_fog_wofire_BEIS_0.5ISO/Output/'
export Geofile='/scratch1/BMC/rcm2/mli/nyc18_WPS/WPSV4.0/geo_em.d01.nc' # TO GET THE WRF boundaries

echo "---> Data Locations < ---"
echo "TROPOMI data are in:    " $Basedir_tropomi
echo "WRF-Chem data are in:   " $Basedir_wrfoutput
echo "Geophysical data are in:" $Geofile
echo "Out data locations:     " $Baseoutdir

year=2018
month=6
logdir=logs


for day in $(seq 1 1 15) 
do 
    echo "Processing TROPOMI data in " $year $month $day
    python CMP_WRFChem_TROPOMI_NO2Col_Daily_12km_Conservative_v3.py $year $month $day > $logdir/log_$year"_"$month"_"$day
done

