#!/bin/bash -l

set -x 

case='2021REB'
indir=/scratch1/BMC/rcm2/jhe/model_output/sunvex21/run_CONUS_fv21_BEIS_0.75xISO_RACM_VCP_bdy/Output/
outdir1=/scratch2/BMC/rcm1/jhe/postp/MELODIES-MONET/AQS/simdata/$case/
outdir2=/scratch2/BMC/rcm1/jhe/postp/MELODIES-MONET/AQS/outdir_12kmCONUS_${case}/

year=2021
month_list=(04 05 06 07)
daymax_list=(30 31 30 31)
endday_list=(31 30 31 30)

length=${#month_list[@]} # number of elements in the array
echo $length

mkdir $outdir1
mkdir $outdir2

for (( ii = 0; ii < length; ii++ ))
do
    month=${month_list[$ii]}
    daymax=${daymax_list[$ii]}
    endday=${endday_list[$ii]}

    cd $outdir1
    mkdir $month
    cd $month
    for day in $(seq -f "%02g" 1 1 $daymax)
    do
        ln -sf $indir/$month$day/wrfout_d01_${year}-${month}-${day}_06:00:00 .
        ln -sf $indir/$month$day/wrfout_d01_${year}-${month}-${day}_12:00:00 .
        ln -sf $indir/$month$day/wrfout_d01_${year}-${month}-${day}_18:00:00 .
  
        if [ ${day#0} -eq 1 ]
        then
          pmonth=$((${month#0}-1))
          smonth=$(printf '%02d' $pmonth)
          ln -sf $indir/$smonth$endday/wrfout_d01_${year}-${month}-01_00:00:00 .
        fi

        if [ ${day#0} -eq $daymax ]
        then
          nmonth=$((${month#0}+1))
          echo $nmonth
          fmonth=$(printf '%02d' $nmonth)
          ln -sf $indir/$month$day/wrfout_d01_${year}-$fmonth-01_00:00:00 .
        else
          nday=$((${day#0}+1))
          echo $nday
          fday=$(printf '%02d' $nday)
          ln -sf $indir/$month$day/wrfout_d01_${year}-${month}-${fday}_00:00:00 .
        fi
 
    done
    cd ../..

    cd $outdir2
    mkdir $month
    mkdir $month/pm25

    cd ..

done

