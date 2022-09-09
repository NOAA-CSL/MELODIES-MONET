#/bin/bash -li

DATE=/bin/date

#--------------------------------------------------------------------------------------------------------------------------------------------------------
# || Model Name           ||      Initialization Times (*=used) ||  Run Length  ||  Relative Day to compare to AirNow || Color Indicator  ||     Marker
#---------------------------------------------------------------------------------------------------------------------------------------------------------
# || RAP-Chem (1)         ||         6Z*, 18Z*, day = N - 2    ||   36 hours   ||    N (06Z), N-1 (18Z)               || dodger blue      ||     x/x
# || WRF-Chem (firex)     ||         0                         ||   48 Hours   ||    N - 1                            || green            ||      x
# || WRF-Chem (AQ-WATCH)  ||         0                         ||   48 Hours   ||    N - 1                            || greenyellow      ||      x
# || NAQFC (cmaq_opr)     ||         0,6*,12*,18               ||   6,48,48,6  ||    N - 1                            || red              ||     x/x
# || NAQFC (cmaq_exp)     ||            ""                     ||      ""      ||    N (06Z), N-1 (12Z)               || magenta          ||     x/x
# || RAQMS                ||         0                         ||              ||                                     ||                  ||
# || WACCM                ||
# || HRRR-Smoke           ||            
# || onlinecmaq_b1        ||         12*                       ||   48         ||    N-1
# || onlinecmaq_c3
# || RAP-Smoke
# || HRRR-Chem
# || .../

# Date Variables YYYYMMDD format
# The HH should only ever be 00Z; i.e., we only run one comparsion per day starting at the top of the day
todays_date=`${DATE} +%Y%m%d`                     # today, YYYYMMDD format
ytdays_date=`${DATE} +%Y%m%d -d " 24 hours ago"`  # yesterday
t2days_date=`${DATE} +%Y%m%d -d " 48 hours ago"`  # two days ago
t3days_date=`${DATE} +%Y%m%d -d " 72 hours ago"`  # 3 days ago
t4days_date=`${DATE} +%Y%m%d -d " 96 hours ago"`  # 4 days ago
t5days_date=`${DATE} +%Y%m%d -d " 120 hours ago"`  # 5 days ago

# Today variables 
YYYY_today=`${DATE} +%Y`
MM_today=`${DATE} +%m`
DD_today=`${DATE} +%d`
HH_today=`${DATE} +%H`
# Yesterday Variables
YYYY_ytday=`${DATE} +%Y -d " 24 hours ago"`
MM_ytday=`${DATE} +%m -d " 24 hours ago"`
DD_ytday=`${DATE} +%d -d " 24 hours ago"`
HH_ytday=`${DATE} +%H -d " 24 hours ago"`
# Two days ago variables
YYYY_t2day=`${DATE} +%Y -d " 48 hours ago"`
MM_t2day=`${DATE} +%m -d " 48 hours ago"`
DD_t2day=`${DATE} +%d -d " 48 hours ago"`
HH_t2day=`${DATE} +%H -d " 48 hours ago"`
# 3 days ago
YYYY_t3day=`${DATE} +%Y -d " 72 hours ago"`
MM_t3day=`${DATE} +%m -d " 72 hours ago"`
DD_t3day=`${DATE} +%d -d " 72 hours ago"`
HH_t3day=`${DATE} +%H -d " 72 hours ago"`
# 4 days ago
YYYY_t4day=`${DATE} +%Y -d " 96 hours ago"`
MM_t4day=`${DATE} +%m -d " 96 hours ago"`
DD_t4day=`${DATE} +%d -d " 96 hours ago"`
HH_t4day=`${DATE} +%H -d " 96 hours ago"`
# 5 days ago
YYYY_t5day=`${DATE} +%Y -d " 120 hours ago"`
MM_t5day=`${DATE} +%m -d " 120 hours ago"`
DD_t5day=`${DATE} +%d -d " 120 hours ago"`
HH_t5day=`${DATE} +%H -d " 120 hours ago"`

# start/end time is for yesterday, day N-1
start_time_yaml=${YYYY_ytday}-${MM_ytday}-${DD_ytday}-00:00:00
end_time_yaml=`${DATE} -d "${YYYY_ytday}${MM_ytday}${DD_ytday} + 24 hours" +%Y-%m-%d-%H:00:00`
#start_time_yaml=${YYYY_t4day}-${MM_t4day}-${DD_t4day}-00:00:00
#end_time_yaml=`${DATE} -d "${YYYY_t4day}${MM_t4day}${DD_t4day} + 24 hours" +%Y-%m-%d-%H:00:00`

# Dates also need to be for yesterday -- see here I think it only makes sense to go from 00-23
start_time_reformat=`${DATE} -d "${YYYY_ytday}${MM_ytday}${DD_ytday}" +%Y-%m-%d`
end_time_reformat=`${DATE} -d "${YYYY_ytday}${MM_ytday}${DD_ytday} + 24 hours" +%Y-%m-%d`
#start_time_reformat=`${DATE} -d "${YYYY_t4day}${MM_t4day}${DD_t4day}" +%Y-%m-%d`
#end_time_reformat=`${DATE} -d "${YYYY_t4day}${MM_t4day}${DD_t4day} + 24 hours" +%Y-%m-%d`

# Print to verfiy  we are starting the script for today
echo "Today is $todays_date" 
echo "Attempting to run MONET Analysis over the period from ${start_time_reformat} to ${end_time_reformat}"


# Start hour for each of the chosen forecasts (for now, only once intialization can be analyzed at a time)
#------------------------------------------------
start_HH_rapchem1="18"
start_HH_rapchem2="06"
start_HH_rapchem3="06"

start_HH_wrfchem_firexaq="00"
start_HH_wrfchem_aqwatch="00"
start_HH_wrfchem_firexaq2="00"
start_HH_wrfchem_aqwatch2="00"

start_HH_cmaq_oper06="06"
start_HH_cmaq_oper12="12"
start_HH_cmaq_oper12_2="12"

start_HH_cmaq_expr12="12"
start_HH_cmaq_expr06="06"

start_HH_hrrr_smoke06="06"
start_HH_hrrr_smoke18="18"

start_HH_geos_cf12="12"

start_HH_onlinecmaq_b1="12"
start_HH_onlinecmaq_c3="12"

start_HH_raqms="00"

start_HH_rapsmoke="06"

start_HH_hrrrchem="00"
#------------------------------------------------

#Region Choices? 
#------------------------------------------------
#CONUS=1; R1=1; R2=1; R3=1; R4=1; R5=1; R6=1; R7=1; R8=1; R9=1; R10=1

# .. Species choices ..
#------------------------------------------------
#species_list=("AOD_550" "CO")
species_list=("PM25" "PM10" "OZONE" "NO2" "CO" "TEMP" "PRECIP" "AOD_550")
#species_list=("OZONE" "CO" "TEMP")
#PM25=0; OZONE=0; WS=0; PRSFC=0; PRECIP=0; TEMP=0; CO=0; SO2=0; NO=0; NO2=0
ns=${#species_list[@]}
echo "Will loop through variables: ${species_list[@]}"

# .. Other namelist
tol_hours_missing=0  #tolerance threshold for number of hours acceptable to be missing 
mdl_lw=1.5

# Loop over each of the species and perform the analysis individually
for is in $( seq 0 $(((${ns} - 1))) )
do
	species=${species_list[$is]}

	echo "Now working on analysis for $species"
if [[ "${species}" == "PM10" ]]; then
	CONUS=0; R1=0; R2=0; R3=0; R4=0; R5=0; R6=0; R7=0; R8=0; R9=0; R10=0; site=1
#        CONUS=1; R1=1; R2=0; R3=0; R4=1; R5=1; R6=1; R7=1; R8=1; R9=1; R10=1
elif [[ "${species}" == "NO2" ]]; then
        CONUS=0; R1=0; R2=0; R3=0; R4=0; R5=0; R6=0; R7=0; R8=0; R9=0; R10=0; site=1
#        CONUS=1; R1=0; R2=1; R3=1; R4=0; R5=1; R6=0; R7=0; R8=0; R9=1; R10=0
elif [[ "${species}" == "CO" ]]; then
	CONUS=0; R1=0; R2=0; R3=0; R4=0; R5=0; R6=0; R7=0; R8=0; R9=0; R10=0; site=1
#        CONUS=1; R1=0; R2=0; R3=1; R4=0; R5=1; R6=0; R7=0; R8=0; R9=1; R10=0
elif [[ "${species}" == "OZONE" ]]; then
        CONUS=0; R1=0; R2=0; R3=0; R4=0; R5=0; R6=0; R7=0; R8=0; R9=0; R10=0; site=1
#        CONUS=1; R1=1; R2=1; R3=1; R4=1; R5=1; R6=1; R7=1; R8=0; R9=1; R10=1
else
	CONUS=0; R1=0; R2=0; R3=0; R4=0; R5=0; R6=0; R7=0; R8=0; R9=0; R10=0; site=1
#        CONUS=1; R1=1; R2=1; R3=1; R4=1; R5=1; R6=1; R7=1; R8=1; R9=1; R10=1	
fi

# choose which models to use
do_rapchem1=1
do_rapchem2=1
do_rapchem3=0

if [[ "${species}" != "PRECIP" ]]; then
do_wrfchem_firexaq=1
do_wrfchem_aqwatch=1
do_wrfchem_firexaq2=0
do_wrfchem_aqwatch2=0
else
do_wrfchem_firexaq=0
do_wrfchem_aqwatch=0
do_wrfchem_firexaq2=0
do_wrfchem_aqwatch2=0
fi

if [[ "${species}" != "AOD_550" && "${species}" != "CO" && "${species}" != "PM10" && "${species}" != "TEMP" && "${species}" != "PRECIP" ]]; then
do_cmaq_oper12=1
do_cmaq_oper06=1
do_cmaq_oper12_2=0
do_cmaq_expr12=1
do_cmaq_expr06=1
else
do_cmaq_oper12=0
do_cmaq_oper06=0
do_cmaq_oper12_2=0
do_cmaq_expr12=0
do_cmaq_expr06=0
fi

if [[ "${species}" != "PM25" && "${species}" != "TEMP" && "${species}" != "PM10" && "${species}" != "PRECIP" ]]; then
do_raqms=0
else
do_raqms=0
fi

do_waccm=0

if [[ "${species}" == "PM25" || "${species}" == "AOD_550" || "${species}" == "PRECIP" || "${species}" == "TEMP" ]]; then
#if [[ "${species}" == "PM25" || "${species}" == "AOD_550" ]]; then
do_hrrrsmoke06=1
do_hrrrsmoke18=1
else
do_hrrrsmoke06=0
do_hrrrsmoke18=0
fi

do_gefsaero=0

if [[ "${species}" != "PM10" && "${species}" != "TEMP"  && "${species}" != "PRECIP" ]]; then
do_geos_cf12=0
else
do_geos_cf12=0
fi

if [[ "${species}" != "AOD_550" && "${species}" != "PM10" && "${species}" != "TEMP" && "${species}" != "PRECIP" ]]; then
do_onlinecmaq_b1=1
do_onlinecmaq_c3=1
else
do_onlinecmaq_b1=0
do_onlinecmaq_c3=0
fi

if [[ "${species}" == "PM25" || "${species}" == "PRECIP" || "${species}" == "TEMP" ]]; then
#if [[ "${species}" == "PM25" ]]; then
do_rapsmoke=1
else
do_rapsmoke=0
fi

do_hrrrchem=1

if [[ ${species} == "AOD_550" ]]; then
do_stats=0
ts_select_time="'time'"
else	
do_stats=0
ts_select_time="'time_local'"
fi

# Loop over the plot types, the spatial plots take  a lot of time and may fail 
#for ip in $( seq 0 1 )
#do 
#	if [[ ${ip} == 0 ]]; then
#		timeseries=1
#		taylor=1
#		spatial_bias=0
#		spatial_overlay=0
#		boxplot=1
#	elif [[ ${ip} == 1 ]]; then
#		timeseries=0
#                taylor=0
#                spatial_bias=1
#                spatial_overlay=0
#                boxplot=0
#	elif [[ ${ip} == 2 ]]; then
#		timeseries=0
#                taylor=0
#                spatial_bias=0
#                spatial_overlay=1
#                boxplot=0
#	fi

timeseries=1
#taylor=1
#spatial_bias=0
#spatial_overlay=0
#boxplot=1

######################################################################
## Add in regions for site analysis ##
# first need to create the appropriate test5 file

# First remove any test5.nc
rm -f test5.nc

#call and run reformat_airnow_rapchemtest.py with appropriate variables
if [[ "${species}" == "AOD_550" ]]; then
python reformat_aeronet_rapchemtest.py ${start_time_reformat} ${end_time_reformat}
else
python reformat_airnow_rapchemtest.py ${start_time_reformat} ${end_time_reformat}
fi


# After we run the reformat script, move the test5 file (if it worked) to todays date and then relink it

# First check to see if a file (should be a link) named "test5.nc" exists
 if [[ -e test5.nc ]]; then
         if [[ "${species}" == "AOD_550" ]]; then
               mv test5.nc test5.aeronet.nc.${todays_date}
         else
               mv test5.nc test5.nc.${todays_date}   #moving the test5.nc generated from reformat to test5.nc.{todays_date}
         fi
 else
        echo "test5.nc is missing!"
        exit 1
 fi
 if [[ "${species}" == "AOD_550" ]]; then
       ln -sf test5.aeronet.nc.${todays_date} test5.nc
 else
       ln -sf test5.nc.${todays_date} test5.nc   #creating symbolic test5.nc file pointing to test5.nc.{todays_date}
 fi

## now need to run the python script to get a list of sites with available data

if [ ${site} -eq 1 ]; then
	if [[ "${species}" == "AOD_550" ]]; then
	echo "running site analysis for AOD"
	python site_analysis_aeronet.py ${species} > sitefile.txt
	else
	python site_analysis_airnow.py ${species} > sitefile.txt
	fi
readarray -t site_list < sitefile.txt
#rm sitefile.txt

nsite=${#site_list[@]}
#echo "Including data from the following sites: ${site_list}"
fi


#build array of region choices to go into control.yaml as RGN_list 
if [ ${CONUS} -eq 1 ]; then
	echo "Adding plots/stats for CONUS"
  if [ ${#rgn_list[@]} -eq 0 ]; then
         rgn_list[0]="'CONUS'"
	 rgn_type[0]="'all'"
  else
         rgn_list[ ${#rgn_list[@]} +1 ]=",'CONUS'"
	 rgn_type[ ${#rgn_type[@]} +1 ]=",'all'"
  fi
fi

if [[ "${species}" != "AOD_550" ]]; then # No EPA Regions for AERONET
for r in $(seq 1 10)
do
	key="R${r}"
	eval do_region='$'$key
	if [ ${do_region} -eq 1 ]; then
        if [ ${#rgn_list[@]} -eq 0 ]; then
         rgn_list[0]="'$key'"
         rgn_type[0]="'epa_region'"
         else
         rgn_list[ ${#rgn_list[@]} +1 ]=",'$key'"
         rgn_type[ ${#rgn_type[@]} +1 ]=",'epa_region'"
        fi
	fi
done
fi

#add the sites to the region list
if [ ${site} -eq 1 ]; then
        echo "Adding plots/stats for sites"
	if [[ "${species}" == "AOD_550" ]]; then
        echo "adding sites for AOD"
	for isite in $( seq 0 $(((${nsite} - 1))) )
	do
        	 if [ ${#rgn_list[@]} -eq 0 ]; then
                	rgn_list[0]="${site_list[$isite]}"
                 	rgn_type[0]="'siteid'"
        	 else
                	# echo ${site_list[$isite]} 
                 	rgn_list[ ${#rgn_list[@]} +1 ]=",${site_list[$isite]}"
                 	rgn_type[ ${#rgn_type[@]} +1 ]=",'siteid'"
        	 fi
	done
	else
	for isite in $( seq 0 $(((${nsite} - 1))) )
	do
        	 if [ ${#rgn_list[@]} -eq 0 ]; then
                	 rgn_list[0]="${site_list[$isite]}"
                	 rgn_type[0]="'site'"
         	else
                	# echo ${site_list[$isite]} 
                 	rgn_list[ ${#rgn_list[@]} +1 ]=",${site_list[$isite]}"
                 	rgn_type[ ${#rgn_type[@]} +1 ]=",'site'"
        	 fi
	done
	fi
fi

echo "Including regions: ${rgn_list[*]}"

start_time_rapchem1=${YYYY_t2day}${MM_t2day}${DD_t2day} 
start_time_rapchem2=${YYYY_ytday}${MM_ytday}${DD_ytday} 
start_time_rapchem3=${YYYY_t2day}${MM_t2day}${DD_t2day}
start_time_wrfchem_firexaq=${YYYY_ytday}${MM_ytday}${DD_ytday} 
start_time_wrfchem_aqwatch=${YYYY_ytday}${MM_ytday}${DD_ytday} 
start_time_wrfchem_firexaq2=${YYYY_t2day}${MM_t2day}${DD_t2day}
start_time_wrfchem_aqwatch2=${YYYY_t2day}${MM_t2day}${DD_t2day}
start_time_cmaq_oper12_2=${YYYY_t3day}${MM_t3day}${DD_t3day} 
start_time_cmaq_oper06=${YYYY_t2day}${MM_t2day}${DD_t2day} 
start_time_cmaq_oper12=${YYYY_t2day}${MM_t2day}${DD_t2day}
start_time_cmaq_expr12=${YYYY_t2day}${MM_t2day}${DD_t2day}
start_time_cmaq_expr06=${YYYY_t2day}${MM_t2day}${DD_t2day} 
start_time_hrrr_smoke06=${YYYY_t2day}${MM_t2day}${DD_t2day} 
start_time_hrrr_smoke18=${YYYY_t2day}${MM_t2day}${DD_t2day} 
start_time_geos_cf12=${YYYY_t2day}${MM_t2day}${DD_t2day}
start_time_onlinecmaq_b1=${YYYY_t2day}${MM_t2day}${DD_t2day} 
start_time_onlinecmaq_c3=${YYYY_t2day}${MM_t2day}${DD_t2day}
start_time_raqms=${YYYY_ytday}${MM_ytday}${DD_ytday}
start_time_rapsmoke=${YYYY_ytday}${MM_ytday}${DD_ytday} 
start_time_hrrrchem=${YYYY_ytday}${MM_ytday}${DD_ytday} 

# Set the data directories based on the time
rapchem1_datadir=/wrk/csd4/rahmadov/RAP-Chem/rapchem_covid/${start_time_rapchem1}${start_HH_rapchem1}
rapchem2_datadir=/wrk/csd4/rahmadov/RAP-Chem/rapchem_covid/${start_time_rapchem2}${start_HH_rapchem2}
rapchem3_datadir=/wrk/csd4/rahmadov/RAP-Chem/rapchem_covid/${start_time_rapchem3}${start_HH_rapchem3}
wrfchem_firexaq_datadir=/wrk/csd4/rahmadov/RAP-Chem/ncar_firexaq/${start_time_wrfchem_firexaq}${start_HH_wrfchem_firexaq}
wrfchem_aqwatch_datadir=/wrk/csd4/rahmadov/RAP-Chem/ncar_aqwatch/${start_time_wrfchem_aqwatch}${start_HH_wrfchem_aqwatch}
wrfchem_firexaq2_datadir=/wrk/csd4/rahmadov/RAP-Chem/ncar_firexaq/${start_time_wrfchem_firexaq2}${start_HH_wrfchem_firexaq2}
wrfchem_aqwatch2_datadir=/wrk/csd4/rahmadov/RAP-Chem/ncar_aqwatch/${start_time_wrfchem_aqwatch2}${start_HH_wrfchem_aqwatch2}
cmaq_oper12_2_datadir=/wrk/csd4/rahmadov/RAP-Chem/cmaq_oper/${start_time_cmaq_oper12_2}${start_HH_cmaq_oper12_2}
cmaq_oper12_datadir=/wrk/csd4/rahmadov/RAP-Chem/cmaq_oper/${start_time_cmaq_oper12}${start_HH_cmaq_oper12}
cmaq_expr12_datadir=/wrk/csd4/rahmadov/RAP-Chem/cmaq_expr/${start_time_cmaq_expr12}${start_HH_cmaq_expr12}
cmaq_oper06_datadir=/wrk/csd4/rahmadov/RAP-Chem/cmaq_oper/${start_time_cmaq_oper06}${start_HH_cmaq_oper06}
cmaq_expr06_datadir=/wrk/csd4/rahmadov/RAP-Chem/cmaq_expr/${start_time_cmaq_expr06}${start_HH_cmaq_expr06}
hrrr_smoke06_datadir=/wrk/csd4/rahmadov/RAP-Chem/hrrr_smoke/${start_time_hrrr_smoke06}${start_HH_hrrr_smoke06}
hrrr_smoke18_datadir=/wrk/csd4/rahmadov/RAP-Chem/hrrr_smoke/${start_time_hrrr_smoke18}${start_HH_hrrr_smoke18}
geos_cf12_datadir=/wrk/csd4/rahmadov/RAP-Chem/geos-cf/${start_time_geos_cf12}${start_HH_geos_cf12}
onlinecmaq_b1_datadir=/wrk/csd4/rahmadov/RAP-Chem/online_cmaq_expr_b1/${start_time_onlinecmaq_b1}${start_HH_onlinecmaq_b1}
onlinecmaq_c3_datadir=/wrk/csd4/rahmadov/RAP-Chem/online_cmaq_expr_v70c3/${start_time_onlinecmaq_c3}${start_HH_onlinecmaq_c3}
raqms_datadir=/wrk/csd4/rahmadov/RAP-Chem/raqms/${start_time_raqms}${start_HH_raqms}
rapsmoke_datadir=/wrk/csd4/rahmadov/RAP-Chem/rap_smoke/${start_time_rapsmoke}${start_HH_rapsmoke}
hrrrchem_datadir=/wrk/csd4/rahmadov/RAP-Chem/hrrrchem_covid/${start_time_hrrrchem}${start_HH_hrrrchem}

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

# CHECKING FOR MODEL DATA
# Check to see the model data is available and set some flags 
# Initiialie the flags
has_rapchem1=0
has_rapchem2=0
has_rapchem3=0
has_wrfchem_firexaq=0
has_wrfchem_aqwatch=0
has_wrfchem_firexaq2=0
has_wrfchem_aqwatch2=0
has_cmaq_oper12=0
has_cmaq_oper12_2=0
has_cmaq_expr12=0
has_cmaq_oper06=0
has_cmaq_expr06=0
has_raqms=0
has_waccm=0
has_hrrrsmoke06=0
has_hrrrsmoke18=0
has_gefsaero=0
has_geos_cf=0
has_onlinecmaq_b1=0
has_onlinecmaq_c3=0
has_rapsmoke=0
has_hrrrchem=0

missing_rapchem1=0
 if [[ ${do_rapchem1} -eq 1 ]]; then
  for i in $(seq 0 29) 
  do
	  datestr_rapchem=`${DATE} -d "${start_time_rapchem1} ${start_HH_rapchem1} + $i hours" +%Y-%m-%d_%H`
	  rapchem_testfile=${rapchem1_datadir}/wrfout_d01_${datestr_rapchem}_00_00_surface
	  if [[ ! -e ${rapchem_testfile} ]]; then
		  missing_rapchem1=$((missing_rapchem1+1))
	  fi
  done
  if [[ ${missing_rapchem1} -gt ${tol_hours_missing} ]]; then
	  echo "Missing ${missing_rapchem1} files for RAP-Chem, more than allowed (${tol_hours_missing}), not processing RAP-Chem (1)"  
  else
	  echo "Found sufficient files to process RAP-Chem (1)"
	has_rapchem1=1
  fi
 fi
## ---------------------------------
#
 missing_rapchem2=0
 if [[ ${do_rapchem2} -eq 1 ]]; then
  for i in $(seq 0 17)
  do
          datestr_rapchem=`${DATE} -d "${start_time_rapchem2} ${start_HH_rapchem2} + $i hours" +%Y-%m-%d_%H`
          rapchem_testfile=${rapchem2_datadir}/wrfout_d01_${datestr_rapchem}_00_00_surface
          if [[ ! -e ${rapchem_testfile} ]]; then
		  echo "Missing ${rapchem_testfile}"
                  missing_rapchem2=$((missing_rapchem2+1))
          fi
  done
  if [[ ${missing_rapchem2} -gt ${tol_hours_missing} ]]; then
          echo "Missing ${missing_rapchem2} files for RAP-Chem, more than allowed (${tol_hours_missing}), not processing RAP-Chem (2)"
  else
          echo "Found sufficient files to process RAP-Chem (2)"
        has_rapchem2=1
  fi
 fi
## ---------------------------------
#

 missing_rapchem3=0
 if [[ ${do_rapchem3} -eq 1 ]]; then
  for i in $(seq 0 30)
  do
          datestr_rapchem=`${DATE} -d "${start_time_rapchem3} ${start_HH_rapchem3} + $i hours" +%Y-%m-%d_%H`
          rapchem_testfile=${rapchem3_datadir}/wrfout_d01_${datestr_rapchem}_00_00_surface
          if [[ ! -e ${rapchem_testfile} ]]; then
                  echo "Missing ${rapchem_testfile}"
                  missing_rapchem3=$((missing_rapchem2+1))
          fi
  done
  if [[ ${missing_rapchem3} -gt ${tol_hours_missing} ]]; then
          echo "Missing ${missing_rapchem3} files for RAP-Chem, more than allowed (${tol_hours_missing}), not processing RAP-Chem (3)"
  else
          echo "Found sufficient files to process RAP-Chem (3)"
        has_rapchem3=1
  fi
 fi
## ---------------------------------
#

 missing_wrfchem_firexaq=0
 if [[ ${do_wrfchem_firexaq} -eq 1 ]]; then
  for i in $(seq 0 23) 
  do
	  datestr_wrfchem=`${DATE} -d "${start_time_wrfchem_firexaq} ${start_HH_wrfchem_firexaq} + $i hours" +%Y-%m-%d_%H`
	  wrfchem_testfile=${wrfchem_firexaq_datadir}/wrfout_hourly_d01_${datestr_wrfchem}:00:00_surface
	  if [[ ! -e ${wrfchem_testfile} ]]; then
		  missing_wrfchem_firexaq=$((missing_wrfchem_firexaq+1))
	  fi
  done
  if [[ ${missing_wrfchem_firexaq} -gt ${tol_hours_missing} ]]; then
	echo "Missing ${missing_wrfchem_firexaq} files for WRF-Chem (FIREX-AQ), more than allowed (${tol_hours_missing}), not processing WRF-Chem (FIREX-AQ)"  
   else
	   echo "Found sufficent files to process NCAR_FIREX-AQ"
	has_wrfchem_firexaq=1
  fi
 fi
#-------------------------------------------------
#
 missing_wrfchem_aqwatch=0
 if [[ ${do_wrfchem_aqwatch} -eq 1 ]]; then
  for i in $(seq 0 23)
  do
          datestr_wrfchem=`${DATE} -d "${start_time_wrfchem_aqwatch} ${start_HH_wrfchem_aqwatch} + $i hours" +%Y-%m-%d_%H`
          wrfchem_testfile=${wrfchem_aqwatch_datadir}/wrfout_hourly_d01_${datestr_wrfchem}:00:00_surface
          if [[ ! -e ${wrfchem_testfile} ]]; then
                  missing_wrfchem_aqwatch=$((missing_wrfchem_aqwatch+1))
          fi
  done
if [[ ${missing_wrfchem_aqwatch} -gt ${tol_hours_missing} ]]; then
        echo "Missing ${missing_wrfchem_aqwatch} files for WRF-Chem (AQ-WATCH), more than allowed (${tol_hours_missing}), not processing WRF-Chem (AQ-WATCH)"  
   else
	   echo "Found sufficient files to process NCAR_AQ-WATCH"
        has_wrfchem_aqwatch=1
fi
fi


## ---------------------------------
#
 missing_wrfchem_firexaq2=0
 if [[ ${do_wrfchem_firexaq2} -eq 1 ]]; then
  for i in $(seq 0 47) 
  do
	  datestr_wrfchem=`${DATE} -d "${start_time_wrfchem_firexaq2} ${start_HH_wrfchem_firexaq2} + $i hours" +%Y-%m-%d_%H`
	  wrfchem_testfile=${wrfchem_firexaq2_datadir}/wrfout_hourly_d01_${datestr_wrfchem}:00:00_surface
	  if [[ ! -e ${wrfchem_testfile} ]]; then
		  missing_wrfchem_firexaq2=$((missing_wrfchem_firexaq2+1))
	  fi
  done
  if [[ ${missing_wrfchem_firexaq2} -gt ${tol_hours_missing} ]]; then
	echo "Missing ${missing_wrfchem_firexaq2} files for WRF-Chem (FIREX-AQ, N-1), more than allowed (${tol_hours_missing}), not processing WRF-Chem (FIREX-AQ, N-1)"  
   else
	   echo "Found sufficent files to process NCAR_FIREX-AQ_N-1"
	has_wrfchem_firexaq2=1
  fi
 fi
#-------------------------------------------------
#
 missing_wrfchem_aqwatch2=0
 if [[ ${do_wrfchem_aqwatch2} -eq 1 ]]; then
  for i in $(seq 0 47)
  do
          datestr_wrfchem=`${DATE} -d "${start_time_wrfchem_aqwatch2} ${start_HH_wrfchem_aqwatch2} + $i hours" +%Y-%m-%d_%H`
          wrfchem_testfile=${wrfchem_aqwatch2_datadir}/wrfout_hourly_d01_${datestr_wrfchem}:00:00_surface
          if [[ ! -e ${wrfchem_testfile} ]]; then
                  missing_wrfchem_aqwatch2=$((missing_wrfchem_aqwatch2+1))
          fi
  done
if [[ ${missing_wrfchem_aqwatch2} -gt ${tol_hours_missing} ]]; then
        echo "Missing ${missing_wrfchem_aqwatch2} files for WRF-Chem (AQ-WATCH, N-1), more than allowed (${tol_hours_missing}), not processing WRF-Chem (AQ-WATCH, N-1)"  
   else
	   echo "Found sufficient files to process NCAR_AQ-WATCH_N-1"
        has_wrfchem_aqwatch2=1
fi
fi
#-------------------------------------------------
if [[ $do_cmaq_oper12 -eq 1 ]]; then
#checking for cmaq_oper
cmaq_oper_testfile=${cmaq_oper12_datadir}/aqm.t12z.aconc.ncf
# Now check for the files, setting flag to 1 if at least one rapchem file is present
if [[ -r ${cmaq_oper_testfile} ]]; then
	has_cmaq_oper12=1
	echo "Found NAQFC Operational file for 12Z, day N -1 "
else
	echo "Cannot find NAQFC (operational, 12Z, day N - 1 ) file, not processing"
fi
fi

if [[ $do_cmaq_oper12_2 -eq 1 ]]; then
#checking for cmaq_oper
cmaq_oper_testfile=${cmaq_oper12_2_datadir}/aqm.t12z.aconc.ncf
# Now check for the files, setting flag to 1 if at least one rapchem file is present
if [[ -r ${cmaq_oper_testfile} ]]; then
        has_cmaq_oper12_2=1
        echo "Found NAQFC Operational file for 12Z, day N -2 "
else
        echo "Cannot find NAQFC (operational, 12Z, day N - 2 ) file, not processing"
fi
fi

if [[ ${do_cmaq_oper06} -eq 1 ]]; then
cmaq_oper_testfile=${cmaq_oper06_datadir}/aqm.t06z.aconc.ncf
# Now check for the files, setting flag to 1 if at least one rapchem file is present
if [[ -r ${cmaq_oper_testfile} ]]; then
        has_cmaq_oper06=1
	echo "Found NAQFC Operational file for 06Z, day N - 1"
else
        echo "Cannot find NAQFC (operational, 06Z, day N - 1 ) file, not processing"
fi
fi

if [[ ${do_cmaq_expr12} -eq 1 ]]; then
#checking for cmaq_expr
cmaq_expr_testfile=${cmaq_expr12_datadir}/aqm.t12z.aconc_sfc.ncf
# Now check for the files, setting flag to 1 if at least one rapchem file is present
if [[ -r ${cmaq_expr_testfile} ]]; then
	has_cmaq_expr12=1
	echo "Found NAQFC Experimental file for 12Z, day N -1"
else
	echo "Cannot find NAQFC (experimental, 12Z) file, not processing"
fi
fi

if [[ ${do_cmaq_expr06} -eq 1 ]]; then
cmaq_expr_testfile=${cmaq_expr06_datadir}/aqm.t06z.aconc_sfc.ncf
# Now check for the files, setting flag to 1 if at least one rapchem file is present
if [[ -r ${cmaq_expr_testfile} ]]; then
        has_cmaq_expr06=1
	echo "Found NAQFC Experimental file for 06Z, day N-1"
else
        echo "Cannot find NAQFC (experimental, 06Z) file, not processing"
fi
fi

# HRRR-SMOKE
if [[ ${do_hrrrsmoke06} -eq 1 ]]; then
	hrrrsmoke_testfile=${hrrr_smoke06_datadir}/smoke.all.${YYYY_t2day}${MM_t2day}${DD_t2day}.nc 
if [[ -r ${hrrrsmoke_testfile} ]]; then
	has_hrrrsmoke06=1
	echo "Found HRRR-Smoke file for 06Z, day N-1"
else
	echo "Cannot find HRRR-Smoke 06Z forecast file, not processing"
fi
fi
# HRRR-SMOKE
if [[ ${do_hrrrsmoke18} -eq 1 ]]; then
        hrrrsmoke_testfile=${hrrr_smoke18_datadir}/smoke.all.${YYYY_t2day}${MM_t2day}${DD_t2day}.nc
if [[ -r ${hrrrsmoke_testfile} ]]; then
        has_hrrrsmoke18=1
        echo "Found HRRR-Smoke file for 18Z, day N-1"
else
        echo "Cannot find HRRR-Smoke 18Z forecast file, not processing"
fi      
fi

if [[ ${do_geos_cf12} -eq 1 ]]; then
	geos_cf_testfile=${geos_cf12_datadir}/geos-cf.${YYYY_t2day}${MM_t2day}${DD_t2day}12.nc
	echo "geos_cf_testfile = ${geos_cf_testfile}"
	if [[ -r ${geos_cf_testfile} ]]; then
		has_geos_cf12=1
		echo "Found GEOS-CF for 12Z, day N-1"
	else
                echo "Cannot find GEOS-CF 12Z forecast file, not processing"
        fi
fi

# onlinecmaq_b1 
if [[ ${do_onlinecmaq_b1} -eq 1 ]]; then
#checking for onlinecmaq_b1
for i in {01..72}
  do
          onlinecmaq_b1_testfile=${onlinecmaq_b1_datadir}/aqm.t12z.chem_3d.f0${i}.nc
          if [[ ! -e ${onlinecmaq_b1_testfile} ]]; then
                  missing_onlinecmaq_b1=$((missing_onlinecmaq_b1+1))
          fi
  done

  if [[ ${missing_onlinecmaq_b1} -gt ${tol_hours_missing} ]]; then
        echo "Missing ${missing_onlinecmaq_b1} files for onlinecmaq_b1, more than allowed (${tol_hours_missing}), not processing onlinecmaq_b1"  
   else
           echo "Found sufficent files to process onlinecmaq_b1"
        has_onlinecmaq_b1=1
  fi
fi

# onlinecmaq_c3
if [[ ${do_onlinecmaq_c3} -eq 1 ]]; then
#checking for onlinecmaq_c3
for i in {01..72}
  do
          onlinecmaq_c3_testfile=${onlinecmaq_c3_datadir}/aqm.t12z.chem_3d.f0${i}.nc
          if [[ ! -e ${onlinecmaq_c3_testfile} ]]; then
                  missing_onlinecmaq_c3=$((missing_onlinecmaq_c3+1))
          fi
  done

  if [[ ${missing_onlinecmaq_c3} -gt ${tol_hours_missing} ]]; then
        echo "Missing ${missing_onlinecmaq_c3} files for onlinecmaq_c3, more than allowed (${tol_hours_missing}), not processing onlinecmaq_c3"  
   else
           echo "Found sufficent files to process onlinecmaq_c3"
        has_onlinecmaq_c3=1
  fi
fi

#RAQMS
missing_raqms=0
if [[ ${do_raqms} -eq 1 ]]; then
for i in $(seq 0 7)
  do
          datestr_raqms=`${DATE} -d "${start_time_raqms} ${start_HH_raqms} + $((6*i)) hours" +%m_%d_%Y_%H`
          raqms_testfile=${raqms_datadir}/uwhyb_${datestr_raqms}Z.chem.assim.nc
          if [[ ! -e ${raqms_testfile} ]]; then
                  missing_raqms=$((missing_raqms+1))
          fi
  done
if [[ ${missing_raqms} -gt ${tol_hours_missing} ]]; then
        echo "Missing ${missing_raqms} files for RAQMS, more than allowed (${tol_hours_missing}), not processing RAQMS"
else
           echo "Found sufficient files to process RAQMS"
        has_raqms=1
fi
fi

# RAP-SMOKE
if [[ ${do_rapsmoke} -eq 1 ]]; then
        rapsmoke_testfile=${rapsmoke_datadir}/smoke.all.${YYYY_t2day}${MM_t2day}${DD_t2day}.nc
if [[ -r ${rapsmoke_testfile} ]]; then
        has_rapsmoke=1
        echo "Found RAP-Smoke file"
else
        echo "Cannot find RAP-Smoke forecast file, not processing"
fi
fi

# HRRR-CHEM
missing_hrrrchem=0
 if [[ ${do_hrrrchem} -eq 1 ]]; then
  for i in $(seq 0 23)
  do
          datestr_hrrrchem=`${DATE} -d "${start_time_hrrrchem} ${start_HH_hrrrchem} + $i hours" +%Y-%m-%d_%H`
          hrrrchem_testfile=${hrrrchem_datadir}/wrfout_d01_${datestr_hrrrchem}_00_00_surface_meiyu
	  if [[ ! -e ${hrrrchem_testfile} ]]; then
                  missing_hrrrchem=$((missing_hrrrchem+1))
          fi
  done
  if [[ ${missing_hrrrchem} -gt ${tol_hours_missing} ]]; then
          echo "Missing ${missing_hrrrchem} files for HRRR-Chem, more than allowed (${tol_hours_missing}), not processing HRRR-Chem"
  else
          echo "Found sufficient files to process HRRR-Chem"
        has_hrrrchem=1
  fi
 fi

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

#... Create model list arrays dependent on model choices 
# rapchem
if [[ ${do_rapchem2} -eq 1 && ${has_rapchem2} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then
         mdl_list[0]="RAP-Chem_06Z"
  else
         mdl_list[ ${#mdl_list[@]} +1 ]="RAP-Chem_06Z"
  fi
fi
if [[ ${do_rapchem1} -eq 1 && ${has_rapchem1} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then
	 mdl_list[0]="RAP-Chem_18Z"
  else
	 mdl_list[ ${#mdl_list[@]} +1 ]="RAP-Chem_18Z"
  fi
fi
if [[ ${do_rapchem3} -eq 1 && ${has_rapchem3} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then
         mdl_list[0]="RAP-Chem_06Z_N-1"
  else
         mdl_list[ ${#mdl_list[@]} +1 ]="RAP-Chem_06Z_N-1"
  fi
fi
#wrfchem-firexaq 
if [[ ${do_wrfchem_firexaq} -eq 1 && ${has_wrfchem_firexaq} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then
          mdl_list[0]="NCAR_FIREX-AQ" #"wrfchem_v4.0"
  else
	 mdl_list[ $((${#mdl_list[@]}+1)) ]="NCAR_FIREX-AQ" #"wrfchem_v4.0"
  fi
fi
#wrfchem - aq-watch
if [[ ${do_wrfchem_aqwatch} -eq 1 && ${has_wrfchem_aqwatch} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then
          mdl_list[0]="NCAR_AQ-WATCH" #"wrfchem_v4.0"
  else
         mdl_list[ $((${#mdl_list[@]}+1)) ]="NCAR_AQ-WATCH" #"wrfchem_v4.0"
  fi
fi
#wrfchem-firexaq2 
if [[ ${do_wrfchem_firexaq2} -eq 1 && ${has_wrfchem_firexaq2} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then
          mdl_list[0]="NCAR_FIREX-AQ_N-1" #"wrfchem_v4.0"
  else
         mdl_list[ $((${#mdl_list[@]}+1)) ]="NCAR_FIREX-AQ_N-1" #"wrfchem_v4.0"
  fi
fi
#wrfchem - aq-watch2
if [[ ${do_wrfchem_aqwatch2} -eq 1 && ${has_wrfchem_aqwatch2} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then
          mdl_list[0]="NCAR_AQ-WATCH_N-1" #"wrfchem_v4.0"
  else
         mdl_list[ $((${#mdl_list[@]}+1)) ]="NCAR_AQ-WATCH_N-1" #"wrfchem_v4.0"
  fi
fi
#cmaq_oper
if [[ ${do_cmaq_oper06} -eq 1 && ${has_cmaq_oper06} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then
         mdl_list[0]="NAQFC_CMAQ_oper_06Z"
  else
         mdl_list[ ${#mdl_list[@]} +1 ]="NAQFC_CMAQ_oper_06Z"
  fi
fi
if [[ ${do_cmaq_oper12} -eq 1 && ${has_cmaq_oper12} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then
         mdl_list[0]="NAQFC_CMAQ_oper_12Z"
  else
         mdl_list[ ${#mdl_list[@]} +1 ]="NAQFC_CMAQ_oper_12Z"
  fi
fi
if [[ ${do_cmaq_oper12_2} -eq 1 && ${has_cmaq_oper12_2} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then
         mdl_list[0]="NAQFC_CMAQ_oper_12Z_N-2"
  else
         mdl_list[ ${#mdl_list[@]} +1 ]="NAQFC_CMAQ_oper_12Z_N-2"
  fi
fi
#cmaq_expr
if [[ ${do_cmaq_expr06} -eq 1 && ${has_cmaq_expr06} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then
         mdl_list[0]="NAQFC_CMAQ_expr_06Z"
  else
          mdl_list[ ${#mdl_list[@]} +1 ]="NAQFC_CMAQ_expr_06Z"
  fi
fi
#cmaq_expr
if [[ ${do_cmaq_expr12} -eq 1 && ${has_cmaq_expr12} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then 
	  mdl_list[0]="NAQFC_CMAQ_expr_12Z"
  else
          mdl_list[ ${#mdl_list[@]} +1 ]="NAQFC_CMAQ_expr_12Z"
  fi
fi
# HRRR-Smoke
if [[ ${do_hrrrsmoke06} -eq 1 && ${has_hrrrsmoke06} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then
         mdl_list[0]="HRRR-Smoke_06Z"
  else
          mdl_list[ ${#mdl_list[@]} +1 ]="HRRR-Smoke_06Z"
  fi
fi
if [[ ${do_hrrrsmoke18} -eq 1 && ${has_hrrrsmoke18} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then
         mdl_list[0]="HRRR-Smoke_18Z"
  else
          mdl_list[ ${#mdl_list[@]} +1 ]="HRRR-Smoke_18Z"
  fi
fi

if [[ ${do_geos_cf12} -eq 1 && ${has_geos_cf12} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then
         mdl_list[0]="GEOS-CF_12Z"
  else
         mdl_list[ ${#mdl_list[@]} +1 ]="GEOS-CF_12Z"
  fi
fi

#onlinecmaq_b1
if [[ ${do_onlinecmaq_b1} -eq 1 && ${has_onlinecmaq_b1} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then
          mdl_list[0]="online-CMAQ(UFS,v7.0b1)"
  else
          mdl_list[ ${#mdl_list[@]} +1 ]="online-CMAQ(UFS,v7.0b1)"
  fi
fi

#onlinecmaq_c3
if [[ ${do_onlinecmaq_c3} -eq 1 && ${has_onlinecmaq_c3} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then
          mdl_list[0]="online-CMAQ(UFS,v7.0c3)"
  else
          mdl_list[ ${#mdl_list[@]} +1 ]="online-CMAQ(UFS,v7.0c3)"
  fi
fi

#raqms
if [[ ${do_raqms} -eq 1 && ${has_raqms} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then
          mdl_list[0]="RAQMS"
  else
          mdl_list[ ${#mdl_list[@]} +1 ]="RAQMS"
  fi
fi

# RAP-Smoke
if [[ ${do_rapsmoke} -eq 1 && ${has_rapsmoke} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then
         mdl_list[0]="RAP-Smoke"
  else
          mdl_list[ ${#mdl_list[@]} +1 ]="RAP-Smoke"
  fi
fi

# HRRR-CHEM
if [[ ${do_hrrrchem} -eq 1 && ${has_hrrrchem} -eq 1 ]]; then
  if [ ${#mdl_list[@]} -eq 0 ]; then
         mdl_list[0]="HRRR-Chem_00Z"
  else
         mdl_list[ ${#mdl_list[@]} +1 ]="HRRR-Chem_00Z"
  fi
fi

# .. Observations ..
#for now we just have airnow, can create namelist and code similar to above to add in more i.e balloon soundings 
if [[ ${species} == "AOD_550" ]]; then
      obs_list[0]="aeronet"
else
      obs_list[0]="airnow"
fi

# check that at least one model is being used 
if [ ${#mdl_list[@]} -eq 0 ]; then
	echo "Model list empty, please choose at least one (e.g., rapchem)"
	exit 1
elif [ ${#obs_list[@]} -eq 0 ]; then
	echo "Obs list empty, please choose at least one (e.g., airnow)"
	exit 1
else
	echo "Proceeding using the obs: ${obs_list[*]} and models ${mdl_list[*]}"
fi 

# Construct obs-mdl pair for plot groups
knt=0
for iobs in ${obs_list[@]}; do
for imdl in ${mdl_list[@]}; do
	if [ $knt -eq 0 ]; then
		obs_mdl_list[knt++]="'${iobs}_${imdl}'"
        else
		obs_mdl_list[knt++]=",'${iobs}_${imdl}'"
	fi
done
done


#Building Control Yaml file based of namelists above
# General Description:
# Any key that is specific for a plot type will begin with ts for timeseries, ty for taylor
# Opt: Specifying the variable or variable group is optional
# For now all plots except time series average over the analysis window. 
# Seting axis values - If set_axis = True in data_proc section of each plot_grp the yaxis for the plot will be set based on the values specified in the obs section for each variable. If set_axis is set to False, then defaults will be used. 'vmin_plot' and 'vmax_plot' are needed for 'timeseries', 'spatial_overlay', and 'boxplot'. 'vdiff_plot' is needed for 'spatial_bias' plots and'ty_scale' is needed for 'taylor' plots. 'nlevels' or the number of levels used in the contour plot can also optionally be provided for spatial_overlay plot. If set_axis = True and the proper limits are not provided in the obs section, a warning will print, and the plot will be created using the default limits.
# ... Construct the control.yaml file/namelist
echo "Bulding control.yaml"
rm -f control.yaml.${todays_date}
cat << EOF > control.yaml.${todays_date}
analysis:
  start_time: '${start_time_yaml}'
  end_time:  '${end_time_yaml}'    #UTC #
  debug: True
EOF


# Now insert which models

cat << EOF >> control.yaml.${todays_date}
model:
EOF

#RAP Chem 
if [[ ${do_rapchem1} -eq 1 && ${has_rapchem1} -eq 1 ]]; then
echo "Including RAP-Chem 18Z"
cat << EOF >> control.yaml.${todays_date}
  RAP-Chem_18Z: # model label
    files: ${rapchem1_datadir}/*
    mod_type: 'wrfchem'
    mod_kwargs:
      surf_only_nc: True
    radius_of_influence: 12000 #meters
    mapping: #model species name : obs species name
EOF
if [[ "${species}" == "AOD_550" ]]; then
cat << EOF >> control.yaml.${todays_date}
      aeronet:
EOF
else
cat << EOF >> control.yaml.${todays_date}
      airnow:
EOF
fi
if [[ "${species}" == "CO" ]]; then
	cat << EOF  >> control.yaml.${todays_date} 
        co: CO
EOF
fi
if [[ "${species}" == "NO2" ]]; then
	cat << EOF >> control.yaml.${todays_date}
        no2: NO2
EOF
fi
if [[ "${species}" == "PM25" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM2_5_DRY: "PM2.5"
EOF
fi
if [[ "${species}" == "PM10" ]]; then
	cat << EOF >> control.yaml.${todays_date}
        PM10: "PM10"
EOF
fi
if [[ "${species}" == "OZONE" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        o3: "OZONE"
EOF
fi
if [[ "${species}" == "AOD_550" ]]; then
	cat << EOF >> control.yaml.${todays_date}
        AOD550: "aod_550nm"
EOF
fi
if [[ "${species}" == "TEMP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        T2: "TEMP"
EOF
fi
if [[ "${species}" == "PRECIP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        total_precip_timestep: "PRECIP"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt
      color: 'blue'
      marker: 'o'
      linestyle: '--'
      markersize: 8
      linewidth: ${mdl_lw}
EOF
fi
#RAP Chem 
if [[ ${do_rapchem2} -eq 1 && ${has_rapchem2} -eq 1 ]]; then
echo "Including RAP-Chem 06Z"
cat << EOF >> control.yaml.${todays_date}
  RAP-Chem_06Z: # model label
    files: ${rapchem2_datadir}/*
    mod_type: 'wrfchem'
    mod_kwargs:
      surf_only_nc: True
    radius_of_influence: 12000 #meters
    mapping: #model species name : obs species name
EOF
if [[ "${species}" == "AOD_550" ]]; then
cat << EOF >> control.yaml.${todays_date}
      aeronet:
EOF
else
cat << EOF >> control.yaml.${todays_date}
      airnow:
EOF
fi
if [[ "${species}" == "CO" ]]; then
        cat << EOF  >> control.yaml.${todays_date}
        co: CO
EOF
fi
if [[ "${species}" == "NO2" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        no2: NO2
EOF
fi
if [[ "${species}" == "PM25" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM2_5_DRY: "PM2.5"
EOF
fi
if [[ "${species}" == "PM10" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM10: "PM10"
EOF
fi
if [[ "${species}" == "OZONE" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        o3: "OZONE"
EOF
fi
if [[ "${species}" == "AOD_550" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        AOD550: "aod_550nm"
EOF
fi
if [[ "${species}" == "TEMP" ]]; then
	cat << EOF >> control.yaml.${todays_date}
        T2: "TEMP"
EOF
fi
if [[ "${species}" == "PRECIP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        total_precip_timestep: "PRECIP"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt
      color: 'dodgerblue'
      marker: '^'
      linestyle: '--'
      markersize: 8
      linewidth: ${mdl_lw}
EOF
fi
#RAP Chem
if [[ ${do_rapchem3} -eq 1 && ${has_rapchem3} -eq 1 ]]; then
echo "Including RAP-Chem 06Z (N-1)"
cat << EOF >> control.yaml.${todays_date}
  RAP-Chem_06Z_N-1: # model label
    files: ${rapchem3_datadir}/*
    mod_type: 'wrfchem'
    mod_kwargs:
      surf_only_nc: True
    radius_of_influence: 12000 #meters
    mapping: #model species name : obs species name
EOF
if [[ "${species}" == "AOD_550" ]]; then
cat << EOF >> control.yaml.${todays_date}
      aeronet:
EOF
else
cat << EOF >> control.yaml.${todays_date}
      airnow:
EOF
fi
if [[ "${species}" == "CO" ]]; then
        cat << EOF  >> control.yaml.${todays_date}
        co: CO
EOF
fi
if [[ "${species}" == "NO2" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        no2: NO2
EOF
fi
if [[ "${species}" == "PM25" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM2_5_DRY: "PM2.5"
EOF
fi
if [[ "${species}" == "PM10" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM10: "PM10"
EOF
fi
if [[ "${species}" == "OZONE" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        o3: "OZONE"
EOF
fi
if [[ "${species}" == "AOD_550" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        AOD550: "aod_550nm"
EOF
fi
if [[ "${species}" == "TEMP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        T2: "TEMP"
EOF
fi
if [[ "${species}" == "PRECIP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        total_precip_timestep: "PRECIP"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt
      color: 'slateblue'
      marker: '^'
      linestyle: '--'
      markersize: 8
      linewidth: ${mdl_lw}
EOF
fi
#WRFChem-firexaq
if [[ ${do_wrfchem_firexaq} -eq 1 && ${has_wrfchem_firexaq} -eq 1 ]]; then
echo "Including NCAR's FIREX-AQ"
cat << EOF >> control.yaml.${todays_date}
  NCAR_FIREX-AQ: # model label
    files: ${wrfchem_firexaq_datadir}/wrfout*
    mod_type: 'wrfchem'
    mod_kwargs:
      surf_only_nc: True
    radius_of_influence: 12000 #meters
    mapping: #model species name : obs species name
EOF
if [[ "${species}" == "AOD_550" ]]; then
cat << EOF >> control.yaml.${todays_date}
      aeronet:
EOF
else
cat << EOF >> control.yaml.${todays_date}
      airnow:
EOF
fi
if [[ "${species}" == "CO" ]]; then
        cat << EOF  >> control.yaml.${todays_date}
        co: 'CO'
EOF
fi
if [[ "${species}" == "NO2" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        no2: NO2
EOF
fi
if [[ "${species}" == "PM25" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM2_5_DRY: "PM2.5"
EOF
fi
if [[ "${species}" == "PM10" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM10: "PM10"
EOF
fi
if [[ "${species}" == "OZONE" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        o3: "OZONE"
EOF
fi
if [[ "${species}" == "AOD_550" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        AOD_550: "aod_550nm"
EOF
fi
if [[ "${species}" == "TEMP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        T2: "TEMP"
EOF
fi
if [[ "${species}" == "PRECIP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PRECIP: "PRECIP"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt
      color: 'darkgreen'
      marker: 's'
      linestyle: '--'
      markersize: 8
      linewidth: ${mdl_lw}
    vert: False
EOF
fi
#----------------------------------------------------------------------------
if [[ ${do_wrfchem_firexaq2} -eq 1 && ${has_wrfchem_firexaq2} -eq 1 ]]; then
echo "Including NCAR's FIREX-AQ (N-1)"
cat << EOF >> control.yaml.${todays_date}
  NCAR_FIREX-AQ_N-1: # model label
    files: ${wrfchem_firexaq2_datadir}/wrfout*
    mod_type: 'wrfchem'
    mod_kwargs:
      surf_only_nc: True
    radius_of_influence: 12000 #meters
    mapping: #model species name : obs species name
EOF
if [[ "${species}" == "AOD_550" ]]; then
cat << EOF >> control.yaml.${todays_date}
      aeronet:
EOF
else
cat << EOF >> control.yaml.${todays_date}
      airnow:
EOF
fi
if [[ "${species}" == "CO" ]]; then
        cat << EOF  >> control.yaml.${todays_date}
        co: 'CO'
EOF
fi
if [[ "${species}" == "NO2" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        no2: NO2
EOF
fi
if [[ "${species}" == "PM25" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM2_5_DRY: "PM2.5"
EOF
fi
if [[ "${species}" == "PM10" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM10: "PM10"
EOF
fi
if [[ "${species}" == "OZONE" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        o3: "OZONE"
EOF
fi
if [[ "${species}" == "AOD_550" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        AOD_550: "aod_550nm"
EOF
fi
if [[ "${species}" == "TEMP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        T2: "TEMP"
EOF
fi
if [[ "${species}" == "PRECIP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PRECIP: "PRECIP"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt
      color: 'springgreen'
      marker: 's'
      linestyle: '--'
      markersize: 8
      linewidth: ${mdl_lw}
    vert: False
EOF
fi
#WRFChem-aqwatch
if [[ ${do_wrfchem_aqwatch} -eq 1 && ${has_wrfchem_aqwatch} -eq 1 ]]; then
echo "Including NCAR's AQ-WATCH"
cat << EOF >> control.yaml.${todays_date}
  NCAR_AQ-WATCH: # model label
    files: ${wrfchem_aqwatch_datadir}/wrfout*
    mod_type: 'wrfchem'
    mod_kwargs:
      surf_only_nc: True
    radius_of_influence: 12000 #meters
    mapping: #model species name : obs species name
EOF
if [[ "${species}" == "AOD_550" ]]; then
cat << EOF >> control.yaml.${todays_date}
      aeronet:
EOF
else
cat << EOF >> control.yaml.${todays_date}
      airnow:
EOF
fi
if [[ "${species}" == "CO" ]]; then
        cat << EOF  >> control.yaml.${todays_date}
        co: 'CO'
EOF
fi
if [[ "${species}" == "NO2" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        no2: NO2
EOF
fi
if [[ "${species}" == "PM25" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM2_5_DRY: "PM2.5"
EOF
fi
if [[ "${species}" == "PM10" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM10: "PM10"
EOF
fi
if [[ "${species}" == "OZONE" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        o3: "OZONE"
EOF
fi
if [[ "${species}" == "AOD_550" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        AOD_550: "aod_550nm"
EOF
fi
if [[ "${species}" == "TEMP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        T2: "TEMP"
EOF
fi
if [[ "${species}" == "PRECIP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PRECIP: "PRECIP"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt
      color: 'limegreen'
      marker: '^'
      linestyle: '--'
      markersize: 8
      linewidth: ${mdl_lw}
    vert: False
EOF
fi
#WRFChem-aqwatch
if [[ ${do_wrfchem_aqwatch2} -eq 1 && ${has_wrfchem_aqwatch2} -eq 1 ]]; then
echo "Including NCAR's AQ-WATCH (N-1)"
cat << EOF >> control.yaml.${todays_date}
  NCAR_AQ-WATCH_N-1: # model label
    files: ${wrfchem_aqwatch2_datadir}/wrfout*
    mod_type: 'wrfchem'
    mod_kwargs:
      surf_only_nc: True
    radius_of_influence: 12000 #meters
    mapping: #model species name : obs species name
EOF
if [[ "${species}" == "AOD_550" ]]; then
cat << EOF >> control.yaml.${todays_date}
      aeronet:
EOF
else
cat << EOF >> control.yaml.${todays_date}
      airnow:
EOF
fi
if [[ "${species}" == "CO" ]]; then
        cat << EOF  >> control.yaml.${todays_date}
        co: 'CO'
EOF
fi
if [[ "${species}" == "NO2" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        no2: NO2
EOF
fi
if [[ "${species}" == "PM25" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM2_5_DRY: "PM2.5"
EOF
fi
if [[ "${species}" == "PM10" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM10: "PM10"
EOF
fi
if [[ "${species}" == "OZONE" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        o3: "OZONE"
EOF
fi
if [[ "${species}" == "AOD_550" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        AOD_550: "aod_550nm"
EOF
fi
if [[ "${species}" == "TEMP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        T2: "TEMP"
EOF
fi
if [[ "${species}" == "PRECIP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PRECIP: "PRECIP"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt
      color: 'forestgreen'
      marker: '^'
      linestyle: '--'
      markersize: 8
      linewidth: ${mdl_lw}
    vert: False
EOF
fi
#CMAQ Oper
if [[ ${do_cmaq_oper12} -eq 1 && ${has_cmaq_oper12} -eq 1 ]]; then
echo "Including NAQFC_CMAQ_oper_12Z"
cat << EOF >> control.yaml.${todays_date}
  NAQFC_CMAQ_oper_12Z: # model label
    files: ${cmaq_oper12_datadir}/*
    mod_type: 'cmaq'
    radius_of_influence: 12000 #meters
    #variables: #Opt
    mapping: #model species name : obs species name
      airnow:
EOF
if [[ "${species}" == "CO" ]]; then
	cat << EOF >> control.yaml.${todays_date}
        CO: CO
EOF
fi
if [[ "${species}" == "NO2" ]]; then
	cat << EOF >> control.yaml.${todays_date}
        NO2: NO2
EOF
fi
if [[ "${species}" == "PM25" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM25_TOT: "PM2.5"
EOF
fi
if [[ "${species}" == "PM10" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM10: "PM10"
EOF
fi
if [[ "${species}" == "OZONE" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        O3: "OZONE"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt 
      color: 'red'
      marker: 'o'
      linestyle: '--'
      markersize: 8
      linewidth: ${mdl_lw}
EOF
fi
#CMAQ Oper
if [[ ${do_cmaq_oper12_2} -eq 1 && ${has_cmaq_oper12_2} -eq 1 ]]; then
echo "Including NAQFC_CMAQ_oper_12Z"
cat << EOF >> control.yaml.${todays_date}
  NAQFC_CMAQ_oper_12Z_N-2: # model label
    files: ${cmaq_oper12_2_datadir}/*
    mod_type: 'cmaq'
    radius_of_influence: 12000 #meters
    #variables: #Opt
    #  CO:
    #    unit_scale: 1000.0
    #    unit_scale_method: '*'
    mapping: #model species name : obs species name
      airnow:
EOF
if [[ "${species}" == "CO" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        CO: CO
EOF
fi
if [[ "${species}" == "NO2" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        NO2: NO2
EOF
fi
if [[ "${species}" == "PM25" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM25_TOT: "PM2.5"
EOF
fi
if [[ "${species}" == "PM10" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM10: "PM10"
EOF
fi
if [[ "${species}" == "OZONE" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        O3: "OZONE"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt 
      color: 'coral'
      marker: 'o'
      linestyle: '--'
      markersize: 8
      linewidth: ${mdl_lw}
EOF
fi
#CMAQ expr
if [[ ${do_cmaq_expr12} -eq 1 && ${has_cmaq_expr12} -eq 1 ]]; then
echo "Including NAQFC_CMAQ_expr_12Z"
cat << EOF >> control.yaml.${todays_date}
  NAQFC_CMAQ_expr_12Z: # model label
    files: ${cmaq_expr12_datadir}/*
    mod_type: 'cmaq'
    radius_of_influence: 12000 #meters
    mapping: #model species name : obs species name
      airnow:
EOF
if [[ "${species}" == "CO" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        CO: CO
EOF
fi
if [[ "${species}" == "NO2" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        NO2: NO2
EOF
fi
if [[ "${species}" == "PM25" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM25_TOT: "PM2.5"
EOF
fi
if [[ "${species}" == "PM10" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM10: "PM10"
EOF
fi
if [[ "${species}" == "OZONE" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        O3: "OZONE"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt
      color: 'magenta'
      marker: 'x'
      linestyle: '--'
      markersize: 8
      linewidth: ${mdl_lw}
EOF
fi
#CMAQ Oper
if [[ ${do_cmaq_oper06} -eq 1 && ${has_cmaq_oper06} -eq 1 ]]; then
echo "Including NAQFC_CMAQ_oper_06Z"
cat << EOF >> control.yaml.${todays_date}
  NAQFC_CMAQ_oper_06Z: # model label
    files: ${cmaq_oper06_datadir}/*
    mod_type: 'cmaq'
    radius_of_influence: 12000 #meters
    mapping: #model species name : obs species name
      airnow:
EOF
if [[ "${species}" == "CO" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        CO: CO
EOF
fi
if [[ "${species}" == "NO2" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        NO2: NO2
EOF
fi
if [[ "${species}" == "PM25" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM25_TOT: "PM2.5"
EOF
fi
if [[ "${species}" == "PM10" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM10: "PM10"
EOF
fi
if [[ "${species}" == "OZONE" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        O3: "OZONE"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt 
      color: 'darkred'
      marker: '^'
      linestyle: '--'
      markersize: 8
      linewidth: ${mdl_lw}
EOF
fi

#CMAQ expr
if [[ ${do_cmaq_expr06} -eq 1 && ${has_cmaq_expr06} -eq 1 ]]; then
echo "Including NAQFC_CMAQ_expr_06Z"
cat << EOF >> control.yaml.${todays_date}
  NAQFC_CMAQ_expr_06Z: # model label
    files: ${cmaq_expr06_datadir}/*
    mod_type: 'cmaq'
    radius_of_influence: 12000 #meters
    mapping: #model species name : obs species name
      airnow:
EOF
if [[ "${species}" == "CO" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        CO: CO
EOF
fi
if [[ "${species}" == "NO2" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        NO2: NO2
EOF
fi
if [[ "${species}" == "PM25" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM25_TOT: "PM2.5"
EOF
fi
if [[ "${species}" == "PM10" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM10: "PM10"
EOF
fi
if [[ "${species}" == "OZONE" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        O3: "OZONE"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt
      color: 'violet'
      marker: '^'
      linestyle: '--'
      markersize: 8
      linewidth: ${mdl_lw}
EOF
fi
#HRRR Smoke
if [[ ${do_hrrrsmoke06} -eq 1 && ${has_hrrrsmoke06} -eq 1 ]]; then
echo "Including HRRR_Smoke_06Z"
cat << EOF >> control.yaml.${todays_date}
  HRRR-Smoke_06Z: # model label
    files: ${hrrr_smoke06_datadir}/*
    mod_type: 'hrrr'
    mod_kwargs:
      surf_only_nc: True
    radius_of_influence: 12000 #meters
    mapping: #model species name : obs species name
EOF
if [[ "${species}" == "AOD_550" ]]; then
cat << EOF >> control.yaml.${todays_date}
      aeronet:
EOF
else
cat << EOF >> control.yaml.${todays_date}
      airnow:
EOF
fi
if [[ "${species}" == "AOD_550" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        AOD_550: "aod_550nm"
EOF
fi
if [[ "${species}" == "PM25" ]]; then
	cat << EOF >> control.yaml.${todays_date}
        PM2_5_DRY: "PM2.5"
EOF
fi
if [[ "${species}" == "PRECIP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        APCP_surface: "PRECIP"
EOF
fi
if [[ "${species}" == "TEMP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        TMP_2maboveground: "TEMP"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt
      color: 'orange'
      marker: '^'
      markersize: 8
      linestyle: '--'
      linewidth: ${mdl_lw}
EOF
fi

if [[ ${do_hrrrsmoke18} -eq 1 && ${has_hrrrsmoke18} -eq 1 ]]; then 
echo "Including HRRR_Smoke_18Z"
cat << EOF >> control.yaml.${todays_date}
  HRRR-Smoke_18Z: # model label
    files: ${hrrr_smoke18_datadir}/*
    mod_type: 'hrrr'
    mod_kwargs:
      surf_only_nc: True
    radius_of_influence: 12000 #meters
    mapping: #model species name : obs species name
EOF
if [[ "${species}" == "AOD_550" ]]; then
cat << EOF >> control.yaml.${todays_date}
      aeronet:
EOF
else
cat << EOF >> control.yaml.${todays_date}
      airnow:
EOF
fi
if [[ "${species}" == "AOD_550" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        AOD_550: "aod_550nm"
EOF
fi
if [[ "${species}" == "PM25" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM2_5_DRY: "PM2.5"
EOF
fi
if [[ "${species}" == "PRECIP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        APCP_surface: "PRECIP"
EOF
fi
if [[ "${species}" == "TEMP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        TMP_2maboveground: "TEMP"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt
      color: 'darkorange'
      marker: '^'
      markersize: 8
      linestyle: '--'
      linewidth: ${mdl_lw}
EOF
fi

if [[ ${do_geos_cf12} -eq 1 && ${has_geos_cf12} -eq 1 ]]; then
echo "Including GEOS-CF_12Z"
cat << EOF >> control.yaml.${todays_date}
  GEOS-CF_12Z: # model label
    files: ${geos_cf12_datadir}/geos-cf.*.nc
#    files: ${geos_cf12_datadir}/test/test.nc
    mod_type: 'hrrr'
    mod_kwargs:
      surf_only_nc: True
    radius_of_influence: 12000 #meters
    mapping: #model species name : obs species name
EOF
if [[ "${species}" == "AOD_550" ]]; then
cat << EOF >> control.yaml.${todays_date}
      aeronet:
EOF
else
cat << EOF >> control.yaml.${todays_date}
      airnow:
EOF
fi
if [[ "${species}" == "CO" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        CO: CO
EOF
fi
if [[ "${species}" == "NO2" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        NO2: NO2
EOF
fi
if [[ "${species}" == "PM25" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM25_RH35_GCC: "PM2.5"
EOF
fi
if [[ "${species}" == "OZONE" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        O3: "OZONE"
EOF
fi
if [[ "${species}" == "AOD_550" ]]; then
	cat << EOF >> control.yaml.${todays_date}
        AOD_550: "aod_550nm"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt
      color: 'cyan'
      marker: '^'
      linestyle: '--'
      markersize: 8
      linewidth: ${mdl_lw}
EOF
fi

#onlinecmaq_b1
if [[ ${do_onlinecmaq_b1} -eq 1 && ${has_onlinecmaq_b1} -eq 1 ]]; then
echo "Including onlinecmaq_b1"
cat << EOF >> control.yaml.${todays_date}
  online-CMAQ(UFS,v7.0b1): # model label
    files: ${onlinecmaq_b1_datadir}/aqm.t12z.chem_3d*
    mod_type: 'rrfs'
    radius_of_influence: 12000 #meters
    mapping: #model species name : obs species name
      airnow:
EOF
if [[ "${species}" == "CO" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        co: CO
EOF
fi
if [[ "${species}" == "NO2" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        no2: NO2
EOF
fi
if [[ "${species}" == "PM25" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM25_TOT: "PM2.5"
EOF
fi
if [[ "${species}" == "OZONE" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        o3: "OZONE"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt
      color: 'blueviolet'
      marker: 'x'
      linestyle: '--'
      markersize: 8
      linewidth: ${mdl_lw}
EOF
fi

#onlinecmaq_c3
if [[ ${do_onlinecmaq_c3} -eq 1 && ${has_onlinecmaq_c3} -eq 1 ]]; then
echo "Including onlinecmaq_c3"
cat << EOF >> control.yaml.${todays_date}
  online-CMAQ(UFS,v7.0c3): # model label
    files: ${onlinecmaq_c3_datadir}/aqm.t12z.chem_3d*
    mod_type: 'rrfs'
    radius_of_influence: 12000 #meters
    mapping: #model species name : obs species name
      airnow:
EOF
if [[ "${species}" == "CO" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        co: CO
EOF
fi
if [[ "${species}" == "NO2" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        no2: NO2
EOF
fi
if [[ "${species}" == "PM25" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM25_TOT: "PM2.5"
EOF
fi
if [[ "${species}" == "OZONE" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        o3: "OZONE"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt
      color: 'mediumorchid'
      marker: 'o'
      linestyle: '--'
      markersize: 8
      linewidth: ${mdl_lw}
EOF
fi

#RAQMS
if [[ ${do_raqms} -eq 1 && ${has_raqms} -eq 1 ]]; then
echo "Including RAQMS"
cat << EOF >> control.yaml.${todays_date}
  RAQMS: # model label
    files: ${raqms_datadir}/*
    mod_type: 'raqms'
    radius_of_influence: 12000 #meters
    variables: #Opt
      ico:
        unit_scale: 1000000000.0
        unit_scale_method: '*'
      ino2:
        unit_scale: 1000000000.0
        unit_scale_method: '*'
      PM25_RH35_GCC:
        unit_scale: 1000000000.0
        unit_scale_method: '*'
      o3vmr:
        unit_scale: 1000000000.0
        unit_scale_method: '*'
    mapping: #model species name : obs species name
EOF
if [[ "${species}" == "AOD_550" ]]; then
cat << EOF >> control.yaml.${todays_date}
      aeronet:
EOF
else
cat << EOF >> control.yaml.${todays_date}
      airnow:
EOF
fi
if [[ "${species}" == "CO" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        ico: CO
EOF
fi
if [[ "${species}" == "NO2" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        ino2: NO2
EOF
fi
if [[ "${species}" == "PM25" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM25_RH35_GCC: "PM2.5"
EOF
fi
if [[ "${species}" == "OZONE" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        o3vmr: "OZONE"
EOF
fi
if [[ "${species}" == "AOD_550" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        aod: "aod_550nm"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt 
      color: 'hotpink'
      marker: 'o'
      linestyle: '--'
      markersize: 8
      linewidth: ${mdl_lw}
EOF
fi

#RAP Smoke
if [[ ${do_rapsmoke} -eq 1 && ${has_rapsmoke} -eq 1 ]]; then
echo "Including RAP_Smoke"
cat << EOF >> control.yaml.${todays_date}
  RAP-Smoke: # model label
    files: ${rapsmoke_datadir}/*
    mod_type: 'hrrr'
    mod_kwargs:
      surf_only_nc: True
    radius_of_influence: 12000 #meters
    mapping: #model species name : obs species name
      airnow:
EOF
if [[ "${species}" == "PM25" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM2_5_DRY: "PM2.5"
EOF
fi
if [[ "${species}" == "PRECIP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PRECIP: "PRECIP"
EOF
fi
if [[ "${species}" == "TEMP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        T2: "TEMP"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt
      color: 'peru'
      marker: '^'
      markersize: 8
      linestyle: '--'
      linewidth: ${mdl_lw}
EOF
fi

#HRRR Chem
if [[ ${do_hrrrchem} -eq 1 && ${has_hrrrchem} -eq 1 ]]; then
echo "Including HRRR-Chem 00Z"
cat << EOF >> control.yaml.${todays_date}
  HRRR-Chem_00Z: # model label
    files: ${hrrrchem_datadir}/*
    mod_type: 'wrfchem'
    mod_kwargs:
      surf_only_nc: True
    radius_of_influence: 12000 #meters
    mapping: #model species name : obs species name
EOF
if [[ "${species}" == "AOD_550" ]]; then
cat << EOF >> control.yaml.${todays_date}
      aeronet:
EOF
else
cat << EOF >> control.yaml.${todays_date}
      airnow:
EOF
fi
if [[ "${species}" == "CO" ]]; then
        cat << EOF  >> control.yaml.${todays_date}
        co: CO
EOF
fi
if [[ "${species}" == "NO2" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        no2: NO2
EOF
fi
if [[ "${species}" == "PM25" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM2_5_DRY: "PM2.5"
EOF
fi
if [[ "${species}" == "PM10" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        PM10: "PM10"
EOF
fi
if [[ "${species}" == "OZONE" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        o3: "OZONE"
EOF
fi
if [[ "${species}" == "AOD_550" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        AOD550: "aod_550nm"
EOF
fi
if [[ "${species}" == "TEMP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        T2: "TEMP"
EOF
fi
if [[ "${species}" == "PRECIP" ]]; then
        cat << EOF >> control.yaml.${todays_date}
        total_precip_timestep: "PRECIP"
EOF
fi
cat << EOF >> control.yaml.${todays_date}
    projection: None
    plot_kwargs: #Opt
      color: 'teal'
      marker: 'o'
      linestyle: '--'
      markersize: 8
      linewidth: ${mdl_lw}
EOF
fi


#Observation Types 
if [[ ${species} == "AOD_550" ]]; then
cat << EOF >> control.yaml.${todays_date}

obs:
  aeronet: # obs label
    filename: test5.nc
    obs_type: pt_sfc
    variables: #Opt 
EOF
else
cat << EOF >> control.yaml.${todays_date}

obs:
  airnow: # obs label
    use_airnow: True
    filename: test5.nc
    obs_type: pt_sfc
    variables: #Opt
EOF
fi
if [[ ${species} == "AOD_550" ]]; then
echo "Including AOD_550"
cat << EOF >> control.yaml.${todays_date}
      aod_550nm:
        unit_scale: 1
        unit_scale_method: '*' # Multiply = '*' , Add = '+', subtract = '-', divide = '/'
        nan_value: -1.0 # Set this value to NaN
        ylabel_plot: 'Aeronet 550nm AOD' #Optional to set ylabel so can include units and/or instr etc.
        vmin_plot: 0.0 #Opt Min for y-axis during plotting. To apply to a plot, change restrict_yaxis = True.
        vmax_plot: 1.0 #Opt Max for y-axis during plotting. To apply to a plot, change restrict_yaxis = True.
        vdiff_plot: 0.2 #Opt +/- range to use in bias plots. To apply to a plot, change restrict_yaxis = True.
        nlevels_plot: 23 #Opt number of levels used in colorbar for contourf plot.
EOF
fi
if [[ "${species}" == "PM25" ]]; then
echo "Including PM25 in the Analysis"

cat << EOF >> control.yaml.${todays_date}

      PM2.5:
        unit_scale: 1
        unit_scale_method: '*' # Multiply = '*' , Add = '+', subtract = '-', divide = '/'
        nan_value: -1.0 # Set this value to NaN
        #The obs_min, obs_max, and nan_values are set to NaN first and then the unit conversion is applied.
        ylabel_plot: 'PM2.5 (ug/m3)' #Optional to set ylabel so can include units and/or instr etc.
        ty_scale: 2.0 #Opt
        vmin_plot: 0.0 #Opt Min for y-axis during plotting. To apply to a plot, change restrict_yaxis = True.
        vmax_plot: 22.0 #Opt Max for y-axis during plotting. To apply to a plot, change restrict_yaxis = True.
        vdiff_plot: 15.0 #Opt +/- range to use in bias plots. To apply to a plot, change restrict_yaxis = True.
        nlevels_plot: 23 #Opt number of levels used in colorbar for contourf plot.

EOF
fi
if [[ "${species}" == "PM10" ]]; then
echo "Including PM10 in the Analysis"

cat << EOF >> control.yaml.${todays_date}

      PM10:
        unit_scale: 1
        unit_scale_method: '*' # Multiply = '*' , Add = '+', subtract = '-', divide = '/'
        #obs_min: 0 # set all values less than this value to NaN
        #obs_max: 100 # set all values greater than this value to NaN
        nan_value: -1.0 # Set this value to NaN
        #The obs_min, obs_max, and nan_values are set to NaN first and then the unit conversion is applied.
        ylabel_plot: 'PM10 (ug/m3)' #Optional to set ylabel so can include units and/or instr etc.
        ty_scale: 2.0 #Opt
        vmin_plot: 0.0 #Opt Min for y-axis during plotting. To apply to a plot, change restrict_yaxis = True.
        vmax_plot: 22.0 #Opt Max for y-axis during plotting. To apply to a plot, change restrict_yaxis = True.
        vdiff_plot: 15.0 #Opt +/- range to use in bias plots. To apply to a plot, change restrict_yaxis = True.
        nlevels_plot: 23 #Opt number of levels used in colorbar for contourf plot.

EOF
fi
if [[ "${species}" == "OZONE" ]]; then
echo "including OZONE"
cat << EOF >> control.yaml.${todays_date}

      OZONE:
        unit_scale: 1 #Opt Scaling factor 
        unit_scale_method: '*' #Opt Multiply = '*' , Add = '+', subtract = '-', divide = '/'
        nan_value: -1.0 # Opt Set this value to NaN
        ylabel_plot: 'Ozone (ppbv)'
        vmin_plot: 15.0 #Opt Min for y-axis during plotting. To apply to a plot, change restrict_yaxis = True.
        vmax_plot: 55.0 #Opt Max for y-axis during plotting. To apply to a plot, change restrict_yaxis = True.
        vdiff_plot: 20.0 #Opt +/- range to use in bias plots. To apply to a plot, change restrict_yaxis = True.
        nlevels_plot: 21 #Opt number of levels used in colorbar for contourf plot.


EOF
fi

if [[ "${species}" == "WS" ]]; then
echo "WS"
cat << EOF >> control.yaml.${todays_date}

      WS:
        unit_scale: 0.514  # convert obs knots-->m/s
        unit_scale_method: '*'
        obs_min: 0.2 # m/s

EOF
fi

if [[ "${species}" == "PRSFC" ]]; then
echo "PRSFC"
cat << EOF >> control.yaml.${todays_date}

      PRSFC:
        unit_scale: 0.01  # convert model Pascals-->millibars
        unit_scale_method: '*'

EOF
fi

if [[ "${species}" == "PRECIP" ]]; then
echo "PRECIP"
cat << EOF >> control.yaml.${todays_date}

      PRECIP:
        unit_scale: 0.1  # convert obs mm-->cm
        unit_scale_method: '*'

EOF
fi

if [[ "${species}" == "TEMP" ]]; then
echo "TEMP"
cat << EOF >> control.yaml.${todays_date}

      TEMP:
        unit_scale: 273.16
        unit_scale_method: '+'
        nan_value: -1.0
        ylabel_plot: '2-m Temperature (K)'

EOF
fi

if [[ "${species}" == "CO" ]]; then
echo "CO"
cat << EOF >> control.yaml.${todays_date}

      CO:
        unit_scale: 1000. #Convert from ppmv to ppbv.
        unit_scale_method: '*' # Multiply = '*' , Add = '+', subtract = '-', divide = '/'
        nan_value: -1.0 # Set this value to NaN
        #The obs_min, obs_max, and nan_values are set to NaN first and then the unit conversion is applied.
        ylabel_plot: 'CO (ppbv)' #Optional to set ylabel so can include units and/or instr etc.
        vmin_plot: 50.0 #Opt Min for y-axis during plotting. To apply to a plot, change restrict_yaxis = True.
        vmax_plot: 750.0 #Opt Max for y-axis during plotting. To apply to a plot, change restrict_yaxis = True.
        vdiff_plot: 400.0 #Opt +/- range to use in bias plots. To apply to a plot, change restrict_yaxis = True
        nlevels_plot: 15 #Opt number of levels used in colorbar for contourf plot.

EOF
fi


if [[ "${species}" == "SO2" ]]; then
echo "SO2"
cat << EOF >> control.yaml.${todays_date}

      SO2:
        nan_value: -1.0 # Set this value to NaN
        ylabel_plot: 'SO2 (ppbv)' #Optional to set ylabel so can include units and/or instr etc.

EOF
fi

if [[ "${species}" == "NO" ]]; then
echo "NO"
cat << EOF >> control.yaml.${todays_date}

      'NO':
        nan_value: -1.0 # Set this value to NaN
        ylabel_plot: 'NO (ppbv)' #Optional to set ylabel so can include units and/or instr etc.
        vmin_plot: 0.0 #Opt Min for y-axis during plotting. To apply to a plot, change restrict_yaxis = True.
        vmax_plot: 20.0 #Opt Max for y-axis during plotting. To apply to a plot, change restrict_yaxis = True.
        vdiff_plot: 15.0 #Opt +/- range to use in bias plots. To apply to a plot, change restrict_yaxis = True.
        nlevels_plot: 21 #Opt number of levels used in colorbar for contourf plot.

EOF
fi

if [[ "${species}" == "NO2" ]]; then
echo "NO2"
cat << EOF >> control.yaml.${todays_date}

      NO2:
        #obs_max: 1 # ppbv
        nan_value: -1.0 # Set this value to NaN
        ylabel_plot: 'NO2 (ppbv)' #Optional to set ylabel so can include units and/or instr etc.
        vmin_plot: 0.0 #Opt Min for y-axis during plotting. To apply to a plot, change restrict_yaxis = True.
        vmax_plot: 20.0 #Opt Max for y-axis during plotting. To apply to a plot, change restrict_yaxis = True.
        vdiff_plot: 15.0 #Opt +/- range to use in bias plots. To apply to a plot, change restrict_yaxis = True.
        nlevels_plot: 21 #Opt number of levels used in colorbar for contourf plot.

EOF
fi


#Plot Types 

cat << EOF >> control.yaml.${todays_date}
plots:
EOF
# Now insert which plot groups  
#Timeseries plot group 1 
if [[ ${timeseries} -eq 1 ]]; then
echo "Now doing timeseries, plot group 1" 
cat << EOF >> control.yaml.${todays_date} 
  plot_grp1:
    type: 'timeseries' # plot type
    fig_kwargs: #Opt to define figure options
      figsize: [12,6] # figure size if multiple plots
    default_plot_kwargs: # Opt to define defaults for all plots. Model kwargs overwrite these.
      linewidth: 3.0
      markersize: 10.
    text_kwargs: #Opt
      fontsize: 18.
    domain_type: [${rgn_type[*]}] #List of domain types: 'all' or any domain in obs file. (e.g., airnow: epa_region, state_name, siteid, etc.)
    domain_name: [${rgn_list[*]}] #List of domain names. If domain_type = all domain_name is used in plot title.
    data: [${obs_mdl_list[*]}]
    data_proc:
      rem_obs_nan: True # True: Remove all points where model or obs variable is NaN. False: Remove only points where model variable is NaN.
      ts_select_time: ${ts_select_time} #Time used for avg and plotting: Options: 'time' for UTC or 'time_local'
      ts_avg_window: 'H' # Options: None for no averaging or list pandas resample rule (e.g., 'H', 'D')
      set_axis: False #If select True, add vmin_plot and vmax_plot for each variable in obs.

EOF
fi 

#Taylor Plot Group 2 
if [[ ${taylor} -eq 1 ]]; then 
echo "Now doing Taylor plot, plot group 2"
cat << EOF >> control.yaml.${todays_date}

  plot_grp2:
    type: 'taylor' # plot type
    fig_kwargs: #Opt to define figure options
      figsize: [8,8] # figure size if multiple plots
    default_plot_kwargs: # Opt to define defaults for all plots. Model kwargs overwrite these.
      linewidth: 2.0
      markersize: 10.
    text_kwargs: #Opt
      fontsize: 16.
    domain_type: [${rgn_type[*]}] #List of domain types: 'all' or any domain in obs file. (e.g., airnow: epa_region, state_name, siteid, etc.)
    domain_name: [${rgn_list[*]}] #List of domain names. If domain_type = all domain_name is used in plot title.
    data: [${obs_mdl_list[*]}]
    data_proc:
      rem_obs_nan: True # True: Remove all points where model or obs variable is NaN. False: Remove only points where model variable is NaN.
      set_axis: True #If select True, add ty_scale for each variable in obs.

EOF
fi

#Spatial Bias Plot Group 3 
if [[ ${spatial_bias} -eq 1 ]]; then
echo "Now doing spatial bias, plot group 3"
cat << EOF >> control.yaml.${todays_date}

  plot_grp3:
    type: 'spatial_bias' # plot type
    fig_kwargs: #For all spatial plots, specify map_kwargs here too.
      states: True
      figsize: [10, 5] # figure size 
    text_kwargs: #Opt
      fontsize: 16.
    domain_type: [${rgn_type[*]}] #List of domain types: 'all' or any domain in obs file. (e.g., airnow: epa_region, state_name, siteid, etc.) 
    domain_name: [${rgn_list[*]}] #List of domain names. If domain_type = all domain_name is used in plot title.
    data: [${obs_mdl_list[*]}]
    data_proc:
      rem_obs_nan: True # True: Remove all points where model or obs variable is NaN. False: Remove only points where model variable is NaN.
      set_axis: True #If select True, add vdiff_plot for each variable in obs.

EOF
fi

#Spatial Overlay 
if [[ ${spatial_overlay} -eq 1 ]]; then 
echo "Now doing spatial overlay, plot group 4"
cat << EOF >> control.yaml.${todays_date}

  plot_grp4:
    type: 'spatial_overlay' # plot type
    fig_kwargs: #For all spatial plots, specify map_kwargs here too.
      states: True
      figsize: [10, 5] # figure size
    text_kwargs: #Opt
      fontsize: 16.
    domain_type: [${rgn_type[*]}] #List of domain types: 'all' or any domain in obs file. (e.g., airnow: epa_region, state_name, siteid, etc.)
    domain_name: [${rgn_list[*]}] #List of domain names. If domain_type = all domain_name is used in plot title.
    data: [${obs_mdl_list[*]}]
    data_proc:
      rem_obs_nan: True # True: Remove all points where model or obs variable is NaN. False: Remove only points where model variable is NaN.
      set_axis: True #If select True, add vmin_plot and vmax_plot for each variable in obs.

EOF
fi


#Boxplot Plot Group 5 
if [[ ${boxplot} -eq 1 ]]; then
echo "Now doing boxplot, plot group 5"
cat << EOF >> control.yaml.${todays_date}

  plot_grp5:
    type: 'boxplot' # plot type
    fig_kwargs: #Opt to define figure options
      figsize: [8, 6] # figure size 
    text_kwargs: #Opt
      fontsize: 20.
    domain_type: [${rgn_type[*]}] #List of domain types: 'all' or any domain in obs file. (e.g., airnow: epa_region, state_name, siteid, etc.) 
    domain_name: [${rgn_list[*]}] #List of domain names. If domain_type = all domain_name is used in plot title.
    data: [${obs_mdl_list[*]}]
    data_proc:
      rem_obs_nan: True # True: Remove all points where model or obs variable is NaN. False: Remove only points where model variable is NaN.
      set_axis: False #If select True, add vmin_plot and vmax_plot for each variable in obs.
EOF
fi

if [ ${do_stats} -eq 1 ]; then
echo "Adding stats"
cat << EOF >> control.yaml.${todays_date}

stats:
  #Stats require positive numbers, so if you want to calculate temperature use Kelvin!
  #Wind direction has special calculations for AirNow if obs name is 'WD'
  stat_list: ['STDO', 'STDP', 'MdnNB', 'NO','NOP','NP','MO','MP', 'MdnO', 'MdnP', 'RM', 'RMdn', 'MB', 'MdnB', 'NMB', 'NMdnB', 'FB', 'NME', 'R2', 'RMSE', 'IOA', 'AC'] #List stats to calculate. Dictionary of definitions included in plots/proc_stats.py Only stats listed below are currently working.
  #stat_list: ['MB']
  #Full calc list ['STDO', 'STDP', 'MdnNB', 'NO','NOP','NP','MO','MP', 'MdnO', 'MdnP', 'RM', 'RMdn', 'MB', 'MdnB', 'NMB', 'NMdnB', 'FB', 'NME', 'R2', 'RMSE', 'IOA', 'AC']
  round_output: 2 #Opt, defaults to rounding to 3rd decimal place.
  output_table: True #Always outputs a .txt file. Optional to also output as a table.
  output_table_kwargs: #Opt For all spatial plots, specify map_kwargs here too.
    figsize: [8, 11] # figure size 
    fontsize: 12.
    xscale: 1.1
    yscale: 1.1
    edges: 'horizontal'
  domain_type: [${rgn_type[*]}] #List of domain types: 'all' or any domain in obs file. (e.g., airnow: epa_region, state_name, siteid, etc.) 
  domain_name: [${rgn_list[*]}] #List of domain names. If domain_type = all domain_name is used in plot title.
  data: [${obs_mdl_list[*]}] # make this a list of pairs in obs_model where the obs is the obs label and model is the model_label
EOF
fi
# First remove any file named control.yaml
rm -f control.yaml
# Now link to the one we have contructed for today
ln -s control.yaml.${todays_date} control.yaml

#once the test5.nc file is written and linked we can 
#call and run Monet-analysis-example-plots-wrf-rapchemtest
#this will pair the observational data with model data and generate the chosen plots with chosen models and species
if [[ "${species}" == "AOD_550" ]]; then
   python Monet-analysis-example-plots-wrf-rapchemtest_aeronet.py
else
 #  python -m trace -t  Monet-analysis-example-plots-wrf-rapchemtest.py ${do_stats} >> rrfsdebug.txt 2>&1
   python Monet-analysis-example-plots-wrf-rapchemtest.py ${do_stats}
fi

output_directory=/wrk/csd4/rahmadov/RAP-Chem/MELODIES-MONET-analysis/run_sites/plots_out/${todays_date}
mkdir -p ${output_directory}
mv /wrk/csd4/rahmadov/RAP-Chem/MELODIES-MONET-analysis/run_sites/plot_*.png ${output_directory} 
#mv /wrk/csd4/rahmadov/RAP-Chem/MELODIES-MONET-analysis/run_sites/stats*.csv ${output_directory}
#mv /wrk/csd4/rahmadov/RAP-Chem/MELODIES-MONET-analysis/run_sites/stats*.png ${output_directory}

#done # Plot loop
# plots will be outputted to /wrk/csd4/rahmadov/RAP-Chem/MONET-analysis/MONET-analysis/monet_analysis 
# taking outputted plots and moving them to a directory 

#renaming the plots by replacing spaces with underscores
for f in  ${output_directory}/*
do
  new="${f// /_}"
  if [ "$new" != "$f" ]
  then
    if [ -e "$new" ]
    then
      echo not renaming \""$f"\" because \""$new"\" already exists
    else
      echo moving "$f" to "$new"
    mv "$f" "$new"
  fi
fi
done

#create a folium map
if [ ${site} -eq 1 ]; then
	if [[ "${species}" == "AOD_550" ]]; then
       		python folium_site_map_aeronet.py ${species} ${start_time_reformat} ${end_time_reformat} ${todays_date}
	else
		python folium_site_map_airnow.py ${species} ${start_time_reformat} ${end_time_reformat} ${todays_date}
	fi

fi

#move the html and site file to the plot directory
mv /wrk/csd4/rahmadov/RAP-Chem/MELODIES-MONET-analysis/run_sites/sitefile.txt ${output_directory}/sitefile.txt.${species}.${todays_date}
mv /wrk/csd4/rahmadov/RAP-Chem/MELODIES-MONET-analysis/run_sites/${todays_date}_*.html ${output_directory}

unset rgn_list
unset rgn_type
unset mdl_list
unset obs_mdl_list
mv control.yaml.${todays_date} ${output_directory}/control.yaml.${species}.${todays_date}

done # Species Loop



