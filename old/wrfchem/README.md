# naqfc_verify_scripts
New command-line python scripts for regional NAQFC AIRNOW-CMAQ and global FV3-AERONET and FV3-Open-AQ (limited) verification using MONETv2.

Currently supported to work with MONETv2.1.4+.

Invoking each script with -help will show the list of arguments available.

AIRNOW chemical and meteorological verification is avaialable for the following species:

OZONE, PM2.5, PM10, CO, NO, NO2, SO2, NOX, NO2Y, TEMP, WS, WD, SRAD, BARPR, PRECIP, RHUM

OPEN-AQ chemical verification is avaialable for the following species:

pm25_ugm3, pm10_ugm3

AERONET chemical verification is available for the following species:

aod_550nm

MONETv2 is found at: https://github.com/noaa-oar-arl/MONET

See the naqfc_test.sh script for simple example of running some of the scripts (many options/arguments are available).
