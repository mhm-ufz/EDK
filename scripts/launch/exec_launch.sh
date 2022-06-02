#!/bin/bash
set -e    # will stop script when error occurs / return_value != 0
set -x

###############################################
### Script to execute the edk launch script ###
###############################################

# set main paths and variables here
# further settings for job submission etc. in the launch_edk.sh

PWD=$(dirname $(readlink -f "${BASH_SOURCE[0]}"))
DWDPATH="/data/klimabuero/GDM/tmp/dwd/" # path containing the station data
EDKSTATIC="/data/klimabuero/GDM/data/EDK/" # path containing the variogram parameters and dem
OUTPATH_MAIN="${PWD}/2020/" # outpath, can be decided!
start_date="2020-01-01"
end_date="2020-05-02"

METEOVARs=("pre" "tavg" "tmin" "tmax")
METEOVARs=("pre")
for var in "${METEOVARs[@]}"; do  
    bash ${PWD}/launch_edk.sh "${DWDPATH}/${var}/" "${EDKSTATIC}"  "${start_date}" "${end_date}" "${var}" "${OUTPATH_MAIN}/${var}"
done
