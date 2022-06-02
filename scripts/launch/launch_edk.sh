#!/bin/bash

set -e
# set -x
############################################################
##script for copying or creating folders and updating data##
############################################################
# created by Matthias Zink 2011
# adapted by Friedrich Boeing Apr. 2020

stationpath=${1}
datapath=${2}
start_date=${3}
end_date=${4}
var=${5}         # name of input variable
outpath=${6}


PWD=$(dirname $(readlink -f "${BASH_SOURCE[0]}"))
EDK_PROG_PATH="/data/duerre/GDM_new/local/Progs/edk_latlon/edk"
# decide to run the job on the frontend (interactive='yes') or on the grid engine ('no')
interactive='no'
# define estimated runtime and memeory requirements for the grid engine
rtime='10:00:00'
vmem='12G' 
mail='friedrich.boeing@ufz.de'
hmem='true'
OMP_NTHREAD='10'
qsub_manual='false'
# Define directories and presettings
start_year=$(echo ${start_date} | cut -d '-' -f 1)
start_month=$(echo ${start_date} | cut -d '-' -f 2)
start_day=$(echo ${start_date} | cut -d '-' -f 3)
end_year=$(echo ${end_date} | cut -d '-' -f 1)
end_month=$(echo ${end_date} | cut -d '-' -f 2)
end_day=$(echo ${end_date} | cut -d '-' -f 3)

# directory of input - DEM
dem="${datapath}/dem_hicam_1km_g.nc"
# path and filename of station look-up table
lut="${stationpath}/StationLuT_${var}.txt"
# cellsize = cellfactor * Gridcellsize_DEM (e.g. 100m)
cellfactor=1
# data converting factor --> value(double) = value(integer) * dataconv
dataconv=1.0e-1
# VARIOGRAM-estimate please see notes below
vario=".FALSE."
# variogarm type (spherical or exponential)
var_type=-9
# number of variogram parameters
var_param=3
# path and filename of first estimation if vario= ".false."
# or calculated variogram parameters    if vario= ".true."
vario_fname="${datapath}/var_param_de_${var}.txt"
# parameters for variogram estimation
# width [m] for binning the datapoints of the empirical variogram
binsize=1.0e3
# maximum search distance for data pairs
hmax=8.0e4
# create variables needed for processing
case ${var} in
    'pre')
    datatype="Precipitation"
    FlMethod=2   # see below
    no_data=-9   # nodata value
    buffer=1.0e5 # maximal distance of stations taken into consideration for interpolation [m]
    var_name='pre'
    var_unit='mm d-1'
    var_long_name='daily sum of precipitation'
    corrNeg='.True.'
    distZ='.True.'
    ;;
    'tavg') 
    datatype="average_Temp"
    FlMethod=2
    no_data=-999
    buffer=1.5e5
    var_name='tavg'
    var_unit='degC'
    var_long_name='mean daily temperature'
    corrNeg='.False.'
    distZ='.False.'
    ;;   
    'tmax') 
    datatype="max_Temp"
    FlMethod=2
    no_data=-999
    buffer=1.5e5
    var_name='tmax'
    var_unit='degC'
    var_long_name='maximum daily temperature'
    corrNeg='.False.'
    distZ='.False.'
    ;;   
    'tmin') 
    datatype="min_Temp"
    FlMethod=2
    no_data=-999
    buffer=1.5e5
    var_name='tmin'
    var_unit='degC'
    var_long_name='minimum daily temperature'
    corrNeg='.False.'
    distZ='.False.'
    ;;
    'rh') 
    datatype="Relative Humidity"
    FlVar=3
    FlMethod=2
    no_data=-9
    buffer=1.5e5
    ;;
    'windspeed') 
    datatype="horizontal windspeed"
    FlVar=4
    FlMethod=1
    no_data=-9
    buffer=1.5e5
    ;;
    *) 'ERROR: Please enter a valid var!'
esac

# *************************************************************
# ***NOTES***
# RUN FLAGS - 
# FlagVarType - FlVar    1 PRE
#                        2 TEM
#                        3 RU
#                        4 windspeed
#
# FalgMthType - FlMehtod 2 - EDK
#                        1 - OK                         
#                        0 - no interpolation
# Variogram Type - var_type (vType)
# var_type     1  : composed:   nugget + spherical   + sill
#             2  : composed:   nugget + exponential + sill
#
# VARIOGRAM-estimate: true:  calculation of variogram parameters
#                     false: read parameters from file
# 
# *************************************************************

outpath=${outpath}/

# check if the blocks do exist
if [ -d ${outpath} ] ; then 
    echo "Path: " ${outpath}" does already exists"
else
    mkdir -p ${outpath}
fi

# create main file in the output folder
mainfile=${outpath}/'edk.nml'
#mainfile='main.dat'
  # create mainfile containing filedir's etc.
[[ -f "${outpath}/edk.nml" ]] && rm "${outpath}/edk.nml"
cat > "${outpath}/edk.nml" << EOF
!EDK-Grid***Main configuration file for Precipitation for GER***
&mainVars ! namelist mainpaths
! -------------------------------------------------------------------
! ----------------- INPUT SPECIFICATIONS ----------------------------
!
! DataPathIn contains the path to the meteorological input files
! it can be a netcdf file, then a 3-dim field is expected and
! further coordinate and variable names are expected
DataPathIn               = "${stationpath}"
!
! value in the meteorological input files that is used as missing value
noDataValue              = ${no_data}
! Specifications of variable and coordinate names in case DataPathIn
! is a netcdf file
ncIn_variable_name       = "${var_name}"
ncIn_dem_variable_name   = "dem"
! name of coordinates, these must be in [m]
ncIn_yCoord_name         = "northing"
ncIn_xCoord_name         = "easting"
!
! Correct interpolated negative values to zero 
! (set for precipitation to true !)
correctNeg               = ${corrNeg}
!
! Set values further away than the distance threshold to zero
distZero                 = ${distZ}
!
cellFactor               = ${cellfactor}

fNameDEM                 = "${dem}"
! name of coordinates (northing and easting) in [m]
ncOut_dem_yCoord_name    = "northing"
ncOut_dem_xCoord_name	 = "easting"
! name of coordinates (latitude and longitude) in deg
ncOut_dem_Latitude	 = "latitude"
ncOut_dem_Longitude	 = "longitude"
! name of the DEM variable in the netcdf file
ncOut_dem_variable_name  = "dem"
!
DataPathOut              = "${outpath}"
FileOut                  = "${var_name}.nc"
!
! Name of Look Up Table of Station data
fNameSTA                 = "${lut}"
DataConvertFactor        = ${dataconv}
! The value that should be added to the data in the netcdf file 
OffSet			 = 0
![FINAL VALUE = NETCDF VALUE * DataConvertFactor + OffSet]


!
! -------------------------------------------------------------------
! ------------ PROCESSING PERIOD ------------------------------------ 
yStart                   = ${start_year}
mStart                   = ${start_month}
dStart                   = ${start_day}
yEnd                     = ${end_year}
mEnd                     = ${end_month}
dEnd                     = ${end_day}
! Number of Time Buffering Days (Divides the EDK processing into chunks accross the time dimension) 
tBuffer                  = 30

!
! -------------------------------------------------------------------
! ------------ INTERPOLATION METHOD ---------------------------------
!
! InterMth = 2 -> EDK
! InterMth = 1 -> OK
! InterMth = 0 -> No interpolation
InterMth                 = ${FlMethod}
! maximum search distance for interpolation [m]
maxDist                  = ${buffer}
!
! -------------------------------------------------------------------
! ----------------- VARIOGRAM ESTIMATION ----------------------------
flagVario                = ${vario}
!
! number of variogram parameters
nParam                   = ${var_param}
! type of theoretical variogramm (1=spherical, 2 = exponential)
vType                    = ${var_type}
!
! file name where to store the variogram parameters 
! (if flagVario is false, variogram parameters )
! (for interpolation are read in here:         )
fNameVario               = "${vario_fname}" 
dh                       = ${binsize} ! binsize for variogram [m]
hMax                     = ${hmax} ! max distance h for variogram estimation [m]
! -------------------------------------------------------------------
! --------------- NC OUTPUT SPECIFICATION ----------------------------
author_name              = 'Friedrich Boeing'
projection_name          = 'EPSG: 4326'
invert_y                 = .True. ! (set True if working with mHM input data!)

variable_name            = "${var_name}"
variable_unit            = "${var_unit}"
variable_long_name       = "${var_long_name}"
variable_standard_name	 = 'precipitation_flux'
variable_calendar_type	 = 'proleptic_gregorian'
! -------------------------------------------------------------------
! -------------------------------------------------------------------
/ !END*******Main Config***********
EOF





joblabel="edk_"${var}_${start_date}_${end_date}
outfile=${var}_${start_date}_${end_date}'.output'

# create Readme file in output directory
printf "Interpolation Started "             > ${outpath}/README
date "+on %m/%d/%Y at %H:%M:%S %n"         >> ${outpath}/README
cat ${mainfile}                            >> ${outpath}/README
printf "\n \n>****Variogram Parameters \n" >> ${outpath}/README
cat ${vario_fname}                          >> ${outpath}/README

if [ ${interactive} == "no" ] ; then
    # create script for submission to the cluster
    fSubmit="${outpath}/qsub_tmp.sh"
    if [ -f ${fSubmit} ] ; then rm ${fSubmit} ; fi  
    cat > ${fSubmit} << EOF
#!/bin/bash
#-------------------------------
#\$ -S /bin/bash
#\$ -N ${joblabel}
#\$ -o ${outpath}/${outfile}.\$JOB_ID
#\$ -j y
#\$ -l h_rt=${rtime}
#\$ -l h_vmem=${vmem}
#\$ -l highmem=${hmem}
#\$ -cwd
#\$ -pe smp ${OMP_NTHREAD}
#\$ -m ea      # mail notification in case of job ended abortion 
#\$ -M ${mail} # mail address
#-------------------------------
ulimit -s unlimited

export OMP_NUM_THREADS=${OMP_NTHREAD}
cd ${outpath} || exit


ml foss/2018b
ml netCDF-Fortran

time  ./edk

# rm edk
# rm edk.nml
EOF
    ln -fs ${EDK_PROG_PATH} ${outpath}
    cd ${outpath}
    if [[ $qsub_manual == "true" ]]; then
        echo "..!manual job submission to cluster!..."
    else
        qsub qsub_tmp.sh
        echo "Finished submission"
    fi
else
    # run program on frontend, lauch it from the output directory
    ln -fs ${EDK_PROG_PATH} ${outpath}
    cd ${outpath}

    export OMP_NUM_THREADS=${OMP_NTHREAD}
    ml foss/2018b
    ml netCDF-Fortran
    time ./edk
    # rm edk
    # rm edk.nml
    echo "Finished submission"
fi
  # change back to root directory before next job starts/is submitted
cd ${PWD}


exit 0
