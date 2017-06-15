#!/bin/bash

set -e

  ############################################################
  ##script for copying or creating folders and updating data##
  #######################mz2011-02-01#########################
  # decide to run the job on the frontend (interactive='yes') or on the grid engine ('no')
interactive='no'
# define estimated runtime and memeory requirements for the grid engine
rtime='32:00:00' # 540 stations (tavg) ca 7h per year
vmem='4G' #temp 2G pre 4G
# Define directories and presettings
# prefix describing data, in- and output directories
prefix='pre'
# directory of the edk program - executeable
rundir=$(pwd)'/'
# directory of input - data files
inpath="/data/stohyd/data/processed/Germany/DWD/${prefix}/1900-2014_02_23/"
# directory of input - DEM
DEMfile="/data/stohyd/data/processed/Germany/DEM/filled_final1.asc"
#DEMfile="/home/zink/data/edk/dem.test"
# path and filename of station look-up table
LuT="/data/stohyd/data/processed/Germany/DWD/${prefix}/1900-2014_02_23/StationLuT_${prefix}.txt"
# target directory
outpath="/work/zink/edk_crossloc/vario/" #${i}/"
# cellsize = cellfactor * Gridcellsize_DEM (e.g. 100m)
cellfactor=40
# data converting factor --> value(double) = value(integer) * dataconv
dataconv=1.0e-1
# datat format of the output either 'bin' or 'nc'
dataformat="bin"
# flag for interpolation - ".true."=interpolate or ".false."=skip interpolation
edkesti=".FALSE."
# starting year of interpolation
edkyStart=1950
edkmStart=01
edkdStart=01
# ending date of interpolation
edkyEnd=2010
edkmEnd=12
edkdEnd=31
# VARIOGRAM-estimate please see notes below
vario=".TRUE."
# variogarm type (spherical or exponential)
varType=2
# number of variogram parameters
varParam=3
# path and filename of first estimation if vario= ".false."
# or calculated variogram parameters    if vario= ".true."
#varioFname="/home/zink/source_code/edk/variogram_param_apriori.txt"
varioFname="/data/stohyd/data/processed/Germany/DWD/${prefix}/var_param_de_${prefix}.txt"
# parameters for variogram estimation
# width [m] for binning the datapoints of the empirical variogram
binsize=1.0e3
# maximum search distance for data pairs
hMAx=8.0e4
# name of the catchment
catchment="GER"
# create variables needed for processing
case ${prefix} in
  'pre')
    datatype="Precipitation"
    FlVar=1      # see below
    FlMethod=1   # see below
    no_data=-9   # nodata value
    buffer=1.0e5 # maximal distance of stations taken into consideration for interpolation [m]
    var_name='pre'
    var_unit='mm d-1'
    var_long_name='daily sum of precipitation'
    ;;
  'tavg') 
    datatype="average_Temp"
    FlVar=2
    FlMethod=1
    no_data=-999
    buffer=1.5e5
    var_name='tavg'
    var_unit='degC'
    var_long_name='mean daily temperature'
    ;;   
  'tmax') 
    datatype="max_Temp"
    FlVar=2
    FlMethod=1
    no_data=-999
    buffer=1.5e5
    var_name='tmax'
    var_unit='degC'
    var_long_name='maximum daily temperature'
    ;;   
  'tmin') 
    datatype="min_Temp"
    FlVar=2
    FlMethod=1
    no_data=-999
    buffer=1.5e5
    var_name='tmin'
    var_unit='degC'
    var_long_name='minimum daily temperature'
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
  *) 'ERROR: Please enter a valid prefix!'
esac

# *************************************************************
# ***NOTES***
# RUN FLAGS - 
# FlagVarType - FlVar    1 PRE
#                        2 TEM
#                        3 RU
#                        4 windspeed
#
# FalgMthType - FlMehtod 1 - EDK
#                        2 - OK                         
#
# Variogram Type - varType (vType)
# varType     1  : composed:   nugget + spherical   + sill
#             2  : composed:   nugget + exponential + sill
#
# VARIOGRAM-estimate: true:  calculation of variogram parameters
#                     false: read parameters from file
# 
# *************************************************************
cd ${rundir}

joblabel="edk_"${prefix}_${edkyStart}_${edkyEnd}
outfile=${prefix}_${edkyStart}_${edkyEnd}'.output'

  # check if the blocks do exist
echo "Processing year : ${edkyStart}"
if [ -f "${outpath}${edkyStart}.bin" ] ; then 
  echo "***CYCLING: year.bin already existing"
  continue
      # exit 2
fi

outpath=${outpath}/

  # check if the blocks do exist
if [ -d ${outpath} ] ; then 
  echo "Path: " ${outpath}" does alreaddy exists"
  continue # exit 2
else
  mkdir -p ${outpath}
fi

  # create main file in the output folder
mainfile=${outpath}/'main.dat'
  #mainfile='main.dat'

  # create mainfile containing filedir's etc.
if [ -f ${mainfile} ] ; then rm ${mainfile} ; fi
cat > ${mainfile} << EOF
!EDK-Grid***Main configuration file for ${datatype} for ${catchment}***
&mainVars ! namelist mainpaths
flagMthTyp               = ${FlMethod}
flagVarTyp               = ${FlVar}
noDataValue              = ${no_data}
DataPathIn               = "${inpath}"
fNameDEM                 = "${DEMfile}"
DataPathOut              = "${outpath}"
fNameSTA                 = "${LuT}"
cellFactor               = ${cellfactor}
DataConvertFactor        = ${dataconv}
outputformat             = "${dataformat}"
variable_name            = "${var_name}"
variable_unit            = "${var_unit}"
variable_long_name       = "${var_long_name}"
flagEDK                  = ${edkesti}
yStart                   = ${edkyStart}
mStart                   = ${edkmStart}
dStart                   = ${edkdStart}
yEnd                     = ${edkyEnd}
mEnd                     = ${edkmEnd}
dEnd                     = ${edkdEnd}
maxDist                  = ${buffer}
flagVario                = ${vario}
vType                    = ${varType}
nParam                   = ${varParam}
fNameVario               = "${varioFname}"                 ! path+filename 
dh                       = ${binsize}                      ! binsize for variogram
hMax                     = ${hMAx}                         ! max distance h
/ !END*******Main Config***********
EOF

# create Readme file in output directory
printf "Interpolation Started "             > ${outpath}/README
date "+on %m/%d/%Y at %H:%M:%S %n"         >> ${outpath}/README
cat ${mainfile}                            >> ${outpath}/README
printf "\n \n>****Variogram Parameters \n" >> ${outpath}/README
cat ${varioFname}                          >> ${outpath}/README

if [ ${interactive} == "no" ] ; then
      # create script for submission to the cluster
  fSubmit="${outpath}/qsub_tmp.sh"
  if [ -f ${fSubmit} ] ; then rm ${fSubmit} ; fi  
  cat > ${fSubmit} << EOF
#!bin/bash
#-------------------------------
#\$ -S /bin/bash
#\$ -N ${joblabel}
#\$ -o ${outpath}/${outfile}.\$JOB_ID
#\$ -j y
#\$ -l h_rt=${rtime}
#\$ -l h_vmem=${vmem}
#\$ -cwd
#\$ -binding linear:1
#-------------------------------
ulimit -s unlimited
export OMP_NUM_THREADS=1

cd ${outpath}
time  ./edk

rm edk
rm main.dat
EOF
      #rm ./edk
      # copy executeable to the output directory to run it from there (uniqe executeable)
  cp edk ${outpath}
  cd ${outpath}      
  qsub qsub_tmp.sh
else
      # run program on frontend, lauch it from the output directory
  cp edk ${outpath}
  cd ${outpath}
  time ./edk
  rm edk
  rm main.dat
fi
  # change back to root directory before next job starts/is submitted
cd ${rundir}

echo "Finished submission"

exit 0
