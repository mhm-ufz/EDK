#!/bin/bash
set -e
# qsub -t n-m launch_edk_crossloc_de.sh
#-------------------------------
#$ -S /bin/bash
#$ -N edk_cross
#$ -o $JOB_NAME-$TASK_ID.out
#$ -j y
#$ -l h_rt=00:45:00
#$ -l h_vmem=1G
#$ -cwd
#$ -binding linear:1
#-------------------------------
ulimit -s unlimited
export OMP_NUM_THREADS=1


############################################################
##script for copying or creating folders and updating data##
#######################mz2011-02-01#########################
# decide to run the job on the frontend (interactive='yes') or on the grid engine ('no')
# define estimated runtime and memeory requirements for the grid engine
# Define directories and presettings
# prefix describing data, in- and output directories
prefix='pre'
# directory of the edk program - executeable
rundir=$(pwd)'/'
# directory of input - data files
inpath="/data/stohyd/data/processed/Germany/DWD/${prefix}/${prefix}_1900-2015/"
# path and filename of station look-up table
LuT="/data/stohyd/data/processed/Germany/DWD/${prefix}/pre_1900-2015/StationLuT_${prefix}.txt" # 5394 stations
# target directory
pathout="/work/zink/edk_crossloc/${prefix}/"
# cellsize = cellfactor * Gridcellsize_DEM (e.g. 100m)
cellfactor=1
# data converting factor --> value(double) = value(integer) * dataconv
dataconv=1.0e-1
# datat format of the output either 'bin' or 'nc'
dataformat="nc"
# flag for interpolation - ".true."=interpolate or ".false."=skip interpolation
edkesti=".TRUE."
# starting year of interpolation
edkyStart=1950
edkmStart=01
edkdStart=01
# ending date of interpolation
edkyEnd=2010
edkmEnd=12
edkdEnd=31
# VARIOGRAM-estimate please see notes below
vario=".FALSE."
# variogarm type (spherical or exponential)
varType=2
# number of variogram parameters
varParam=3
# path and filename of first estimation if vario= ".false."
# or calculated variogram parameters    if vario= ".true."
#varioFname="/home/zink/source_code/fortran/edk_nc/variogram_param_apriori.txt"
varioFname="/data/stohyd/data/processed/Germany/DWD/${prefix}/var_param_de_${prefix}.txt"
#varioFname="/data/stohyd/data/processed/Germany/DWD/${prefix}/var_param_de_${prefix}.txt"
# parameters for variogram estimation
# width [m] for binning the datapoints of the empirical variogram
binsize=1.0e3
# maximum search distance for data pairs
hMAx=8.0e4
# name of the catchment
catchment="test"
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

no_stat=$(( $(grep 'Number_of_Stations:' ${LuT} | cut -d ':' -f 2)))

stationlist=($(cat ${LuT} | sed '1,2d' | tr -s ' ' | cut -d ' ' -f 2))

#for ((SGE_TASK_ID=1 ; SGE_TASK_ID<=2 ; SGE_TASK_ID++ )) ; do


istation=${stationlist[$((${SGE_TASK_ID}-1))]}
echo 'processing station:' ${istation}

cd ${rundir}
  
outpath="${pathout}exclude_${istation}/"  
# check if the blocks do exist
if [ -d ${outpath} ] ; then 
  echo "Path: " ${outpath}" does alreaddy exists"
  exit 2
else
  mkdir -p ${outpath}
fi
echo ${outpath}

# shorten station look up table by one station - leave one out
LuT_new="${outpath}Station_LuT.txt"
no_stat_new=$(( $(grep 'Number_of_Stations:' ${LuT} | cut -d ':' -f 2)-1 ))
#echo "new ${no_stat_new}"
echo "Number_of_Stations: ${no_stat_new}" > ${LuT_new}
sed "/^ \+${istation}/d" ${LuT}  >> ${LuT_new}
stat_details=$(sed -n "/^ \+${istation}/p" ${LuT})
stat_east=$(echo ${stat_details} | cut -d ' ' -f 2)
stat_nort=$(echo ${stat_details} | cut -d ' ' -f 3)
stat_heig=$(echo ${stat_details} | cut -d ' ' -f 4)
sed -i'' '2d' ${LuT_new}

# create dem file in which only station height is active cell
DEMfile="${outpath}station_height_${istation}.asc"
cat > ${DEMfile} << EOF  
ncols                  1
nrows                  1
xllcorner   ${stat_east}
yllcorner   ${stat_nort}
cellsize               1
NODATA_value       -9999
${stat_heig}
EOF
  
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
fNameSTA                 = "${LuT_new}"
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
EOF

# create Readme file in output directory
printf "Interpolation Started "             > ${outpath}/README
date "+on %m/%d/%Y at %H:%M:%S %n"         >> ${outpath}/README
cat ${mainfile}                            >> ${outpath}/README
printf "\n \n>****Variogram Parameters \n" >> ${outpath}/README
cat ${varioFname}                          >> ${outpath}/README

cp ../edk ${outpath}

cd ${outpath}
time  ./edk

# load module system
source /etc/profile.d/000-modules.sh
# load python environment
module load /global/apps/chs-virtualenv/chspython/2.7.6

export PYTHONPATH=/home/zink/lib/PYTHON_chs_lib/:$PYTHONPATH

python /home/zink/source_code/fortran/edk_nc/cross_loc_test/interpolation_error_crosscheck_single_stations.py \
       -d ${inpath} -l ${istation} -n ${pathout}

# change back to root directory before next job starts/is submitted
cd ${rundir}

#done

echo "Finished submission"

exit 0
