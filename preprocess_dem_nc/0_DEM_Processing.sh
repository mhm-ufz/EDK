#!/bin/bash
# Bash script to preprocess DEM for EDK 
# Author: Akash Koppa
# Description: a single bash script to prepare the target DEM and the meteorological forcings for EDK in netcdf format
# Workflow: 1) Target DEM:  Upscale source DEM to target resolution ---> subset required domain ---> change variable names ---> assign missing value ---> add  easting and northing ---> Ready for EDK
#           2) Meteorological Forcing: ---> subset to target domain ---> add easting and northing ---> add DEM ---> Read for EDK


#------------------------USER DEFINED SETTINGS-------------------------#
## input and output paths
# absolute path to target DEM 
source_dem="/data/hicam/data/processed/de_hicam/morph/dem.asc"
# absolute path to meteorological files (currently supports mHM inputs)
path_rr="/data/hicam/EOBS/raw/latlon/rr_ens_mean_0.1deg_reg_v20.0e.nc"  # precipitation
path_tg="/data/hicam/EOBS/raw/latlon/tg_ens_mean_0.1deg_reg_v20.0e.nc"  # average temperature
path_tn="/data/hicam/EOBS/raw/latlon/tn_ens_mean_0.1deg_reg_v20.0e.nc"  # minimum temperature
path_tx="/data/hicam/EOBS/raw/latlon/tx_ens_mean_0.1deg_reg_v20.0e.nc"  # maximum temperature 
# output directory
output_dir="/work/koppa/hicam_edk/test"

# domain specifications
xmin="4.0" 
xmax="19.5"
ymin="45.0"
ymax="55.5"  # domain extent, format xmin xmax, ymin ymax 
target_res="0.015625"   # target resolution
met_res="0.10"
crs_reqd="+init=epsg:3035" # provide the epsg of the m-m coordinate system

# meteorological variable details (name of the variable in the netcdf file)
var_rr="rr"
var_tg="tg"
var_tn="tn"
var_tx="tx"



# --------------------PROCESSING STEPS---------------------------#
# load required modules (currently using  eve modules).
# Note: Change appropriately 
# load netcdf and R modules
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 NCO/4.7.8 CDO/1.9.5 GDAL/2.2.3-Python-3.6.6 ncview/2.1.7 R/3.6.0

## process target DEM
# convert the ascii file to .tif 
gdal_translate ${source_dem} "${output_dir}/dem_target.tif"

# upscale the DEM to the required resolution. Also convert to netcdf
gdalwarp -tr ${target_res} ${target_res} -r "average" -of netCDF "${output_dir}/dem_target.tif" "${output_dir}/dem_target.nc"

# use ncks to cut domain 
ncks -d lon,${xmin},${xmax} -d lat,${ymin},${ymax} -O "${output_dir}/dem_target.nc" "${output_dir}/dem_target.nc" 

# change the names of primary variable, lat and lon
ncrename -v Band1,dem -d lat,latitude -d lon,longitude -O "${output_dir}/dem_target.nc" "${output_dir}/dem_target.nc"
ncrename -v lat,latitude -O "${output_dir}/dem_target.nc" "${output_dir}/dem_target.nc" 
ncrename -v lon,longitude -O "${output_dir}/dem_target.nc" "${output_dir}/dem_target.nc" 

# set missing value to -9999
cdo setmissval,-9999 "${output_dir}/dem_target.nc" "${output_dir}/dem_target.nc"

# add northing and easting to the target variable using the R script
Rscript latlon2northeast.R "${output_dir}/dem_target.nc" "dem" ${crs_reqd}

## process input meteorological forcing
# subset the meteorological forcing
ncks -d longitude,${xmin},${xmax} -d latitude,${ymin},${ymax} "${path_rr}" "${output_dir}/rr_temp.nc" 
ncks -d longitude,${xmin},${xmax} -d latitude,${ymin},${ymax} "${path_tg}" "${output_dir}/tg_temp.nc"
ncks -d longitude,${xmin},${xmax} -d latitude,${ymin},${ymax} "${path_tn}" "${output_dir}/tn_temp.nc"
ncks -d longitude,${xmin},${xmax} -d latitude,${ymin},${ymax} "${path_tx}" "${output_dir}/tx_temp.nc"

# set missing value to -9999
cdo setmissval,-9999 "${output_dir}/rr_temp.nc" "${output_dir}/rr_final.nc"
cdo setmissval,-9999 "${output_dir}/tg_temp.nc" "${output_dir}/tg_final.nc"
cdo setmissval,-9999 "${output_dir}/tn_temp.nc" "${output_dir}/tn_final.nc"
cdo setmissval,-9999 "${output_dir}/tx_temp.nc" "${output_dir}/tx_final.nc"

# remove all temporary files
rm "${output_dir}/rr_temp.nc" "${output_dir}/tg_temp.nc" "${output_dir}/tn_temp.nc" "${output_dir}/tx_temp.nc"

# add northing and easting information
Rscript latlon2northeast.R "${output_dir}/rr_final.nc" ${var_rr} ${crs_reqd}
Rscript latlon2northeast.R "${output_dir}/tg_final.nc" ${var_tg} ${crs_reqd}
Rscript latlon2northeast.R "${output_dir}/tn_final.nc" ${var_tn} ${crs_reqd}
Rscript latlon2northeast.R "${output_dir}/tx_final.nc" ${var_tx} ${crs_reqd}

# remap dem to meteorological forcing
gdalwarp -tr ${met_res} ${met_res} -r "average" -of netCDF "${output_dir}/dem_target.tif" "${output_dir}/dem_met.nc"
ncrename -v Band1,dem "${output_dir}/dem_met.nc" "${output_dir}/dem_temp.nc"
cdo setmissval,-9999 "${output_dir}/dem_temp.nc" "${output_dir}/dem_met.nc"
rm "${output_dir}/dem_temp.nc"

# add dem to the meteorological forcing
Rscript add_dem2met.R "${output_dir}/rr_final.nc" "${output_dir}/dem_met.nc"
Rscript add_dem2met.R "${output_dir}/tg_final.nc" "${output_dir}/dem_met.nc"
Rscript add_dem2met.R "${output_dir}/tn_final.nc" "${output_dir}/dem_met.nc"
Rscript add_dem2met.R "${output_dir}/tx_final.nc" "${output_dir}/dem_met.nc"









  




