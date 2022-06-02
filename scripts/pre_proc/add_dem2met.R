# Script to add DEM data to the EOBS netcdf file
# Author - Akash Koppa
# Date - 03/04/2020

# clear workspace 
rm(list=ls())

# read command line arguments
args = commandArgs(trailingOnly=TRUE)
met_file = args[1]
dem_file = args[2]

# load required libraries
library(ncdf4)
library(raster)

# read in the 10km dem for Europe (E-OBS grid)
id = nc_open(dem_file)
dem = ncvar_get(id, varid="dem")
nc_close(id)

# write the DEM variable into the netcdf files 
id = nc_open(met_file,write=TRUE)
xdim = id$dim[['longitude']]
ydim = id$dim[['latitude']]

# create the dimension information for the netcdf files
var_dem = ncvar_def(name = "dem",units="meter",dim=list(xdim, ydim),missval=-9999,longname="Elevation")

# put variables
id = ncvar_add(id, var_dem)
ncvar_put(id, var_dem, dem)
nc_close(id)

