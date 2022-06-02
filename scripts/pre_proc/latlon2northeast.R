#Script to convert lat lon coordinates and add extra "northing" and "easting" variables to the EOBS netcdf files
# Author - Akash Koppa 
# Date - 2/10/2020

# clear workspace
rm(list=ls())

# read in arguments from the command line
args      = commandArgs(trailingOnly=TRUE)
file_reqd  = args[1]
varid_reqd = args[2]
crs_reqd  = args[3] 

# load required libraries
library(ncdf4)
library(raster)

# read in the lat and lon coordinates
id = nc_open(file_reqd)
data = ncvar_get(id, varid = varid_reqd)
northing = ncvar_get(id, varid = "latitude")
easting  = ncvar_get(id, varid = "longitude")
nc_close(id)

# create a two dimensional matrix of latitudes and longitudes
north2d = NULL
for(i in 1:length(easting)){
  north2d = rbind(north2d, northing)
}

east2d = NULL
for (i in 1:length(northing)){
  east2d = cbind(east2d, easting)
}

# convert latitude and longitude  into northing and easting 
coord_latlon = data.frame(Longitude = c(east2d),
                          Latitude  = c(north2d))
coord_latlon = coordinates(coord_latlon)
coord_latlon = SpatialPoints(coord_latlon, CRS("+init=epsg:4326"))
coord_proj   = spTransform(coord_latlon, CRS(crs_reqd))
coord_reqd   = coordinates(coord_proj)
north2d_reqd = matrix(data = coord_reqd[,2],nrow=nrow(data))
east2d_reqd  = matrix(data = coord_reqd[,1],nrow=nrow(data))

# write these variables into the netcdf files
id = nc_open(file_reqd,write=TRUE)
xdim = id$dim[['longitude']]
ydim = id$dim[['latitude']]

var_north2d = ncvar_def(name = "northing",units = "meter",dim = list(xdim, ydim), missval = -9999,longname = paste0("Northing in ", crs_reqd))
var_east2d = ncvar_def(name = "easting",units = "meter",dim = list(xdim, ydim), missval = -9999,longname = paste0("Easting in ", crs_reqd))

# put variables
id = ncvar_add(id, var_north2d)
id = ncvar_add(id, var_east2d)
ncvar_put(id, var_north2d, north2d_reqd)
ncvar_put(id, var_east2d, east2d_reqd)
nc_close(id)
