module mo_write

  implicit none

  PRIVATE

  PUBLIC :: open_netcdf

CONTAINS

  subroutine open_netcdf(nc, var_data, var_time)

    use mo_kind, only: i4, sp, dp
    use mo_netcdf, only: NcDataset, NcDimension, NcVariable
    use mo_string_utils, only: num2str
    use mainVar, only: gridMeteo, yStart, mStart, dStart, DEMNcFlag, DataConvertFactor
    use NetCDFVar, only: fileOut, author_name, variable_name, variable_unit, variable_long_name, &
                         projection_name,invert_y, variable_standard_name, variable_calendar_type, &
                         ncOut_dem_Latitude, ncOut_dem_Longitude

    implicit none

    type(NcDataset),  intent(out) :: nc
    type(NcVariable), intent(out) :: var_time, var_data
    
    type(NcDimension)     :: dim_x, dim_y, dim_time
    type(NcVariable)      :: var_east, var_north, var_lat, var_lon
    integer(i4)           :: i, f
    real(dp), allocatable :: dummy(:, :)
    real(dp), allocatable :: dummy_lat(:)

    ! 1.1 create a file
    nc = NcDataset(trim(fileOut), "w")
   
    if (DEMNcFlag == 1) then
      dim_y    = nc%setDimension(ncOut_dem_Latitude, gridMeteo%ncols)
      dim_x    = nc%setDimension(ncOut_dem_Longitude, gridMeteo%nrows)
      dim_time = nc%setDimension("time", -1)
    else 
    ! create dimensions
      dim_x    = nc%setDimension("x", gridMeteo%ncols)
      dim_y    = nc%setDimension("y", gridMeteo%nrows)
      dim_time = nc%setDimension("time", -1)
    end if

    ! create variables
    var_time  = nc%setVariable('time', "i32", (/dim_time/))
    ! add some variable attributes
    call var_time%setAttribute("units", "days since " // trim(num2str(yStart, form='(I4)')) // "-"// &
        trim(num2str(mStart, form='(I0.2)')) // "-" // &
        trim(num2str(dStart, form='(I0.2)')) // " " // "00:00:00")
    call var_time%setAttribute("calendar", variable_calendar_type)

  if (DEMNcFlag == 1) then

    if (invert_y) then

       allocate(dummy(gridMeteo%nrows, gridMeteo%ncols))
       f = 1
       do i = 1, gridMeteo%ncols
         dummy(:,gridMeteo%ncols - i + 1) = gridMeteo%northing(:,f) 
         f = f + 1 
       end do
       var_north = nc%setVariable('northing', 'f32', (/dim_x, dim_y/))
       call var_north%setAttribute("standard_name", "northing")
       call var_north%setAttribute("units", "m")
       call var_north%setData(dummy)

    else
       var_north = nc%setVariable('northing', 'f32', (/dim_x, dim_y/))
       call var_north%setData(gridMeteo%northing)
       call var_north%setAttribute("standard_name", "northing")
       call var_north%setAttribute("units", "m")
    end if

    var_east = nc%setVariable('easting', 'f32', (/dim_x, dim_y/))
    call var_east%setData(gridMeteo%easting)
    call var_east%setAttribute("standard_name", "easting")
    call var_east%setAttribute("units", "m")
     
    if (invert_y) then
      allocate(dummy_lat(gridMeteo%ncols)) 
      f = 1
      do i = 1, gridMeteo%ncols
         dummy_lat(gridMeteo%ncols - i + 1) = gridMeteo%latitude(f)
         f = f + 1
      end do
      var_lat = nc%setVariable(ncOut_dem_Latitude, 'f32', (/dim_y/))
      call var_lat%setAttribute("standard_name", "latitude")
      call var_lat%setAttribute("units", "degrees_north")
      call var_lat%setData(dummy_lat)
    else
      var_lat = nc%setVariable(ncOut_dem_Latitude, 'f32', (/dim_y/))
      call var_lat%setAttribute("standard_name", "latitude")
      call var_lat%setAttribute("units", "degrees_north")
      call var_lat%setData(gridMeteo%latitude)
    end if
    var_lon = nc%setVariable(ncOut_dem_Longitude, 'f32', (/dim_x/))
    call var_lon%setAttribute("standard_name", "longitude")
    call var_lon%setAttribute("units", "degrees_east")
    call var_lon%setData(gridMeteo%longitude)

  else
    allocate(dummy(gridMeteo%ncols, gridMeteo%nrows))
    if (invert_y) then
       f = 1
       do i = gridMeteo%nrows, 1, -1
          dummy(:, i) = gridMeteo%yllcorner + (real(f,dp) - 0.5_dp) * gridMeteo%cellsize
          f = f + 1
       end do
    else
       do i = 1, gridMeteo%nrows
          dummy(:, i) = gridMeteo%yllcorner + (real(i,dp) - 0.5_dp) * gridMeteo%cellsize
       end do
    end if
    var_north = nc%setVariable('northing',  "f32", (/dim_x, dim_y/))
    call var_north%setAttribute("standard_name", "northing")
    call var_north%setAttribute("units", "m")
    call var_north%setData(dummy)

    do i = 1, gridMeteo%ncols
      dummy(i, :) = gridMeteo%xllcorner + (real(i,dp) - 0.5_dp) * gridMeteo%cellsize
    end do
    var_east  = nc%setVariable('easting', "f32", (/dim_x, dim_y/))
    call var_east%setAttribute("standard_name", "easting")
    call var_east%setAttribute("units", "m")
    call var_east%setData(dummy)
    deallocate(dummy)
    if( allocated(dummy_lat) ) deallocate(dummy_lat)

  end if
  !var_data = nc%setVariable(variable_name, "f32", (/dim_x, dim_y, dim_time/))
  var_data = nc%setVariable(variable_name, "f32",  (/dim_x, dim_y, dim_time/))
  ! add some more variable attributes
  call var_data%setAttribute("units",   trim(variable_unit))
  call var_data%setAttribute("long_name", trim(variable_long_name))
  call var_data%setAttribute("standard_name",trim(variable_standard_name))
  call var_data%setAttribute("scaling", 1.0_dp)
  call var_data%setAttribute("missing_value", -9999._dp)
    

  ! add global attributes
  call nc%setAttribute("Author", trim(author_name))
  call nc%setAttribute("Projection", trim(projection_name))
  
  end subroutine open_netcdf

end module mo_write
