module mo_write

  implicit none

  PRIVATE

  PUBLIC :: open_netcdf

CONTAINS

  subroutine open_netcdf(nc, var_data, var_time)

    use mo_kind, only: i4, sp, dp
    use mo_netcdf, only: NcDataset, NcDimension, NcVariable
    use mo_string_utils, only: num2str
    use mainVar, only: gridMeteo, yStart, mStart, dStart
    use NetCDFVar, only: fileOut, author_name, variable_name, variable_unit, variable_long_name, projection_name,invert_y

    implicit none

    type(NcDataset),  intent(out) :: nc
    type(NcVariable), intent(out) :: var_time, var_data
    
    type(NcDimension)     :: dim_x, dim_y, dim_time
    type(NcVariable)      :: var_east, var_north
    integer(i4)           :: i, f
    real(dp), allocatable :: dummy(:, :)

    ! 1.1 create a file
    nc = NcDataset(trim(fileOut), "w")

    ! create dimensions
    dim_x    = nc%setDimension("x", gridMeteo%ncols)
    dim_y    = nc%setDimension("y", gridMeteo%nrows)
    dim_time = nc%setDimension("time", -1)

    ! create variables
    var_time  = nc%setVariable('time', "i32", (/dim_time/))
    ! add some variable attributes
    call var_time%setAttribute("units", "days since " // trim(num2str(yStart, form='(I4)')) // "-"// &
        trim(num2str(mStart, form='(I0.2)')) // "-" // &
        trim(num2str(dStart, form='(I0.2)')) // "-" // "00:00:00")

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
    
    var_data = nc%setVariable(variable_name, "f64", (/dim_x, dim_y, dim_time/))
    ! add some more variable attributes
    call var_data%setAttribute("units",   trim(variable_unit))
    call var_data%setAttribute("long_name", trim(variable_long_name))
    call var_data%setAttribute("scaling", 0.1_dp)
    call var_data%setAttribute("missing_value", -9999._dp)

    ! add global attributes
    call nc%setAttribute("Author", trim(author_name))
    call nc%setAttribute("Projection", trim(projection_name))
    
  end subroutine open_netcdf

end module mo_write
