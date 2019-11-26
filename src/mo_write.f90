module mo_write

  implicit none

  PRIVATE

  PUBLIC :: open_netcdf

CONTAINS

  subroutine open_netcdf(nc, var_data, var_time)

    use mo_kind, only: i4, sp, dp
    use mo_netcdf, only: NcDataset, NcDimension, NcVariable
    use mo_string_utils, only: num2str
    use mainVar, only: gridMeteo, yStart
    use NetCDFVar, only: fileOut, author_name, variable_name

    implicit none

    type(NcDataset),  intent(out) :: nc
    type(NcVariable), intent(out) :: var_time, var_data
    
    type(NcDimension) :: dim_x, dim_y, dim_time
    type(NcVariable)  :: var_lon, var_lat

    ! 1.1 create a file
    nc = NcDataset(trim(fileOut), "w")

    ! create dimensions
    dim_x    = nc%setDimension("x", gridMeteo%ncols)
    dim_y    = nc%setDimension("y", gridMeteo%nrows)
    dim_time = nc%setDimension("time", -1)

    ! create variables
    var_time = nc%setVariable('time', "i32", (/dim_time/))
    ! var_lat  = nc%setVariable(vname_lat,  "f32", (/dim_x, dim_y/))
    ! var_lon  = nc%setVariable(vname_lon , "f32", (/dim_x, dim_y/))
    var_data = nc%setVariable(variable_name, "f64", (/dim_y, dim_x, dim_time/))

    ! add some variable attributes
    call var_time%setAttribute("units", "days since " // trim(num2str(yStart - 1, form='(I4)')) // "-12-31 12:00:00")

    ! ! write data of static variables
    ! call var_lat%setData(wlat)
    ! call var_lon%setData(wlon)
    
    ! add some more variable attributes
    call var_data%setAttribute("units",   "mm/d")
    call var_data%setAttribute("scaling", 0.1_dp)
    call var_data%setAttribute("missing_value", -9999._dp)

    ! add global attributes
    call nc%setAttribute("Author", trim(author_name))

  end subroutine open_netcdf

end module mo_write
