module mo_write

  implicit none

  PRIVATE

  PUBLIC :: open_netcdf

CONTAINS

  subroutine open_netcdf(fname, ncols, nrows, nc_id, var_id, var_time)

    use mo_kind, only: i4, sp, dp
    use mo_netcdf, only: NcDataset, NcDimension, NcVariable

    implicit none

    character(256),   intent(in)  :: fname
    integer(i4),      intent(in)  :: ncols
    integer(i4),      intent(in)  :: nrows
    type(NcDataset),  intent(out) :: nc
    type(NcVariable), intent(out) :: var_time, var_data
    
    type(NcDimension) :: dim_x, dim_y, dim_time
    type(NcVariable)  :: var_lon, var_lat

    ! 1.1 create a file
    nc = NcDataset(trim(fname), "w")

    ! create dimensions
    dim_x    = nc%setDimension("x", ncols)
    dim_y    = nc%setDimension("y", nrows)
    dim_time = nc%setDimension("time", -1)

    ! create variables
    var_time = nc%setVariable(vname_time, "i32", (/dim_time/))
    ! var_lat  = nc%setVariable(vname_lat,  "f32", (/dim_x, dim_y/))
    ! var_lon  = nc%setVariable(vname_lon , "f32", (/dim_x, dim_y/))
    var_data = nc%setVariable(vname_data, "f64", (/dim_x, dim_y, dim_time/))

    ! add some variable attributes
    call var_time%setAttribute("units", "days since " // trim(num2str(yearstart - 1, form='(I4)')) // "-12-31 12:00:00")

    ! ! write data of static variables
    ! call var_lat%setData(wlat)
    ! call var_lon%setData(wlon)
    
    ! ! append data within a loop
    ! do i=1, ntime
    !   call var_time%setData(wtime(i),     start=(/i/))
    !   call var_data%setData(wdata(:,:,i), start=(/1,1,i/))
    ! end do

    ! add some more variable attributes
    call var_data%setAttribute("units",   "mm/d")
    call var_data%setAttribute("scaling", 0.1_dp)

    ! add global attributes
    call nc%setAttribute("Author", trim(author_name))



  end subroutine open_netcdf

end module mo_write
