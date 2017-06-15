!> \file mo_write_fluxes_states.f90

!> \brief Creates NetCDF output for different fluxes and state variabels of mHM.

!> \details NetCDF is first initialized and later on variables are put to the NetCDF.

!> \authors Matthias Zink
!> \date Apr 2013



MODULE mo_write_fluxes_states

  ! This module creates the output for mHM.

  ! Written Matthias Zink, Apr 2012

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: WriteFluxStateInit
  PUBLIC :: WriteFluxState
  PUBLIC :: CloseFluxState_file

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !      NAME
  !          WriteFluxStateInit

  !>        \brief Initialization of the NetCDF for mHM outputs.

  !>        \details The NetCDF file is set up with its dimensions, variables and variable attributes. Additionally 
  !>                 the arrays for aggregating the output to the output time step are allocated.

  !     INTENT(IN)
  !          None

  !     INTENT(INOUT)
  !>         \param[inout] "real(dp), allocatable :: data_out(:)"        ! Interception

  !     INTENT(OUT)
  !>         \param[out] "integer(i4), intent(out) :: ncid" ! ID of NetCDF to be written in

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !

  !     EXAMPLE
  !        

  !     LITERATURE

  !     HISTORY
  !>        \author Matthias Zink
  !>        \date May 2014


  Subroutine WriteFluxStateInit(ncid)

    use RunControl,            only: DataPathOut            ! output directory
    use mainVar,               only: gridMeteo,           & ! grid properties
                                     yStart, mStart, dStart ! year, month, day of first day of interpolation
    use NetCDFVar,             only: variable_name,       & ! attributes of variable: name
                                     variable_unit,       & ! attributes of variable: unit
                                     variable_long_name,  & ! Attributes of variable: long name
                                     yCoor, xCoor,        & ! kartesian coordinates
                                     lons, lats             ! geographic coordinates
    use mo_julian,             only: NDYIN
    use mo_ncwrite,            only: create_netCDF, write_static_netCDF, V, GAtt
    use mo_set_netcdf_outputs, only: set_netCDF
    ! use mo_message,            only: message
    !   
    implicit none
    
    !
    integer(i4), intent(out)                         :: ncid      ! ID of NetCDF
    !
    !
    ! local
    character(8)                                     :: date
    character(10)                                    :: time
    character(256)                                   :: fName, dummy

    ! Initialize NetCDF dimensions and Variables incl. attributes
    call set_netCDF(1, gridMeteo%ncols, gridMeteo%nrows)

    ! attributes of the output variable
    ! name
    V(6)%name          = trim(variable_name)
    !unit
    V(6)%att(1)%values = trim(variable_unit)
    ! long_name
    V(6)%att(2)%values = trim(variable_long_name)
    !
    !*******************************************************
    ! here comes the attributes which are equal for all vars
    ! up to Var 5 are static Vars
    !*******************************************************
    write(dummy,*) gridMeteo%nodata_value
    ! scale_factor
    V(6)%att(3)%values = "1."
    ! _FillValue
    V(6)%att(4)%values = trim(adjustl(dummy))
    ! missing_value
    V(6)%att(5)%values = trim(adjustl(dummy))
    ! coordinates
    V(6)%att(6)%values = "lat lon"

    !
    !*******************************************************
    ! set Attributes for dimensions
    ! set unit for dimension time 
    ! (e.g. hours since 2008-01-01 00:00:00)
    !*******************************************************
    write(dummy,"('days since ', i4, '-' ,i2.2, '-', i2.2, 1x, '00:00:00')") yStart, mStart, dStart
    V(3)%att(1)%values = trim(dummy)

    !
    ! global attributes
    GAtt(1)%name    = "title"
    GAtt(1)%values  = "External drift Kriging of DWD data"
    !
    GAtt(2)%name    = "creating_date"
    call DATE_AND_TIME(DATE=date, TIME=time)
    dummy = trim(date) // trim(time)
    write(dummy,"(a4,'-',a2,'-',a2,1x,a2,':',a2,':',a2)") date(1:4), &
      date(5:6), date(7:8), time(1:2), time(3:4), time(5:6)
    GAtt(2)%values   = trim(dummy)
    !
    GAtt(3)%name    = "author"
    GAtt(3)%values   = trim("Matthias Zink")
    !
    GAtt(4)%name    = "institution"
    GAtt(4)%values   = trim("Helmholtz Centre for Environmental Research - UFZ,") // &
                       trim("Department of Computational Hydrosystems, Stochastic Hydrology Group")
    ! to create a new netCDF
    fName = trim(DataPathOut) // trim(variable_name) // '.nc'
    !
    ! call message('')
    ! call message('  OUTPUT: Writing NetCDF file in')
    ! call message('     to ', trim(fName))
    !
    call create_netCDF(trim(fName), ncid)   

    !
    !*******************************************************
    ! write static variables  
    ! put coordinates sytem to the NetCDF
    !*******************************************************
    call CoordSytem(gridMeteo%xllcorner, gridMeteo%yllcorner, gridMeteo%cellsize, gridMeteo%ncols, gridMeteo%nrows)
    !
    V(1)%G1_d        => xCoor
    V(2)%G1_d        => yCoor
    V(4)%G2_d        => lons
    V(5)%G2_d        => lats
    call write_static_netCDF(ncid) 
    !
  end Subroutine WriteFluxStateInit


  Subroutine WriteFluxState(curr_wstep, ncid, data_out) 
                             
    use mainVar,    only: gridMeteo     ! Number of horizons to model
    use mo_ncwrite, only: write_dynamic_netCDF, V

    implicit none
    !
    ! local
    integer(i4),               intent(in)           :: curr_wstep ! current write out step
    integer(i4),               intent(in)           :: ncid       ! ID of NetCDF
    ! States L1
    real(dp), dimension(:,:),  intent(in)           :: data_out   ! output data

    !
    integer(i4),                           target   :: tstep      ! time step of output
    real(dp), dimension(:,:), allocatable, target   :: OutPut
    !
    
    ! Write files
    allocate( OutPut(size(data_out, dim=1), size(data_out, dim=2)) ); OutPut = real(gridMeteo%nodata_value, dp)

    OutPut(:,:)  = data_out
    V(6)%G2_d =>  OutPut(:,:)

    ! timestep
    ! "-1" to start at timestep 0
    tstep = curr_wstep-1
    V(3)%G0_i => tstep
    call write_dynamic_netCDF(ncid, curr_wstep)

    nullify ( V(3)%G0_i)
    nullify ( V(6)%G2_d)

    deallocate(OutPut)
    !
  end subroutine WriteFluxState
  
  ! ------------------------------------------------------------------

  !      NAME
  !          CloseFluxState_file

  !>        \brief Close the netcdf file containing flux and states

  !>        \details  Close the netcdf file containing flux and states

  !     INDENT(IN)
  !>        \param[in] "integer(i4)    :: iBasin"   - ID number of basin
  !>        \param[in] "integer(i4)    :: ncid_out" - ID of NetCDF file

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !

  !     EXAMPLE
  !        

  !     LITERATURE

  !     HISTORY
  !>        \author Rohini Kumar & Stephan Thober
  !>        \date August 2013
  !         Modified, Matthias Zink May 2014

 subroutine CloseFluxState_file(ncid_out)

    use mo_ncwrite,          only: close_netCDF

    implicit none

    integer(i4), intent(in)                          :: ncid_out  ! ID of NetCDF

    !write(*,*) ''
    !write(*,*) '  OUTPUT: saved netCDF file.'
    call close_netCDF(ncid_out)
    !
 end subroutine CloseFluxState_file

  ! ------------------------------------------------------------------

  !      NAME
  !          CoordSytem

  !         \brief Setting Gauss Krueger 3 and WGS84 coordinates for output

  !         \details Calculating the coordinates of the projected coordinate system with
  !                  the help of the coordinates of the lower left corner, the cellsize and the number of
  !                  columns and rows. Extracting lats and lons from respective global variables

  !     INTENT(IN)
  !       integer(i4), intent(in) :: iBasin        ! basin index

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !

  !     EXAMPLE
  !        

  !     LITERATURE

  !     HISTORY
  !         \author Matthias Zink
  !         \date Apr 2013

  subroutine CoordSytem(xll, yll, cz, nrows, ncols)
    !
    use NetCDFVar, only : & 
         yCoor, xCoor       , & ! kartesian coordinates
         lons, lats             ! geographic coordinates

    !
    implicit none

    real(dp), intent(in)                :: xll, yll      ! coordinates of the lower left corner of
    !  projected coordinate system
    integer(i4), intent(in)             :: cz            ! cellsize of the  projected coordinate system
    integer(i4) , intent(in)            :: nrows, ncols  ! number row and columns of array V accrording
    !
    real(dp)                :: xc, yc
    integer(i4)             :: i, j

    !
    if (allocated(xCoor) ) deallocate ( xCoor ) ! northing
    if (allocated(yCoor) ) deallocate ( yCoor ) ! easting
    if (allocated(lons))   deallocate ( lons  )
    if (allocated(lats))   deallocate ( lats  )
    allocate ( xCoor(nrows)       )
    allocate ( yCoor(ncols)       ) 
    allocate ( lons(nrows, ncols) )
    allocate ( lats(nrows, ncols) )
    !
    ! def northings and eastings arrays
    xCoor(1)     =  xll + 0.5_dp * cz
    do i = 2, nrows
       xCoor(i)   =  xCoor(i-1) + cz
    end do
    ! inverse for Panoply, ncview display
    yCoor(ncols) =  yll + 0.5_dp * cz
    do j = ncols-1,1,-1 
       yCoor(j)   =  yCoor(j+1) + cz
    end do
    ! find lat and lon (curvilinear orthogonal gridding, See Panoply ref) 
    do i = 1, nrows 
       do j = 1,ncols   
          ! (x,y) -> (lon, lat) 
          call CoorTransInv(xCoor(i), yCoor(j), xc, yc) 
          lons(i,j) = xc 
          lats(i,j) = yc 
       end do
    end do
    !
  end subroutine CoordSytem

END MODULE mo_write_fluxes_states




