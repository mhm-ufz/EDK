!****************************************************************************
!
!  PROGRAM:  EDK
!
!  PURPOSE:  Perform EDK (daily)
!            Store results suitable for mHM
!            special version to interpolate values in blocks
!            *** can be parallelized ***
!
!  NOTE:     CVF do not use RECL in bytes by default.
!            use:
!                in Project settings | Fortran | Fortran data |
!                mark use Bytes as RECL = unit for unformated files.
!  UPDATES
!            Created        Sa   21.03.2006
!            Last Update    Sa   11.06.2010    ! blocks, whole Germany
!            Last Update    Zi   04.02.2012    ! changed to general edk version 
!                                              ! (excluded block seperation)
!****************************************************************************
program ED_Kriging

  use mo_kind                , only: i4, dp
  use mo_julian              , only: NDAYS, NDYIN
  use mo_write_fluxes_states , only: WriteFluxStateInit, WriteFluxState, CloseFluxState_file

  use runControl             , only: outputformat,               & ! outputformat either 'nc' or 'bin'
                                     flagEDK, flagMthTyp,        & ! flag for activate kriging, flag for 'OK' or 'EDK'
                                     flagVarTyp,                 & ! pre or temp
                                     flagVario                     ! flag for activate variogram estimation
  use mainVar                , only: yStart, yEnd, jStart, jEnd, & ! interpolation time periods
                                     grid, gridMeteo,            & ! grid properties of input and output grid
                                     nCell                         ! number of cells
  use kriging
 
  implicit none

  integer(i4)                           :: icell             ! loop varaible for cells
  integer(i4)                           :: jDay              ! loop variable - current julian day
  integer(i4)                           :: doy               ! day of year
  integer(i4)                           :: year, month, day  ! current date
  integer(i4)                           :: netcdfid          ! id of netcdf files
  real(dp), dimension(:,:), allocatable :: tmp_array         ! temporal array for output

  call Timer
  call ReadDataMain
  
  ! read DEM
  call ReadDEM

  ! read look up table for meteorological stations
  call ReadStationLut
  
  call dMatrix
  print*, 'calculated distance matrix'
  ! read whole METEO data
  call ReadDataMeteo
  print*, 'finished reading of meteorological data'
  ! estimate variogram
  call setVario
  ! write variogram  
  if (flagVario) call WriteDataMeteo(0,0,2)
  ! 
  if (flagEDK) then
     ! open netcdf if necessary
     if (outputformat=='nc') call WriteFluxStateInit(netcdfid)

     timeloop: do jday = jStart, jEnd

        call NDYIN(jday, day, month, year)
        doy = jday - NDAYS(1,1,year) + 1

        print *, 'YEAR: ',year, 'DOY: ', doy

        ncellsloop: do iCell = 1, nCell
           ! interploation
           select case (flagMthTyp)
           case (1)
             call EDK(jday,iCell)
             
           case (2)
              call OK(jday,iCell)
           end select
         end do ncellsloop

        ! correct precipitation values
        if (flagVarTyp == 1) then
           where ((cell(:)%z .LT. 0.0_sp) .AND. (cell(:)%z .GT. real(grid%nodata_value, sp)) )
              cell(:)%z = 0.0_sp
           end where
        end if

        ! write output
        select case (trim(adjustl(outputformat)))

        case('bin')
           call WriteDataMeteo(year,doy,1)

        case('nc')
           allocate(tmp_array(gridMeteo%nrows, gridMeteo%ncols)); tmp_array=real(grid%nodata_value, dp)
           tmp_array = real(reshape(cell(:)%z,(/gridMeteo%nrows, gridMeteo%ncols/)), dp)
           call WriteFluxState((jday-jStart+1), netcdfid, transpose(tmp_array))
           deallocate(tmp_array)

        case DEFAULT
           print*, '***ERROR: Output format not known: ', trim(adjustl(outputformat))
           stop
        end select

      end do timeloop

      ! close netcdf if necessary
      if (outputformat=='nc') call CloseFluxState_file(netcdfid)

  end if
  ! deallocate memory !
  call clean
  ! Timer
  call Timer
  
  ! very important for check cases 
  write(*,*) 'Kriging finished!'
  !
end program

