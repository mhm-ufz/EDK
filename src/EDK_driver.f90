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
  use runControl             , only: flagEDK, interMth,        & ! flag for activate kriging, flag for 'OK' or 'EDK'
                                     correctNeg,                 & ! pre or temp
                                     flagVario                     ! flag for activate variogram estimation
  use mainVar                , only: yStart, yEnd, jStart, jEnd, & ! interpolation time periods
                                     grid, gridMeteo,            & ! grid properties of input and output grid
                                     nCell, MetSta
  use kriging                , only: dCS, dS
  use mo_setVario            , only: setVario, dMatrix
  use mo_netcdf              , only: NcDataset, NcVariable
  use mo_write               , only: open_netcdf
  use mo_message             , only: message
  use kriging
  use mo_EDK                 , only: EDK
  USE mo_timer, ONLY : &
          timers_init, timer_start, timer_stop, timer_get              ! Timing of processes
  !$ use mo_string_utils, ONLY : num2str
  !$ use omp_lib, ONLY : OMP_GET_NUM_THREADS           ! OpenMP routines
  
  implicit none

  character(256)                        :: fname
  character(256)                        :: author_name
  character(256)                        :: vname_data
  integer(i4)                           :: icell             ! loop varaible for cells
  integer(i4)                           :: jDay              ! loop variable - current julian day
  integer(i4)                           :: itimer
  integer(i4)                           :: doy               ! day of year
  integer(i4)                           :: year, month, day  ! current date
  !$ integer(i4)                        :: n_threads        ! OpenMP number of parallel threads
  real(dp), dimension(:,:), allocatable :: tmp_array         ! temporal array for output
  real(dp)                              :: param(3)          ! variogram parameters
  type(NcDataset)                       :: nc_out
  type(NcVariable)                      :: nc_data, nc_time

  !$OMP PARALLEL
  !$ n_threads = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  !$ print *, 'Run with OpenMP with ', trim(num2str(n_threads)), ' threads.'

  ! initialize timers
  call timers_init

  ! start timer for reading
  itimer = 1
  call timer_start(itimer)

  call message('')  
  call message(' >>> Reading data')
  call message('')
  call ReadDataMain
  
  ! read DEM
  call ReadDEM

  ! read look up table for meteorological stations
  call ReadStationLut
  
  ! read whole METEO data
  call ReadDataMeteo
  call timer_stop(itimer)
  call message('')
  call message('    in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')
  call message('')

  itimer = 2
  call timer_start(itimer)
  call message(' >>> Calculating distance matrix')
  call message('')
  ! call distance matrix
  call dMatrix
  call timer_stop(itimer)
  call message('')
  call message('    in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')
  call message('')


  itimer = 3
  call timer_start(itimer)
  call message(' >>> Estimate variogram')
  call message('')
  ! estimate variogram
  call setVario(param)
  ! write variogram  
  if (flagVario) call WriteDataMeteo(0,0,2)
  call timer_stop(itimer)
  call message('')
  call message('    in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')
  call message('')


  
  if (interMth .gt. 0) then
    itimer = 4
    call timer_start(itimer)
    call message(' >>> Perform interpolation')
    call message('')
    ! open netcdf if necessary
    call open_netcdf(nc_out, nc_data, nc_time)

    timeloop: do jday = jStart, jEnd

      call NDYIN(jday, day, month, year)
      doy = jday - NDAYS(1,1,year) + 1

      print *, 'YEAR: ',year, 'DOY: ', doy

      
      !$OMP parallel default(shared) &
      !$OMP private(iCell)
      !$OMP do SCHEDULE(STATIC)
      ncellsloop: do iCell = 1, nCell
        ! check DEM
        if (nint(cell(iCell)%h) == grid%nodata_value ) then
          cell(iCell)%z = gridMeteo%nodata_value
          cycle
        end if
        ! interploation
        select case (interMth)
        case (1)
          call OK(jday, iCell)
        case (2)
          call EDK(jday, iCell, dCS, MetSta, dS, cell)
        end select
      end do ncellsloop
      !$OMP end do
      !$OMP end parallel

      ! correct negative values
      if (correctNeg) then
        where ((cell(:)%z .LT. 0.0_sp) .AND. (cell(:)%z .GT. real(grid%nodata_value, sp)) )
          cell(:)%z = 0.0_sp
        end where
      end if

      ! write output
      allocate(tmp_array(gridMeteo%nrows, gridMeteo%ncols)); tmp_array=real(grid%nodata_value, dp)
      tmp_array = real(reshape(cell(:)%z,(/gridMeteo%nrows, gridMeteo%ncols/)), dp)
      call nc_time%setData(jday-jStart + 1,      start=(/jday - jStart + 1/))
      call nc_data%setData(tmp_array, start=(/1, 1, jday - jStart + 1/))

      deallocate(tmp_array)

    end do timeloop

    ! close netcdf if necessary
    call nc_out%close()

    call timer_stop(itimer)
    call message('')
    call message('    in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')
    call message('')
  end if
  ! deallocate memory
  call clean
  
  ! very important for check cases 
  write(*,*) 'Kriging finished!'
  !
end program

