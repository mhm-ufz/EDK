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

  use mo_kind                , only: i4, dp, sp
  use mo_print_message       , only: print_start_message, print_end_message
  use mo_julian              , only: NDAYS, NDYIN, dec2date, julday
  use runControl             , only: flagEDK, interMth,        & ! flag for activate kriging, flag for 'OK' or 'EDK'
    correctNeg,                 & ! pre or temp
    flagVario                     ! flag for activate variogram estimation
  use mainVar                , only: yStart, yEnd, jStart, jEnd, tBuffer, nSta, DEMNcFlag, & ! interpolation time periods
    grid, gridMeteo,            & ! grid properties of input and output grid
    nCell, MetSta, &
    noDataValue
  use kriging                , only: edk_dist, cell
  use mo_setVario            , only: setVario, dMatrix
  use mo_netcdf              , only: NcDataset, NcVariable
  use mo_write               , only: open_netcdf
  use mo_message             , only: message
  use mo_EDK                 , only: EDK
  use mo_ReadData            , only: readData
  use NetCDFVar              , only: invert_y
  USE mo_timer, ONLY : &
    timers_init, timer_start, timer_stop, timer_get              ! Timing of processes
  use mo_string_utils, ONLY : num2str
  !$ use omp_lib, ONLY : OMP_GET_NUM_THREADS           ! OpenMP routines
  implicit none

  character(256)        :: fname
  character(256)        :: author_name
  character(256)        :: vname_data
  integer(i4)           :: i, j, k, t
  integer(i4)           :: icell              ! loop varaible for cells
  integer(i4)           :: jDay               ! loop variable - current julian day
  integer(i4)           :: itimer
  integer(i4)           :: doy                ! day of year
  integer(i4)           :: year, month, day   ! current date
  integer(i4)           :: itime, itemp,sttemp,cnttemp       ! loop variable for chunking write
  integer(i4)           :: jstarttmp, jendtmp ! more loop variables
  integer(i4)           :: n_threads          ! OpenMP number of parallel threads
  integer(i4)           :: ncell_thread
  integer(i4)           :: iThread
  integer(i4)           :: loop_factor
  real(sp), allocatable :: tmp_array(:, :, :) ! temporal array for output
  real(sp), allocatable :: tmp_time(:)        ! temporal array for time output
  real(dp)              :: param(3)           ! variogram parameters
  type(NcDataset)       :: nc_out
  type(NcVariable)      :: nc_data, nc_time
  ! integer(i4), allocatable           :: Nk_old(:)
  ! real(dp), allocatable              :: X(:)

  call print_start_message()

  loop_factor = 10 ! factor for setting openMP loop size

  n_threads = 1
  !$omp PARALLEL
  !$ n_threads = OMP_GET_NUM_THREADS()
  !$omp END PARALLEL
  !$ print *, 'Run with OpenMP with ', trim(num2str(n_threads)), ' threads.'

  ! initialize timers
  call timers_init

  ! start timer for reading
  itimer = 1
  call timer_start(itimer)

  call message('')
  call message(' >>> Reading data')
  call message('')
  call ReadData
  call timer_stop(itimer)
  call message('')
  call message('    in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')
  call message('')
  ! number of cells per thread
  ncell_thread = ceiling(real(nCell, sp) / real(loop_factor * n_threads, sp))
  !print *, 'nCell: ', nCell
  !print *, "ncell_thread: ", ncell_thread
  !print *, 'n_threads: ', n_threads
  ! print *, 'DEMNcFlag', DEMNcFlag

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

  !write(*,*), "jStart = ",jStart

  if (interMth .gt. 0) then
    itimer = 4
    call timer_start(itimer)
    call message(' >>> Perform interpolation')
    call message('')
    ! open netcdf if necessary
    call open_netcdf(nc_out, nc_data, nc_time)

    ! do iCell = 1, nCell
    !   ! initialize cell
    !   allocate(cell(iCell)%Nk_old(nSta))
    !   cell(iCell)%Nk_old = nint(noDataValue)
    ! end do

    if (mod((jEnd - jStart + 1),tBuffer) .eq. 0) then  ! just use mod
      iTime = ((jEnd - jStart + 1)/tBuffer)
    else
      iTime = ((jEnd - jStart + 1)/tBuffer) + 1
    end if
    write(*,*) "Total Number of Threads = ",n_threads
    write(*,*) "Total Number of Time Buffers = ",iTime
    t = 0
    bufferloop: do iTemp = 1, iTime
      write(*,*) "  >>> Started buffer #", iTemp
      jStartTmp = jStart + (iTemp - 1) * tBuffer
      if (iTemp .lt. iTime) then
        jEndTmp = jStartTmp + tBuffer - 1
      else
        jEndTmp = jStartTmp + (jEnd-jStartTmp+1)
      end if   ! use minimum to never exceed jEnd
      jEndTmp = min(jEndTmp, jEnd)

      do iCell = 1, nCell
        ! initialize cell        ! deallocate similarly
        allocate(cell(iCell)%z(jStartTmp:jEndTmp))
        cell(iCell)%z = real(noDataValue, sp)
      end do

      !print *, iTemp, iTime

      !$omp parallel default(shared) &
      !$omp private(iThread, iCell)
      !$omp do SCHEDULE(dynamic)
      do iThread = 1, loop_factor * n_threads
        !  print *, 'thread: ', iThread, " start"

        ncellsloop: do iCell = (iThread - 1) * ncell_thread + 1, min(iThread * ncell_thread, ncell)

          ! check DEM
          if (nint(cell(iCell)%h) == grid%nodata_value ) then
            cell(iCell)%z = gridMeteo%nodata_value
            cycle
          end if
          ! interploation
          call EDK(iCell, jStartTmp, jEndTmp, edk_dist, MetSta, cell(iCell), doOK=(interMth==1))
        end do ncellsloop

        ! print *, 'thread: ', iThread, " end"
      end do
      !$omp end do
      !$omp end parallel

      if (DEMNcFlag == 1) then
        ! write output
        allocate(tmp_array(gridMeteo%nrows, gridMeteo%ncols, jEndTmp - jStartTmp + 1))
        allocate(tmp_time(jEndTmp - jStartTmp + 1))

        k = 0
        if (invert_y) then
          do i = 1, gridMeteo%ncols
            do j = 1, gridMeteo%nrows
              k = k + 1
              tmp_array(j,gridMeteo%ncols - i + 1,:) = cell(k)%z
            end do
          end do
        else
          do i = 1, gridMeteo%ncols
            do j = 1, gridMeteo%nrows
              k = k + 1
              tmp_array(j,i,:) = cell(k)%z
            end do
          end do
        end if

        do i = 1, jEndTmp - jStartTmp + 1
          tmp_time(i) = t
          t = t + 1
        end do

        sttemp = nint(tmp_time(1)+1)
        cnttemp = nint((tmp_time(size(tmp_time)) - sttemp))+2

        !write(*,*),"Final Output ",shape(tmp_array)

        call nc_time%setData(values=tmp_time,start=(/sttemp/),cnt=(/cnttemp/))
        !call nc_data%setData(values=tmp_array,start=(/1,1,sttemp/),cnt=(/size(tmp_array,1),size(tmp_array,2),cnttemp/))
        call nc_data%setData(values=tmp_array,start=(/1,1,sttemp/),cnt=(/size(tmp_array,1),size(tmp_array,2),cnttemp/))

      else
        ! write output
        allocate(tmp_array(gridMeteo%ncols, gridMeteo%nrows, jEndTmp - jStartTmp + 1))
        allocate(tmp_time(jEndTmp - jStartTmp + 1))

        k = 0
        do i = 1, gridMeteo%ncols
          !    do j = 1, gridMeteo%nrows
          !       k = k + 1
          !       tmp_array(i, gridMeteo%nrows - j + 1, :) = cell(k)%z
          !    end do
          ! end do
          if (invert_y) then
            do j = gridMeteo%nrows, 1, -1
              k = k + 1
              tmp_array(i, gridMeteo%nrows - j + 1, :) = cell(k)%z
            end do
          else
            do j = 1, gridMeteo%nrows
              k = k + 1
              tmp_array(i, gridMeteo%nrows - j + 1, :) = cell(k)%z
            end do
          end if
        end do
        !t = 0
        !do i = 1,  jEnd - jStart + 1
        !  tmp_time(i) = t
        !  t = t + 1
        !end do

        do i = 1,  jEndTmp - jStartTmp + 1
          tmp_time(i) = t
          t = t + 1
        end do

        !write(*,*),tmp_time
        sttemp = nint(tmp_time(1)+1)
        cnttemp = nint((tmp_time(size(tmp_time)) - sttemp))+2

        !write(*,*),"Final Output ",shape(tmp_array)

        call nc_time%setData(values=tmp_time,start=(/sttemp/),cnt=(/cnttemp/))
        !call nc_data%setData(values=tmp_array,start=(/1,1,sttemp/),cnt=(/size(tmp_array,1),size(tmp_array,2),cnttemp/))
        call nc_data%setData(values=tmp_array,start=(/1,1,sttemp/),cnt=(/size(tmp_array,1),size(tmp_array,2),cnttemp/))
      end if

      deallocate(tmp_array, tmp_time)
      !deallocate(cell)

      do iCell = 1, nCell
        ! initialize cell
        deallocate(cell(iCell)%z)
        !cell(iCell)%z = noDataValue
      end do

      ! close netcdf if necessary
      !call nc_out%close()   ! outside

    end do bufferloop

    ! close netcdf if necessary
    call nc_out%close()

    call timer_stop(itimer)
    call message('')
    call message('    in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')
    call message('')
  end if
  ! deallocate memory
  call clean
  !
  call print_end_message()
  !
end program ED_Kriging
