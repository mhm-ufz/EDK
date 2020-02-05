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
  use mainVar                , only: yStart, yEnd, jStart, jEnd, tBuffer, nSta, & ! interpolation time periods
      grid, gridMeteo,            & ! grid properties of input and output grid
      nCell, MetSta, &
      noDataValue
  use kriging                , only: dCS, dS, cell
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
  !$ integer(i4)        :: n_threads          ! OpenMP number of parallel threads
  real(sp), allocatable :: tmp_array(:, :, :) ! temporal array for output
  real(sp), allocatable :: tmp_time(:)        ! temporal array for time output
  real(dp)              :: param(3)           ! variogram parameters
  type(NcDataset)       :: nc_out
  type(NcVariable)      :: nc_data, nc_time
  integer(i4), allocatable           :: Nk_old(:)
  real(dp), allocatable              :: X(:)
   
  call print_start_message()
  
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
  call ReadData
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

!write(*,*), "jStart = ",jStart

  if (interMth .gt. 0) then
    itimer = 4
    call timer_start(itimer)
    call message(' >>> Perform interpolation')
    call message('')
    ! open netcdf if necessary
    call open_netcdf(nc_out, nc_data, nc_time)

   ! do iCell = 1, nCell
    !  ! initialize cell
    !  allocate(cell(iCell)%z(jStart:jEnd))
    !  cell(iCell)%z = noDataValue
    !end do

! Akash---------------------------------------------------------------- 
  

    do iCell = 1, nCell
      ! initialize cell
      allocate(cell(iCell)%Nk_old(nSta))
      cell(iCell)%Nk_old = -9999
    end do


  if (mod((jEnd - jStart + 1),tBuffer) .eq. 0) then  ! just use mod 
        iTime = ((jEnd - jStart + 1)/tBuffer)                         
  else   
  iTime = ((jEnd - jStart + 1)/tBuffer) + 1
  end if
  write(*,*),"Total Number of Time Buffers = ",iTime
  t = 0 
  bufferloop: do iTemp = 1, iTime
    !call message('Time Loop Running')
    
    !call open_netcdf(nc_out, nc_data, nc_time) ! dont do this 
     
    !jStartTmp = jStart
    !if (iTemp .lt. iTime) then
    !    jEndTmp = jStartTmp + 100
    !else
    !    jEndTmp = (jEnd-jStart+1) - ((iTemp-1) * 100)
    !end if   ! use minimum to never exceed jEnd
     
    jStartTmp = jStart + (iTemp - 1) * tBuffer
    if (iTemp .lt. iTime) then
        jEndTmp = jStartTmp + tBuffer - 1
    else 
        jEndTmp = jStartTmp + (jEnd-jStartTmp+1) 
    end if   ! use minimum to never exceed jEnd
        jEndTmp = min(jEndTmp, jEnd) 
   
  ! write(*,*), "jStartTmp = ",jStartTmp," jEndTmp = ",jEndTmp,"jEnd = ",jEnd

   
    do iCell = 1, nCell
      ! initialize cell        ! deallocate similarly 
      allocate(cell(iCell)%z(jStartTmp:jEndTmp))
      cell(iCell)%z = noDataValue
    end do



   
!$OMP parallel default(shared) &
!$OMP private(iCell, X, Nk_old)
!$OMP do SCHEDULE(STATIC)
    ncellsloop: do iCell = 1, nCell

      ! check DEM
      if (nint(cell(iCell)%h) == grid%nodata_value ) then
        cell(iCell)%z = gridMeteo%nodata_value
        cycle
      end if
      !write(*,*),"Interpolation Method = ",interMth
      !write(*,*),"Cell = ",iCell
      !write(*,*),"jStartTmp = ",jStartTmp
      !write(*,*),"jEndTmp = ",jEndTmp 
      !write(*,*),"dCS = ",dCS
      !write(*,*),"MetSta = ",MetSta
      !write(*,*),"Flag 1"
    
      ! interploation
      select case (interMth)
      case (1)
        !call EDK(iCell, jStartTmp, jEndTmp, dCS, MetSta, dS, cell, doOK=.True.)
        !call EDK(iCell, jStartTmp, jEndTmp, dCS, MetSta, dS, cell, tempX, tempNkOld, doOK=.True.)
         call EDK(iCell, jStartTmp, jEndTmp, dCS, MetSta, dS, cell, cell(iCell)%W, cell(iCell)%Nk_old, doOK=.True.)
      case (2)
        !call EDK(iCell, jStartTmp, jEndTmp, dCS, MetSta, dS, cell)
        !call EDK(iCell, jStartTmp, jEndTmp, dCS, MetSta, dS, cell, tempX, tempNkOld)
        call EDK(iCell, jStartTmp, jEndTmp, dCS, MetSta, dS, cell, cell(iCell)%W, cell(iCell)%Nk_old)
      end select
      
    !write(*,*),"X after call EDK = ",X
    end do ncellsloop
    !$OMP end do
    !$OMP end parallel
      
   

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
 
   !write(*,*),"array_size = ",size(tmp_array,1)," ",size(tmp_array,2) 
   !write(*,*),"array_size = ",size(tmp_array,3)
   !write(*,*),"tmp_time_1 = ",sttemp
   !write(*,*),"tmp_last = ",cnttemp
   !write(*,*),"Start tmp_time",(/sttemp/)
   !write(*,*),"End tmp_time",(/cnttemp/)
   !write(*,*),"Start tmp_array",(/sttemp,1,1/)
   !write(*,*),"End tmp_array",(/cnttemp,size(tmp_array,2),size(tmp_array,1)/)     
                                 
    !call nc_time%setData(values=tmp_time,start=(/sttemp/),cnt=(/cnttemp/))
    !call nc_data%setData(values=tmp_array,start=(/sttemp,1,1/),cnt=(/cnttemp,size(tmp_array,2),size(tmp_array,1)/))
    call nc_time%setData(values=tmp_time,start=(/sttemp/),cnt=(/cnttemp/))
    call nc_data%setData(values=tmp_array,start=(/1,1,sttemp/),cnt=(/size(tmp_array,1),size(tmp_array,2),cnttemp/))
 

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
! Akash end---------------------------------------------------- 

     do iCell = 1, nCell
      ! initialize cell
      deallocate(cell(iCell)%Nk_old)
      !cell(iCell)%z = noDataValue
    end do

    do iCell = 1, nCell
      ! initialize cell
      if (allocated(cell(iCell)%W)) deallocate(cell(iCell)%W)
      !cell(iCell)%z = noDataValue
    end do

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

