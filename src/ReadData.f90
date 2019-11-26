!*************************************************************************
!    SUBROUTINE Read Database Precipitation
!               Read grid DEM
!               Reads parameters variogram
!    AUTHOR:    Luis E. Samaniego-Eguiguren, UFZ
!    UPDATES
!               Created        Sa   21.03.2006
!               Last Update    Sa   11.06.2010
!**************************************************************************
subroutine ReadDataMain()
  use mainVar,    only: noDataValue, cellFactor, DataConvertFactor, &
      yStart, mStart, dStart, yEnd, mEnd, dEnd
  use mo_kind,    only: i4, dp
  use mo_julian , only: NDAYS

  use kriging,    only: maxDist
  use VarFit,     only: vType, nParam, dh, hMax, beta
  use runControl, only: interMth, fnameDEM, DataPathOut, DataPathIn, fNameSta, correctNeg, &
      distZero, flagVario, fNameVario, flagEDK
  use NetCDFVar,  only: FileOut, author_name, variable_name, variable_unit, variable_long_name
  use mo_message, only: message

  implicit none
  !
  integer(i4)                   :: i, ios
  character(256)                :: dummy, fileName
  !
  !===============================================================
  !  namelist definition
  !===============================================================
  !
  namelist/mainVars/noDataValue, DataPathIn, fNameDEM,                    & 
                    DataPathOut, FileOut, fNameSTA, cellFactor, DataConvertFactor, InterMth, correctNeg,  &
                    distZero, author_name, variable_name, variable_unit, variable_long_name,                             & 
                    yStart, mStart, dStart, yEnd, mEnd, dEnd, maxDist, flagVario, vType, nParam,  &  
                    fNameVario, dh, hMax
  !
  ! -----------------------------------------------------------------------
  !	                               MAIN.DAT
  ! -----------------------------------------------------------------------
  open(unit=10, file='edk.nml', STATUS='OLD', ACTION='read')
  read(10, nml=mainVars)
  close(10)

  ! -----------------------------------------------------------------------
  !                         consistency checks
  ! -----------------------------------------------------------------------
  if ((interMth .gt. 2) .or. (interMth .lt. 0)) then
    call message(' >>> Interpolated Method Value not valid')
    call message('     Please choose one of the following:')
    call message('     0 - no interpolation')
    call message('     1 - ordinary kriging')
    call message('     2 - external drift kriging')
    call message('')
    stop 1
  end if

  
  FileOut = trim(DataPathOut) // trim(FileOut)
  !
  !	*************************************
  ! Read initial control parameters
  !	*************************************

  if (flagVario .AND. flagEDK) then
     print*, '***Warning: Both flags flagVario and flagEDK should not be activated at the same time!'
     ! stop
  end if
  
  ios = 0
  allocate (beta(nParam))

  if ( .not. flagVario ) then
    filename = trim(fNameVARIO)
    open (unit=20, file=filename, STATUS='OLD', ACTION='READ')
    read (20,1) dummy
    read (20, *,iostat=ios) (beta(i), i=1,nParam), i
    close (20)
  end if

  ! if vario is read from file take read in vario type
  if ( .not. flagVario ) then
     vType = i
  end if
  if (.NOT. ((vType .EQ. 1) .OR. (vType .EQ. 2))) then
     print*, "***ERROR: Variogram type not known: ", vType
  end if

  write(*,*) 'variogram parameters: ' 
  write(*,'(3(a12),   a18)') 'nugget', 'sill', 'range', 'Variogramtype' 
  write(*,'(3(f12.2), i18)') (beta(i), i=1,nParam), vType
  
  1 format (a256)
  
end subroutine ReadDataMain


subroutine ReadStationLut
  use mo_kind,    only: i4, dp

  use mainVar,    only: nSta, MetSta, grid
  use kriging,    only: maxDist
  use VarFit,     only: hmax
  use runControl

  implicit none
  !
  integer(i4)                         :: iSta
  integer(i4)                         :: n_stations_all  ! number of stations in file
  character(256)                      :: dummy, filename
  real(dp)                            :: search_distance ! highest distance for considered stations
  real(dp)                            :: domain_bottom   ! upper border of station search domain
  real(dp)                            :: domain_top      ! lower border of station search domain
  real(dp)                            :: domain_left     ! left border of station search domain
  real(dp)                            :: domain_right    ! right border of station search domain
  real(dp), dimension(:), allocatable :: temp_id         ! ID of read in station
  real(dp), dimension(:), allocatable :: temp_x          ! easting coordinate of station
  real(dp), dimension(:), allocatable :: temp_y          ! northing coordinate of station
  real(dp), dimension(:), allocatable :: temp_h          ! heigth of station
  
  !-----------------------------------------------------------------------
  !	                   STATIONS COORDINATES AND ELEVATION
  !-----------------------------------------------------------------------

  ! take only stations which are in search distance
  ! two possibilties for search distance for vario - hMax for EDK - maxDist
  if (flagVario) then
     search_distance = hMax
  else
     search_distance = maxDist
  end if

  filename = trim(fNameSTA)
  open (unit=20, file=filename, status='old', action='read')
  read (20,*)  dummy, n_stations_all
  read (20, 1) dummy

  ! derive borders for station search
  domain_bottom = -search_distance + grid%yllcorner  
  domain_top    =  search_distance + grid%yllcorner + grid%nrows * grid%cellsize
  domain_left   = -search_distance + grid%xllcorner
  domain_right  =  search_distance + grid%xllcorner + grid%ncols * grid%cellsize 

  ! allocate
  allocate(temp_id(n_stations_all)); allocate(temp_x(n_stations_all))
  allocate( temp_y(n_stations_all)); allocate(temp_h(n_stations_all))

  ! initialize
  temp_id = grid%nodata_value; temp_h = grid%nodata_value
  temp_x = grid%nodata_value; temp_y = grid%nodata_value
  ! init
  nSta = 1
 
  ! new format (reads only part of the line)
  do  iSta = 1 , n_stations_all
     read (20, *) temp_id(nSta), temp_x(nSta), temp_y(nSta), temp_h(nSta)

     ! check if station is within search distance
     ! if yes increase station counting by one
     ! otherwise overwrite station in next iteration
     if ((temp_x(nSta) .GE. domain_left)   .AND. (temp_x(nSta) .LE. domain_right) .AND. &
         (temp_y(nSta) .GE. domain_bottom) .AND. (temp_y(nSta) .LE. domain_top)  ) then
        nSta            = nSta + 1
     end if
  end do

  nSta = nSta - 1
  
  close(20)

  ! save relevant stations in MetSta
  allocate(MetSta(nSta))
  do iSta = 1, nSta
     MetSta(iSta)%Id =  temp_id(iSta)
     MetSta(iSta)%x  =  temp_x(iSta)
     MetSta(iSta)%y  =  temp_y(iSta)
     MetSta(iSta)%h  =  temp_h(iSta)
  end do

  deallocate(temp_id, temp_x, temp_y, temp_h)
  
  print*, "Only " , nSta , ' of the given ' , n_stations_all, ' are considered since they are within the search distance.'
  
  ! ! test filtering
  ! open(20, file='station_test.txt', status='unknown', action='write')
  ! write(20,*) dummy
  ! do  iSta = 1 , nSta
  !    write(20,"(i10,' ',f12.2,' ',f12.2,' ',f12.2)")  MetSta(iSta)%Id, MetSta(iSta)%x, MetSta(iSta)%y, MetSta(iSta)%h
  ! end do
  ! close(20)
  ! !pause
  
  ! formats
  1   format (a256)
  
end subroutine ReadStationLut

! -----------------------------------------------------------------------
!	                    DEM Blocks (all must be EQUAL in size!)
! -----------------------------------------------------------------------
subroutine ReadDEM
  use mainVar
  use mo_kind, only        : i4
  use runControl

  implicit none

  integer(i4)             :: i, j
  character(256)          :: dummy, fileName
  !

  fileName = trim(fNameDEM)
  call readHeader(15, fileName)
  if (.not. allocated(G) ) allocate (G(grid%nrows , grid%ncols))
  do i=1,grid%nrows
   read (15, *) (G(i,j)%h, j=1,grid%ncols) 
  end do

  close (15)
  ! Assumed all blocks equal
  nCell = gridMeteo%ncols * gridMeteo%nrows
  !
  print *, 'DEM read ...'
  ! save gridMeteo*
  fileName= trim(dataPathOut)//trim('header.txt')
  call writeHeader(20, fileName)
  !
  ! formats
  100 format (i1.1, a4)
  101 format (i2.2, a4)

end subroutine ReadDEM

!
!
!	*************************************************************************
!
! SUBROUTINE   READ HEADER GRIDs
!
!	*************************************************************************
subroutine readHeader(ic, fName)
  use mainVar
  use mo_kind, only    : i4, dp
  implicit none
  !
  character(50)                  :: dummy
  integer(i4), intent(in)        :: ic                  ! input channel 
  character(256), intent(in)     :: fName
  !
  open (unit=ic, file=fName, status='old', action='read')
  read (ic, *) dummy, grid%ncols                       ! number of columns
  read (ic, *) dummy, grid%nrows                       ! number of rows
  read (ic, *) dummy, grid%xllcorner                   ! x coordinate of the lowerleft corner
  read (ic, *) dummy, grid%yllcorner                   ! y coordinate of the lowerleft corner
  read (ic, *) dummy, grid%cellsize                    ! cellsize x = cellsize y
  read (ic, *) dummy, grid%nodata_value                ! code to define the mask
  !
  if (mod(float(grid%ncols),float(cellFactor)) /= 0.0_dp .or.  &
      mod(float(grid%nrows),float(cellFactor)) /= 0.0_dp     ) then
      print *, 'nrows or ncols are not exactly divisible by cellFactor ...'
      stop
  end if
  !
  gridMeteo%ncols = ceiling(float(grid%ncols)/float(cellFactor))
  gridMeteo%nrows = ceiling(float(grid%nrows)/float(cellFactor))
  !
  gridMeteo%cellsize     = grid%cellsize * cellFactor
  gridMeteo%xllcorner    = grid%xllcorner 
  ! to correct the increment that may be produced by ceiling (ensure consistent output in surfer)
  gridMeteo%yllcorner    = grid%yllcorner  + float(grid%nrows) * float(grid%cellsize) &
                                           - float(gridMeteo%nrows) * float(gridMeteo%cellsize)  
  gridMeteo%nodata_value = grid%nodata_value
  thresholdDist = 7.07d-1 * gridMeteo%cellsize
end subroutine readHeader


!	**************************************************************************
!
! SUBROUTINE      READ METEOROLOGICAL DATABASE
! UPDATES         Sa.2010-07-06     reading formats / conversion units
!
!	**************************************************************************
subroutine ReadDataMeteo
  use mo_kind, only         : i4, dp
  use mo_julian, only       : NDAYS 
  use mainVar
  use runControl
  implicit none
  integer(i4)               :: i, ios
  integer(i4)               :: jDay, jS, jE
  integer(i4)               :: ds, ms, ys
  integer(i4)               :: de, me, ye
  character(256)            :: dummy
  character(12)             :: dateStart
  character(12)             :: dateEnd
  character(256)            :: fileName
  integer(i4)               :: nStaEfec
  real(dp)                  :: p
  !
  ! allocate vectors for interpolation period
  jStart = NDAYS (dStart, mStart, yStart)
  jEnd   = NDAYS (dEnd,   mEnd,   yEnd)

  nStaEfec = 0
  do i = 1, nSta
    allocate (MetSta(i)%z(jStart:jEnd))
    MetSta(i)%z = noDataValue
    !
    ! read yearly data file
    write (dummy, 2) MetSta(i)%Id
    fileName = trim(dataPathIn)//trim(dummy)
    print *, 'read file: '//trim(fileName)
    open (60, file=fileName, status='old', action='read', iostat=ios)
    read (60, *) dummy
    !
    !>>>  -------------------------------------------------------------
          ! this is temporary / better, robust free format MUST be used
          !read(60,1) dummy
          !k = scan(dummy, '/', back=.true.)
          read(60,*) dummy, dateStart , dummy, dateEnd
          dateStart = dateStart(3:12)
          dateEnd   = dateEnd(3:12)
          read(dateStart, 3) ys, ms, ds
          read(dateEnd, 3)   ye, me, de
    !<<<  -------------------------------------------------------------   
    !    
    ! set limits
    jS = NDAYS (ds,ms,ys)
    jE = NDAYS (de,me,ye)
    jDay = jS - 1
    ! check whether the time series is within the period of interest
    if ( (jE >= jStart ) .and.  (jS <= jEnd ) ) then
        do while (.NOT. ios < 0)
          read (60, *,iostat=ios) p
          jDay = jDay + 1
          ! check limits
          if (jDay < jStart ) cycle
          if (jDay > jEnd   ) cycle
          MetSta(i)%z(jDay) = p
       end do
       jDay = jDay - 1 ! because iostat is reading an additional line at the end
       close (60)
       nStaEfec = nStaEfec + 1
       ! unit conversion
       if ( DataConvertFactor /= 1.0_dp ) then
          where ( MetSta(i)%z(:) /= noDataValue  )   MetSta(i)%z(:) = MetSta(i)%z(:) * DataConvertFactor
       end if
       ! check consistency
       if (jDay /= jE ) then
          print *, 'File ', trim(fileName), ' has missing values.',  jE - jDay 
          stop
       end if
    else
       close (60)
    end if
  end do
  print *, nStaEfec, ' of ', nSta, ' have data from ', yStart , ' to ', yEnd
  !
  ! formats
  1   format (a256)
  2   format (i5.5,'.dat')
  3   format (i4, 1x, i2, 1x, i2)
  !
end subroutine ReadDataMeteo

