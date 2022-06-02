module mo_edk_read_data

  use mo_kind, only: dp

  implicit none

  private

  public :: ReadData

  type extend
    real(dp) :: left
    real(dp) :: top
    real(dp) :: right
    real(dp) :: bottom
  end type extend

contains
  !*************************************************************************
  !    SUBROUTINE Read Database Precipitation
  !               Read grid DEM
  !               Reads parameters variogram
  !    AUTHOR:    Luis E. Samaniego-Eguiguren, UFZ
  !    UPDATES
  !               Created        Sa   21.03.2006
  !               Last Update    Sa   11.06.2010
  !**************************************************************************
  subroutine ReadData()
    use mainVar,    only: noDataValue, cellFactor, DataConvertFactor, OffSet, &
        yStart, mStart, dStart, yEnd, mEnd, dEnd, tBuffer,DEMNcFlag
    use mo_kind,    only: i4, dp
    use mo_julian , only: NDAYS

    use kriging,    only: maxDist
    use VarFit,     only: vType, nParam, dh, hMax, beta
    use runControl, only: interMth, fnameDEM, DataPathOut, DataPathIn, fNameSta, correctNeg, &
                          distZero, flagVario, fNameVario, flagEDK
    use NetCDFVar,  only: FileOut, author_name, projection_name, variable_name, variable_unit, &
                          variable_standard_name, variable_calendar_type, &
                          variable_long_name, ncIn_variable_name, ncIn_dem_variable_name, &
                          ncIn_yCoord_name, ncIn_xCoord_name, invert_y,  &
                          ncOut_dem_variable_name, ncOut_dem_yCoord_name, ncOut_dem_xCoord_name, &
                          ncOut_dem_Latitude, ncOut_dem_Longitude
    use mo_message, only: message
    use mo_string_utils, only: divide_string

    implicit none
    !
    logical                       :: isNcFile
    logical                       :: isDEMNcFile
    integer(i4)                   :: i, ios
    character(256)                :: dummy, fileName
    character(256), allocatable   :: DataPathIn_parts(:)
    character(256), allocatable   :: fNameDEM_parts(:)
    !
    !===============================================================
    !  namelist definition
    !===============================================================
    !
    namelist/mainVars/noDataValue, DataPathIn, fNameDEM,                    &
         DataPathOut, FileOut, fNameSTA, cellFactor, DataConvertFactor, OffSet, InterMth, correctNeg,  &
         distZero, author_name, projection_name, variable_name, variable_unit, variable_long_name, &
         variable_standard_name, variable_calendar_type,  &
         yStart, mStart, dStart, yEnd, mEnd, dEnd,tBuffer, maxDist, flagVario, vType, nParam,  &
         fNameVario, dh, hMax, ncIn_variable_name, ncIn_dem_variable_name, ncIn_yCoord_name, &
         ncIn_xCoord_name, invert_y, &
         ncOut_dem_variable_name, ncOut_dem_yCoord_name, ncOut_dem_xCoord_name, ncOut_dem_Latitude,&
         ncOut_dem_Longitude
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

    ! determine whether DEM is an nc file or not
    call divide_string(fNameDEM, '.', fNameDEM_parts)
    isDEMNcFile = (fNameDEM_parts(size(fNameDEM_parts,1)) .eq. 'nc')
    deallocate(fNameDEM_parts)
    DEMNcFlag = 0
    ! read DEM
    if (isDEMNcFile) then
    ! read the netcdf DEM file
    call ReadDEMNc
    DEMNcFlag = 1
    else
    call ReadDEM
    end if

    ! determine whether input path is actually a nc file
    call divide_string(DataPathIn, '.', DataPathIn_parts)
    isNcFile = (DataPathIn_parts(size(DataPathIn_parts,1)) .eq. 'nc')
    deallocate(DataPathIn_parts)

    if (isNcFile) then
      ! read 3d meteo file
      call initMetStaNc

    else
      ! read look up table for meteorological stations
      call ReadStationLut

      ! read whole METEO data
      call ReadDataMeteo
    end if

1   format (a256)

  end subroutine ReadData


  subroutine ReadStationLut
    use mo_kind,    only: i4, dp

    use mainVar,    only: nSta, MetSta, grid, noDataValue
    use kriging,    only: maxDist, edk_dist
    use VarFit,     only: hmax
    use runControl

    implicit none
    !
    integer(i4)                         :: iSta
    integer(i4)                         :: n_stations_all  ! number of stations in file
    character(256)                      :: dummy, filename
    real(dp), dimension(:), allocatable :: temp_id         ! ID of read in station
    real(dp), dimension(:), allocatable :: temp_x          ! easting coordinate of station
    real(dp), dimension(:), allocatable :: temp_y          ! northing coordinate of station
    real(dp), dimension(:), allocatable :: temp_h          ! heigth of station
    type(extend)                        :: domain

    !-----------------------------------------------------------------------
    !	                   STATIONS COORDINATES AND ELEVATION
    !-----------------------------------------------------------------------

    call getSearchDomain(domain)

    filename = trim(fNameSTA)
    open (unit=20, file=filename, status='old', action='read')
    read (20,*)  dummy, n_stations_all
    read (20, 1) dummy

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
      if ((temp_x(nSta) .GE. domain%left)   .AND. (temp_x(nSta) .LE. domain%right) .AND. &
          (temp_y(nSta) .GE. domain%bottom) .AND. (temp_y(nSta) .LE. domain%top)  ) then
        nSta            = nSta + 1
      end if
    end do

    nSta = nSta - 1

    close(20)

    ! save relevant stations in MetSta
    allocate(MetSta(nSta))
    edk_dist%nstat = nSta
    edk_dist%noDataValue = noDataValue
    allocate(edk_dist%stat_pos(nSta, 2))

    do iSta = 1, nSta
      MetSta(iSta)%Id =  temp_id(iSta)
      MetSta(iSta)%x  =  temp_x(iSta)
      MetSta(iSta)%y  =  temp_y(iSta)
      MetSta(iSta)%h  =  temp_h(iSta)
      edk_dist%stat_pos(iSta, :) = [MetSta(iSta)%x, MetSta(iSta)%y]
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

 ! -----------------------------------------------------------------------
 !               DEM Blocks (for Irregular grids in NetCDF format )
 ! -----------------------------------------------------------------------
 subroutine ReadDEMNc
   use mainVar
   use mo_kind, only           : i4
   use runControl
   use mo_netCDF,               only:NcDataset, NcVariable
   use NetCDFVar,               only:ncOut_dem_variable_name, ncOut_dem_yCoord_name, &
                                     ncOut_dem_xCoord_name, ncOut_dem_Latitude, ncOut_dem_Longitude

   implicit none

   integer(i4)                 :: i,j
   character(256)              :: dummy, fileName
   real(dp), allocatable       :: dem_x(:,:)           ! easting
   real(dp), allocatable       :: dem_y(:,:)           ! northing
   real(dp), allocatable       :: dem_data(:,:)        ! dem data
   type(NcDataset)             :: ncin
   type(NcVariable)            :: ncvar

   fileName = trim(fNameDEM)

   ! get northing and easting information
   ncin = NcDataset(fileName, 'r')
   ncvar = ncin%getVariable(ncOut_dem_yCoord_name)
   call ncvar%getData(dem_y)
   ncvar = ncin%getVariable(ncOut_dem_xCoord_name)
   call ncvar%getData(dem_x)

   ! get number of rows and columns
   grid%nrows = size(dem_x, dim=1)
   grid%ncols = size(dem_x, dim=2)

   !write(*,*),"Nrows: ",grid%nrows
   !write(*,*),"Ncols: ",grid%ncols

   ! get latitude and longitude
   if(.not. allocated(gridMeteo%latitude)) allocate(gridMeteo%latitude(grid%ncols))
   if(.not. allocated(gridMeteo%longitude)) allocate(gridMeteo%longitude(grid%nrows))

   ncvar = ncin%getVariable(ncOut_dem_Latitude)
   call ncvar%getData(gridMeteo%latitude)
   ncvar = ncin%getVariable(ncOut_dem_Longitude)
   call ncvar%getData(gridMeteo%longitude)

   ! read in the dem data
   if (.not. allocated(G)) allocate(G(grid%nrows, grid%ncols))
   if (.not. allocated(gridMeteo%easting))  allocate(gridMeteo%easting(grid%nrows, grid%ncols))
   if (.not. allocated(gridMeteo%northing)) allocate(gridMeteo%northing(grid%nrows, grid%ncols))
   !if (.not. allocated(grid%easting))  allocate(grid%easting(grid%nrows, grid%ncols))
   !if (.not. allocated(grid%northing)) allocate(grid%northing(grid%nrows, grid%ncols))


   ncvar = ncin%getVariable(ncOut_dem_variable_name)
   call ncvar%getData(dem_data)
   call ncvar%getAttribute('_FillValue',grid%nodata_value)

   ! store the dem data in G%h
   do i=1,grid%ncols
     do j=1,grid%nrows
       G(j,i)%h = dem_data(j,i)
       gridMeteo%northing(j,i) = dem_y(j,i)
       gridMeteo%easting(j,i)  = dem_x(j,i)
       !grid%northing(i,j) = dem_y(i,j)
       !grid%easting(i,j)  = dem_x(i,j)
     end do
   end do


   !write(*,*),"DEM data: ",shape(G%h)
   !write(*,*),"Northing: ",shape(gridMeteo%northing)
   !write(*,*),"Easting: ",shape(gridMeteo%easting)
   !write(*,*),"Northing,Easting First Row, First Column: ",gridMeteo%northing(1,1),",",gridMeteo%easting(1,1)
   !write(*,*),"Northing,Easting Last Row, Last Column: ",gridMeteo%northing(192,64),",",gridMeteo%easting(192,64)

   gridMeteo%nodata_value = grid%nodata_value
   gridMeteo%nrows = grid%nrows
   gridMeteo%ncols = grid%ncols
   nCell = grid%ncols * grid%nrows
   deallocate(dem_data, dem_y, dem_x)
   ! close netcdf file
   call ncin%close()

 end subroutine ReadDEMNc

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
  ! :UPDATES         Sa.2010-07-06     reading formats / conversion units
  !
  !	**************************************************************************
  subroutine ReadDataMeteo
    use mo_kind, only         : i4, dp
    use mo_julian, only       : NDAYS
    use kriging, only         : edk_dist
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

    allocate(edk_dist%stat_z(nSta, jStart:jEnd))

    nStaEfec = 0
    do i = 1, nSta
      allocate (MetSta(i)%z(jStart:jEnd))
      MetSta(i)%z = noDataValue
      !
      ! read yearly data file
      write (dummy, 2) MetSta(i)%Id
      fileName = trim(dataPathIn)//trim(dummy)
      ! print *, 'read file: '//trim(fileName)
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
        !if ( DataConvertFactor /= 1.0_dp ) then
        where ( MetSta(i)%z(:) /= noDataValue  )   MetSta(i)%z(:) = MetSta(i)%z(:) * DataConvertFactor + OffSet
        !end if
        ! check consistency
        if (jDay /= jE ) then
          print *, 'File ', trim(fileName), ' has missing values.',  jE - jDay
          stop
        end if
      else
        close (60)
      end if
      edk_dist%stat_z(i,:) = MetSta(i)%z ! fill dist container
    end do
    print *, nStaEfec, ' of ', nSta, ' have data from ', yStart , ' to ', yEnd
    !
    ! formats
1   format (a256)
2   format (i5.5,'.dat')
3   format (i4, 1x, i2, 1x, i2)
    !
  end subroutine ReadDataMeteo

  subroutine initMetStaNc

    use mo_kind,        only: i4, dp

    use mainVar,        only: nSta, MetSta, grid, period, noDataValue, &
        yStart, mStart, dStart, yEnd, mEnd, dEnd, jStart, jEnd, DataConvertFactor, OffSet
    use kriging,        only: maxDist, edk_dist
    use VarFit,         only: hmax
    use runControl,     only: interMth, DataPathIn
    use mo_netCDF,      only: NcDataset, NcVariable
    use NetCDFVar,      only: ncIn_variable_name, ncIn_yCoord_name, ncIn_xCoord_name, ncIn_dem_variable_name
    use mo_edk_get_nc_time, only: get_time_vector_and_select

    implicit none
    !
    logical, allocatable  :: mask(:, :)
    integer(i4)           :: iSta, i, j, nrows, ncols
    integer(i4)           :: n_stations_all  ! number of stations in file
    integer(i4)           :: time_start, time_count
    character(256)        :: dummy, filename
    real(dp)              :: missing_value
    real(dp), allocatable :: temp_x(:, :)          ! easting coordinate of station
    real(dp), allocatable :: temp_y(:, :)          ! northing coordinate of station
    real(dp), allocatable :: temp_h(:, :)          ! heigth of station
    real(dp), allocatable :: met_data(:, :, :)     ! meteorological data at station
    type(period)          :: target_period
    type(extend)          :: domain
    type(NcDataset)       :: ncin
    type(NcVariable)      :: ncvar
    type(NcVariable)      :: ncdem

    call getSearchDomain(domain)

    ! setup target period
    call target_period%init(dStart, mStart, yStart, dEnd, mEnd, yEnd)
    jStart = target_period%julStart
    jEnd = target_period%julEnd
    print *, jStart, jEnd

    ! open netcdf file
    ncin = NcDataset(dataPathIn, 'r')

    ncvar = ncin%getVariable(ncIn_yCoord_name)
    call ncvar%getData(temp_y)
    ncvar = ncin%getVariable(ncIn_xCoord_name)
    call ncvar%getData(temp_x)

    !print *,"Meteo Easting Min: ",minval(temp_x)
    !print *,"Meteo Easting Max: ",maxval(temp_x)
    !print *,"Meteo Northing Min: ",minval(temp_y)
    !print *,"Meteo Northing Max: ",maxval(temp_y)

    nrows = size(temp_x, dim=1)
    ncols = size(temp_x, dim=2)

    ! get time variable
    ncvar = ncin%getVariable('time')
    ! read the time vector and get start index and count of selection
    call get_time_vector_and_select(ncvar, dataPathIn, -1, time_start, time_count, target_period)
    ! ! extract data and select time slice
    ! call ncvar%getData(data, start = (/1, 1, time_start/), cnt = (/nRows, nCols, time_count/))

    ! read dem from nc
    if (interMth .eq. 2) then
       ncdem = ncin%getVariable(ncIn_dem_variable_name)
       call ncdem%getData(temp_h, start = (/1, 1/), cnt = (/nRows, nCols/))
    end if

   !write(*,*),"DEM Source: ",shape(temp_h)

    ! read meteo data from nc
    ncvar = ncin%getVariable(ncIn_variable_name)
    call ncvar%getData(met_data, start = (/1, 1, time_start/), cnt = (/nRows, nCols, time_count/))
    !call ncvar%getAttribute('missing_value', missing_value)
    call ncvar%getAttribute('_FillValue',missing_value)

    !write(*,*),"Data Source: ",shape(met_data)

    ! close netcdf file
    call ncin%close()

    ! calculate initial mask based on meteorological data
    mask  = (met_data(:, :, 1) .ne. missing_value)

    n_stations_all = count(mask)
    nSta = 0
    do i = 1, nrows
      do j = 1, ncols
        if (.not. mask(i, j)) cycle

        ! check if station is within search distance
        ! if yes increase station counting by one
        ! otherwise overwrite station in next iteration
        !print *,"Meteo Grid Easting: ",temp_x(i,j)
        !print *,"DEM Grid Easting LL: ",domain%left
        !print *,"DEM Grid Easting TR: ", domain%right
        if ((temp_x(i, j) .GE. domain%left)   .AND. (temp_x(i, j) .LE. domain%right) .AND. &
            (temp_y(i, j) .GE. domain%bottom) .AND. (temp_y(i, j) .LE. domain%top)  ) then
          nSta = nSta + 1
        else
          mask(i, j) = .False.
        end if
      end do
    end do

    print*, ""
    print*, "Only " , nSta , ' of the given ' , n_stations_all, ' are considered since they are within the search distance.'

    ! >>> populate MetSta object
    allocate(MetSta(nSta))
    ! >>> populate the dynamic distance object
    edk_dist%nstat = nSta
    edk_dist%noDataValue = noDataValue
    allocate(edk_dist%stat_pos(nSta, 2))
    allocate(edk_dist%stat_z(nSta, target_period%julStart:target_period%julEnd))

    iSta = 0
    do i = 1, nrows
      do j = 1, ncols
        if (.not. mask(i,j)) cycle
        iSta = iSta + 1
        MetSta(iSta)%Id = iSta
        MetSta(iSta)%x  = temp_x(i, j)
        MetSta(iSta)%y  = temp_y(i, j)
        MetSta(iSta)%h  = temp_h(i, j)
        !
        ! add meteorological data
        allocate(MetSta(iSta)%z(target_period%julStart:target_period%julEnd))
        MetSta(iSta)%z = met_data(i, j, :) * DataConvertFactor + OffSet
        !
        edk_dist%stat_pos(iSta, :) = [MetSta(iSta)%x, MetSta(iSta)%y]
        edk_dist%stat_z(iSta, :) = MetSta(iSta)%z
      end do
   end do

    deallocate(mask, temp_x, temp_y, met_data) !, temp_h)

  end subroutine initMetStaNc

  subroutine getSearchDomain(domain)

    use mainVar,    only: grid,gridMeteo,DEMNcFlag
    use kriging,    only: maxDist
    use VarFit,     only: hMax
    use runControl, only: flagVario

    implicit none

    type(extend), intent(out) :: domain

    real(dp) :: search_distance

    ! take only stations which are in search distance
    ! two possibilties for search distance for vario - hMax for EDK - maxDist
    if (flagVario) then
      search_distance = hMax
    else
      search_distance = maxDist
    end if

    !print *,"DEMNcFlag: ",DEMNcFlag
    if (DEMNcFlag == 0) then
    ! derive borders for station search
      domain%bottom = -search_distance + grid%yllcorner
      domain%top    =  search_distance + grid%yllcorner + grid%nrows * grid%cellsize
      domain%left   = -search_distance + grid%xllcorner
      domain%right  =  search_distance + grid%xllcorner + grid%ncols * grid%cellsize
    else
      domain%bottom = -search_distance + minval(gridMeteo%northing)
      domain%top    =  search_distance + maxval(gridMeteo%northing)
      domain%left   = -search_distance + minval(gridMeteo%easting)
      domain%right  =  search_distance + maxval(gridMeteo%easting)
    end if

    !print *,"DEM Min Northing: ",minval(gridMeteo%northing)
    !print *,"DEM Max Northing: ",maxval(gridMeteo%northing)
    !print *,"DEM Min Easting: ",minval(gridMeteo%easting)
    !print *,"DEM Max Easting: ",maxval(gridMeteo%easting)

  end subroutine getSearchDomain


  !  *************************************************************************
  !
  !  SUBROUTINE   WRITE HEADER GRIDs
  !
  !  *************************************************************************
  subroutine writeHeader(ic, fName)
    use mainVar
    use mo_kind, only         : i4
    implicit none
    integer(i4), intent(in)        :: ic                  ! input channel
    character(256), intent(in)    :: fName
    !
    open (unit=ic, file=fName, status='unknown')
    write (ic, 1)  'ncols       ',    gridMeteo%ncols
    write (ic, 1)  'nrows       ',    gridMeteo%nrows
    write (ic, 2)  'xllcorner   ',    gridMeteo%xllcorner
    write (ic, 2)  'yllcorner   ',    gridMeteo%yllcorner
    write (ic, 1)  'cellsize    ',    gridMeteo%cellsize
    write (ic, 1)  'NODATA_value',    gridMeteo%nodata_value
    close (ic)
    ! formats
    1 format (a12, 2x, i10)
    2 format (a12, 2x, f10.1)
    !
  end subroutine writeHeader
end module mo_edk_read_data
