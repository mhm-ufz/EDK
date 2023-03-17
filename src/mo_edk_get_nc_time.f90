!> \file    mo_edk_get_nc_time.f90
!> \copydoc mo_edk_get_nc_time

!> \brief   Module to get time vector from NetCDF file.
!> \author  Matthias Zink
!> \date    Oct 2012
!> \authors Matthias Cuntz, Juliane Mai
!> \date    Nov 2014
!!          - time int or double
!> \author  Stephan Thober
!> \date    Sep 2015
!!          - added read for hourly data
!> \author  Robert Schweppe
!> \date    Nov 2017
!!          - restructured routine, reads vector now
!> \author  Maren Kaluza
!> \date    May 2018
!!          - fixed bug in time reading
!> \copyright Copyright 2005-\today, the CHS Developers, Sabine Attinger: All rights reserved.
!! EDK is released under the LGPLv3+ license \license_note
module mo_edk_get_nc_time

  implicit none

  private

  public :: get_time_vector_and_select

contains

  !> \brief   Extract time vector
  !> \details Extract time vector in unit julian hours and get supposed time step in hours
  subroutine get_time_vector_and_select(var, fname, inctimestep, time_start, time_cnt, target_period)

    use mainVar, only : period, dayhours, daysecs, yeardays
    use mo_julian, only : caldat, dec2date, julday
    use mo_kind, only : dp, i4, i8
    use mo_message, only : message
    use mo_netcdf, only : NcVariable
    use mo_string_utils, only : DIVIDE_STRING

    implicit none

    !> variable of interest
    type(NcVariable), intent(in) :: var

    !> fname of ncfile for error message
    character(256), intent(in) :: fname

    !> flag for requested time step
    integer(i4), intent(in) :: inctimestep

    !> time_start index of time selection
    integer(i4), intent(out) :: time_start

    !> time_count of indexes of time selection
    integer(i4), intent(out) :: time_cnt

    !> reference period
    type(period), intent(in), optional :: target_period

    ! reference time of NetCDF
    integer(i4) :: yRef, dRef, mRef, hRef, jRef

    ! netcdf attribute values
    character(256) :: AttValues

    ! dummies for netcdf attribute handling
    character(256), dimension(:), allocatable :: strArr, date, time

    ! native time step converter in ncfile
    integer(i8) :: time_step_seconds

    ! time vector
    integer(i8), allocatable, dimension(:) :: time_data

    ! period of ncfile, for clipping
    type(period) :: nc_period, clip_period

    integer(i4) :: ncJulSta1, dd, n_time

    integer(i4) :: mmcalstart, mmcalend, yycalstart, yycalend

    integer(i4) :: mmncstart, yyncstart

    ! helper variable for error output
    integer(i4) :: hstart_int, hend_int

    ! helper variable for error output
    character(256) :: error_msg


    call var%getAttribute('units', AttValues)
    ! AttValues looks like "<unit> since YYYY-MM-DD[ HH:MM:SS]"
    ! split at space
    call DIVIDE_STRING(trim(AttValues), ' ', strArr)

    ! determine reference time at '-' and convert to integer
    call DIVIDE_STRING(trim(strArr(3)), '-', date)
    read(date(1), *) yRef
    read(date(2), *) mRef
    read(date(3), *) dRef

    jRef = julday(dd = dRef, mm = mRef, yy = yRef)

    ! if existing also read in the time (only hour so far)
    hRef = 0
    if(size(strArr) .gt. 3) then
      call DIVIDE_STRING(trim(strArr(4)), ':', time)
      read(time(1), *) hRef
    end if

    ! determine the step_size
    if (strArr(1) .EQ. 'days') then
      time_step_seconds = int(DaySecs)
    else if (strArr(1) .eq. 'hours') then
      time_step_seconds = int(DaySecs / DayHours)
    else if (strArr(1) .eq. 'minutes') then
      time_step_seconds = int(DaySecs / DayHours / 60._dp)
    else if (strArr(1) .eq. 'seconds') then
      time_step_seconds = 1_i8
    else
      call message('***ERROR: Please provide the input data in (days, hours, minutes, seconds) ', &
              'since YYYY-MM-DD[ HH:MM:SS] in the netcdf file. Found: ', trim(AttValues))
      stop
    end if

    ! get the time vector
    call var%getData(time_data)
    ! convert array from units since to seconds
    time_data = time_data * time_step_seconds

    ! check for length of time vector, needs to be at least of length 2, otherwise step width check fails
    if (size(time_data) .le. 1) then
      call message('***ERROR: length of time dimension needs to be at least 2 in file: ' // trim(fname))
      stop
    end if

    ! check for equal timesteps and timestep must not be multiple of native timestep
    error_msg = '***ERROR: time_steps are not equal over all times in file and/or do not conform to' // &
            ' requested timestep in file (' // trim(fname) // ') : '

    ! compare the read period from ncfile to the period required
    ! convert julian second information back to date via conversion to float
    ! the 0.5_dp is for the different reference of fractional julian days, hours are truncated
    n_time = size(time_data)
    call dec2date(time_data(1) / DaySecs - 0.5_dp + jRef + hRef / 24._dp, nc_period%dStart, nc_period%mStart, &
            nc_period%yStart, hstart_int)
    nc_period%julStart = int(time_data(1) / DaySecs + jRef + hRef / 24._dp)
    call dec2date(time_data(n_time) / DaySecs - 0.5_dp + jRef + hRef / 24._dp, nc_period%dEnd, nc_period%mEnd, &
            nc_period%yEnd, hend_int)
    nc_period%julEnd = int(time_data(n_time) / DaySecs + jRef + hRef / 24._dp)

    ! if no target period is present, use the whole time period
    if (present(target_period)) then
      clip_period = target_period
    else
      clip_period = nc_period
    end if

    ! prepare the selection and check for required time_step
    select case(inctimestep)
    case(-1) ! daily
      ! difference must be 1 day
      if (.not. all(abs((time_data(2 : n_time) - time_data(1 : n_time - 1)) / DaySecs - 1._dp) .lt. 1.e-6)) then
        call message(error_msg // trim('daily'))
        stop
      end if
      ncJulSta1 = nc_period%julStart
      time_start = clip_period%julStart - ncJulSta1 + 1_i4
      time_cnt = clip_period%julEnd - clip_period%julStart + 1_i4
    case(-2) ! monthly
      ! difference must be between 28 and 31 days
      if (any(abs((time_data(2 : n_time) - time_data(1 : n_time - 1)) / DaySecs) .gt. 31._dp) .or. &
              any(abs((time_data(2 : n_time) - time_data(1 : n_time - 1)) / DaySecs) .lt. 28._dp)) then
        call message(error_msg // trim('monthly'))
        stop
      end if

      call caldat(clip_period%julStart, dd, mmcalstart, yycalstart)
      call caldat(nc_period%julStart, dd, mmncstart, yyncstart)
      ! monthly timesteps are usually set by month end, so for beginning, we need 1st of month
      ncJulSta1 = julday(1, mmncstart, yyncstart)
      call caldat(clip_period%julEnd, dd, mmcalend, yycalend)
      time_start = (yycalstart * 12 + mmcalstart) - (yyncstart * 12 + mmncstart) + 1_i4
      time_cnt = (yycalend * 12 + mmcalend) - (yycalstart * 12 + mmcalstart) + 1_i4
    case(-3) ! yearly
      ! difference must be between 365 and 366 days
      if (any(abs((time_data(2 : n_time) - time_data(1 : n_time - 1)) / DaySecs) .gt. (YearDays + 1._dp)) .or. &
              any(abs((time_data(2 : n_time) - time_data(1 : n_time - 1)) / DaySecs) .lt. YearDays)) then
        call message(error_msg // 'yearly')
        stop
      end if
      call caldat(clip_period%julStart, dd, mmcalstart, yycalstart)
      call caldat(nc_period%julStart, dd, mmncstart, yyncstart)
      ! yearly timesteps are usually set by year end, so for beginning, we need 1st of year
      ncJulSta1 = julday(1, 1, yyncstart)
      call caldat(clip_period%julEnd, dd, mmcalend, yycalend)
      time_start = yycalstart - yyncstart + 1_i4
      time_cnt = yycalend - yycalstart + 1_i4
    case(-4) ! hourly
      ! difference must be 1 hour
      if (.not. all(abs((time_data(2 : n_time) - time_data(1 : n_time - 1)) / 3600._dp - 1._dp) .lt. 1.e-6)) then
        call message(error_msg // 'hourly')
        stop
      end if
      ncJulSta1 = nc_period%julStart
      time_start = (clip_period%julStart - ncJulSta1) * 24_i4 + 1_i4 ! convert to hours; always starts at one
      time_cnt = (clip_period%julEnd - clip_period%julStart + 1_i4) * 24_i4 ! convert to hours
    case default ! no output at all
      call message('***ERROR: read_forcing_nc: unknown nctimestep switch.')
      stop
    end select

    ! Check if time steps in file cover simulation period
    if (.not. ((ncJulSta1 .LE. clip_period%julStart) .AND. (nc_period%julEnd .GE. clip_period%julEnd))) then
      call message('***ERROR: read_forcing_nc: time period of input data: ', trim(fname), &
              '          is not matching modelling period.')
      stop
    end if

  end subroutine get_time_vector_and_select

end module mo_edk_get_nc_time
