!> \file    mo_edk_modules.f90
!> \brief   Modules with global variables for EDK
!> \details This file contains the MAIN SHARED variables
!!          - `mainvar`: main setup variables
!!          - `runcontrol`: run configuration
!!          - `kriging`: kriging setup variables
!!          - `varfit`: variogram fitting variables
!!          - `netcdfvar`: NetCDF I/O definition variables
!!
!> \author  Luis Samaniego
!> \date    22.03.2006
!> \date    24.03.2006
!> \author  Sebastian Mueller
!> \date    June 2022
!!          - refactored

!> \brief   main variables
!> \copyright Copyright 2005-\today, the CHS Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! EDK is released under the LGPLv3+ license \license_note
module mainVar
  use mo_kind, only: i4, sp, dp
  implicit none
  ! parameters
  integer(i4)                                  :: yStart               !< starting year
  integer(i4)                                  :: mStart               !< starting month
  integer(i4)                                  :: dStart               !< starting day
  integer(i4)                                  :: yEnd                 !< ending year
  integer(i4)                                  :: mEnd                 !< ending month
  integer(i4)                                  :: dEnd                 !< ending day
  integer(i4)                                  :: jStart               !< julian day start
  integer(i4)                                  :: jEnd                 !< julian day end
  integer(i4)                                  :: tBuffer              !< number of days for time buffering
  integer(i4)                                  :: tDays                !< number of days
  integer(i4)                                  :: nSta                 !< number of stations for a block
  integer(i4)                                  :: nCell                !< number of cells to estimate z
  integer(i4)                                  :: cellFactor           !< > 1 , size grid metereological data
  integer(i4)                                  :: DEMNcFlag            !< flag for DEM format 0 = text file, 1 = netCDF
  real(dp)                                     :: DataConvertFactor    !< precipitation & temperature(in 1/10 mm) **** only in NECKAR BASIN *****
  real(dp)                                     :: OffSet               !< constant to be added (Ex: add  273 to convert tavg from C to K )
  real(dp)                                     :: noDataValue          !< no data value of station input data
  real(dp)                                     :: thresholdDist        !< treshold cellsize  distance
  ! constants
  real(dp),  parameter                         :: DayHours = 24.0_dp   !< hours per day
  real(dp),  parameter                         :: YearDays = 365.0_dp  !< days in a year
  real(dp),  parameter                         :: DaySecs = 86400.0_dp !< sec in a day
  real(sp), public, parameter                  :: nodata_sp = -9999._sp !< [-]     global no data value (single precision)

  !> \class cellfiner
  !> \brief cell elevation
  type CellFiner
     real(dp)                                  :: h                    !< elevation (sinks removed) [m]
  end type CellFiner
  type(CellFiner), dimension(:,:), allocatable :: G                    !< Cell characteristics

  !> \class meteostation
  !> \brief Meteo Stations type
  type MeteoStation
    integer(i4)                                :: Id                   !< Id number
    real(dp)                                   :: x                    !< x coordinate
    real(dp)                                   :: y                    !< y coordinate
    real(dp)                                   :: h                    !< elevation
    real(dp), dimension(:), allocatable        :: z                    !< observed daily value (prec. temp, etc)
  end type MeteoStation
  type(MeteoStation),   dimension(:), allocatable :: MetSta            !< Meteo Stations

  !> \class gridgeoref
  !> \brief GRID description
  type gridGeoRef
    integer(i4)                                :: ncols                !< number of columns
    integer(i4)                                :: nrows                !< number of rows
    real(dp)                                   :: xllcorner            !< x coordinate of the lowerleft corner
    real(dp)                                   :: yllcorner            !< y coordinate of the lowerleft corner
    integer(i4)                                :: cellsize             !< cellsize x = cellsize y
    integer(i4)                                :: nodata_value         !< code to define the mask
    real(dp), dimension(:,:), allocatable      :: easting              !< irregular grid easting
    real(dp), dimension(:,:), allocatable      :: northing             !< irregular grid northing
    real(dp), dimension(:), allocatable        :: latitude             !< latitude for the output
    real(dp), dimension(:), allocatable        :: longitude            !< longitude for the output
  end type gridGeoRef
  type (gridGeoRef)                            :: grid                 !< grid
  type (gridGeoRef)                            :: gridMeteo            !< reference of the metereological variables

  !> \class period
  !> \brief PERIOD description
  type period
    integer(i4) :: dStart      !< first day
    integer(i4) :: mStart      !< first month
    integer(i4) :: yStart      !< first year
    integer(i4) :: dEnd        !< last  day
    integer(i4) :: mEnd        !< last  month
    integer(i4) :: yEnd        !< last  year
    integer(i4) :: julStart    !< first julian day
    integer(i4) :: julEnd      !< last  julian day
    integer(i4) :: nObs        !< total number of observations
  CONTAINS
    !> \copydoc mainvar::init
    procedure :: init !< \see mainvar::init
  end type period

contains

  !> \brief   initialize a period
  subroutine init(self, dStart, mStart, yStart, dEnd, mEnd, yEnd)

    use mo_julian, only: julday

    implicit none

    class(period), intent(inout) :: self
    integer(i4), intent(in) :: dStart, mStart, yStart, dEnd, mEnd, yEnd

    self%dStart = dStart
    self%mStart = mStart
    self%yStart = yStart
    self%dEnd = dEnd
    self%mEnd = mEnd
    self%yEnd = yEnd

    self%julStart = julday(dd = dStart, mm = mStart, yy = yStart)
    self%julEnd   = julday(dd = dEnd, mm = mEnd, yy = yEnd)
    self%nObs     = self%julEnd - self%julStart + 1_i4

  end subroutine init

end module mainVar

!> \brief   RUN Control
!> \copyright Copyright 2005-\today, the CHS Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! EDK is released under the LGPLv3+ license \license_note
module runControl
  use mo_kind, only: i4, sp
  ! timer
  character(256)                               :: DataPathIn
  character(256)                               :: DataPathOut
  character(256)                               :: DataPathDEM
  character(256)                               :: fNameDEM            !< DEM file name
  character(256)                               :: fNameSTA            !< Station id and coordinates
  character(256)                               :: fNameVARIO          !< file name of variogram parameters
  character(256)                               :: prefix              !< prefix data file
  integer(i4)                                  :: nBlocks             !< number of interpolation blocks
  logical                                      :: flagEDK             !< estimate EDK (T/F)
  logical                                      :: flagVario           !< estimate variogran (T/F)
  logical                                      :: correctNeg          !< correct negative interpolated values
  logical                                      :: distZero            !< values further away than dist threshold are set to zero
  integer(i4)                                  :: interMth            !< interp. method
end module runControl

!> \brief   kriging variables
!> \copyright Copyright 2005-\today, the CHS Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! EDK is released under the LGPLv3+ license \license_note
module kriging
  use mo_kind, only: i4, sp, dp
  use mo_edk_types, only: dist_t
  real(dp)                                     :: maxDist             !< max distance [m] search stations

  !> \class cellcoarser
  !> \brief cell coarser type
  type CellCoarser
    integer(i4)                                :: nNS                 !< No. Nearest Stations (NS) d <= maxDist
    integer(i4), dimension(:), allocatable     :: listNS              !< list of NS
    real(dp)                                   :: x                   !< x- coordinate
    real(dp)                                   :: y                   !< y- coordinate
    real(dp)                                   :: h                   !< (estimated) elevation [m] (from the nearest cells DEM)
    real(sp),                  allocatable     :: z(:)                !< z values to be interpolated (OUTPUT)
  end type CellCoarser
  type(CellCoarser),  dimension(:), allocatable  :: cell              !< EDK output

  type(dist_t)                                 :: edk_dist            !< distance calculations for EDK
  real(dp)                                     :: xl, xr, yd, yu      !< coordinates of the interpolation block
end module kriging

!> \brief   variogram fitting variables
!> \copyright Copyright 2005-\today, the CHS Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! EDK is released under the LGPLv3+ license \license_note
module VarFit
  use mo_kind, only: i4, dp
  implicit none
  ! parameters
  integer(i4)                                  :: vType               !< variogram type
  integer(i4)                                  :: nParam              !< number of parameters
  integer(i4), dimension(:), allocatable       :: Nh                  !< number of pairs per bin
  integer(i4)                                  :: nbins
  integer(i4)                                  :: nobs                !< number of observations
  real(dp)                                     :: dh                  !< bin size
  real(dp)                                     :: hMax                !< max distance h variogram
  real(dp),  dimension(:,:), allocatable       :: gamma               !< variograms
  real(dp)                                     :: gmax(2)
  real(dp)                                     :: m0
  real(dp)                                     :: v0
  real(dp), dimension(:), allocatable          :: beta                !< parameters of the variogram
  real(dp), dimension(8)                       :: E                   !< efficiency measures
  real(dp), parameter                          :: sRadius = 1.d4      !< searching distance limits
  real(dp), parameter                          :: gradE   = 2.d-1     !< gradient limit = delta E / distance (searching)

end module VarFit

!> \brief   NetCDF IO specifications variables
!> \copyright Copyright 2005-\today, the CHS Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! EDK is released under the LGPLv3+ license \license_note
module NetCDFVar
  use mo_kind, only : dp
  implicit none
  ! Variables for NetCDF writing
  character(256)                                   :: fileOut    !< File Name out
  real(dp), dimension(:),   allocatable, target    :: yCoor               !< GK4 (DHDN3-zone 4) easting
  real(dp), dimension(:),   allocatable, target    :: xCoor               !< GK4 (DHDN3-zone 4) northing
  real(dp), dimension(:,:), allocatable, target    :: lons                !< WGS84 lons
  real(dp), dimension(:,:), allocatable, target    :: lats                !< WGS84 lats
  !
  character(256)                                   :: variable_name !< name of netcdf variable
  character(256)                                   :: variable_unit !< unit of netcdf variable
  character(256)                                   :: variable_long_name !< long name  of netcdf variable
  character(256)                                   :: originator !< originator of netcdf file
  character(256)                                   :: crs !< coordinate reference system (EPSG:XXXX)
  character(256)                                   :: title !< title of netcdf file
  character(256)                                   :: source !< source and methods of original data
  character(256)                                   :: contact !< contact (email address)
  character(256)                                   :: institution !< institution
  character(256)                                   :: variable_standard_name !< standard name of netcdf variable
  character(256)                                   :: variable_calendar_type !< calendar type (time variable attribute)
  logical                                          :: invert_y
  !
  ! netcdf input specifications
  character(256)                                   :: ncIn_variable_name !< netcdf input
  character(256)                                   :: ncIn_dem_variable_name !< netcdf input
  character(256)                                   :: ncIn_yCoord_name !< netcdf input
  character(256)                                   :: ncIn_xCoord_name !< netcdf input
  character(256)                                   :: ncOut_dem_variable_name !< netcdf input
  character(256)                                   :: ncOut_dem_yCoord_name !< netcdf input
  character(256)                                   :: ncOut_dem_xCoord_name !< netcdf input
  character(256)                                   :: ncOut_dem_Latitude !< netcdf input
  character(256)                                   :: ncOut_dem_Longitude !< netcdf input

end module NetCDFVar
