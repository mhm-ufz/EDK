!******************************************************************************                                                              
!  MODULES for EDK
!  AUTHOR  Luis Samaniego (UFZ)
!  DESCRIPTION
!          This file contains the MAIN SHARED variables
!  UPDATES
!          Created        Sa   22.03.2006
!          Last Update    Sa   24.03.2006
!******************************************************************************
module mainVar
  use mo_kind, only: i4, dp
  implicit none
  ! parameters
  character(10), parameter                     :: version = '2.0'
  integer(i4)                                  :: yStart              ! starting year
  integer(i4)                                  :: mStart              ! starting month
  integer(i4)                                  :: dStart              ! starting day
  integer(i4)                                  :: yEnd                ! ending year
  integer(i4)                                  :: mEnd                ! ending month
  integer(i4)                                  :: dEnd                ! ending day
  integer(i4)                                  :: jStart              ! julian day start
  integer(i4)                                  :: jEnd                ! julian day end
  integer(i4)                                  :: tBuffer             ! number of days for time buffering
  integer(i4)                                  :: tDays               ! number of days
  integer(i4)                                  :: nSta                ! number of stations for a block
  integer(i4)                                  :: nCell               ! number of cells to estimate z
  integer(i4)                                  :: cellFactor          ! > 1 , size grid metereological data
  integer(i4)                                  :: DEMNcFlag           ! flag for DEM format 0 = text file, 1 = netCDF  
  real(dp)                                     :: DataConvertFactor   ! precipitation & temperature(in 1/10 mm) **** only in NECKAR BASIN *****
  real(dp)                                     :: noDataValue 
  real(dp)                                     :: thresholdDist        ! treshold cellsize  distance 
  ! constants
  real(dp),  parameter                         :: DayHours = 24.0_dp   ! hours per day
  real(dp),  parameter                         :: YearDays = 365.0_dp  ! days in a year
  real(dp),  parameter                         :: DaySecs = 86400.0_dp ! sec in a day
  ! data input
  type CellFiner
     real(dp)                                  :: h                   ! elevation (sinks removed) [m]  
  end type CellFiner
  type(CellFiner), dimension(:,:), allocatable :: G                   ! Cell characteristics
  type MeteoStation
    integer(i4)                                :: Id                  ! Id number
    real(dp)                                   :: x                   ! x coordinate
    real(dp)                                   :: y                   ! y coordinate
    real(dp)                                   :: h                   ! elevation
    real(dp), dimension(:), allocatable        :: z                   ! observed daily value (prec. temp, etc)
  end type MeteoStation 
  type(MeteoStation),   dimension(:), allocatable :: MetSta 

  ! GRID description
  type gridGeoRef
    integer(i4)                                :: ncols               ! number of columns
    integer(i4)                                :: nrows               ! number of rows
    real(dp)                                   :: xllcorner           ! x coordinate of the lowerleft corner
    real(dp)                                   :: yllcorner           ! y coordinate of the lowerleft corner
    integer(i4)                                :: cellsize            ! cellsize x = cellsize y
    integer(i4)                                :: nodata_value        ! code to define the mask
    real(dp), dimension(:,:), allocatable      :: easting             ! irregular grid easting
    real(dp), dimension(:,:), allocatable      :: northing            ! irregular grid northing
    real(dp), dimension(:), allocatable        :: latitude            ! latitude for the output 
    real(dp), dimension(:), allocatable        :: longitude           ! longitude for the output
  end type gridGeoRef
  type (gridGeoRef)                            :: grid
  type (gridGeoRef)                            :: gridMeteo           ! reference of the metereological variables

  ! -------------------------------------------------------------------
  ! PERIOD description
  ! -------------------------------------------------------------------
  type period
    integer(i4) :: dStart      ! first day
    integer(i4) :: mStart      ! first month
    integer(i4) :: yStart      ! first year
    integer(i4) :: dEnd        ! last  day
    integer(i4) :: mEnd        ! last  month
    integer(i4) :: yEnd        ! last  year
    integer(i4) :: julStart    ! first julian day
    integer(i4) :: julEnd      ! last  julian day
    integer(i4) :: nObs        ! total number of observations
  CONTAINS
    procedure :: init
  end type period

contains

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

  !
  !********************************************
  ! RUN Control
  !********************************************
module runControl
  use mo_kind, only: i4, sp
  ! timer
  character(256)                               :: DataPathIn          !
  character(256)                               :: DataPathOut
  character(256)                               :: DataPathDEM
  character(256)                               :: fNameDEM            ! DEM file name
  character(256)                               :: fNameSTA            ! Station id and coordinates
  character(256)                               :: fNameVARIO          ! file name of variogram parameters
  character(256)                               :: prefix              ! prefix data file
  integer(i4)                                  :: nBlocks             ! number of interpolation blocks
  logical                                      :: flagEDK             ! estimate EDK (T/F)
  logical                                      :: flagVario           ! estimate variogran (T/F)
  logical                                      :: correctNeg          ! correct negative interpolated values
  logical                                      :: distZero            ! values further away than dist threshold are set to zero
  integer(i4)                                  :: interMth            ! interp. method
end module runControl

module kriging
  use mo_kind, only: i4, sp, dp
  use mainVar, only: nSta
  real(dp)                                     :: maxDist             ! max distance [m] search stations
  type CellCoarser
    integer(i4)                                :: nNS                 ! No. Nearest Stations (NS) d <= maxDist
    integer(i4), allocatable                                :: Nk_old(:)              ! old stations (added)
    integer(i4), dimension(:), allocatable     :: listNS              ! list of NS
    real(dp)                                   :: x                   ! x- coordinate
    real(dp)                                   :: y                   ! y- coordinate
    real(dp)                                   :: h                   ! (estimated) elevation [m] (from the nearest cells DEM)
    real(sp),                  allocatable     :: z(:)                ! z values to be interpolated (OUTPUT)
    real(dp),                  allocatable     :: W(:)
  end type CellCoarser
  type(CellCoarser),  dimension(:), allocatable  :: cell             ! EDK output

  real(dp), dimension(:,:), allocatable        :: dCS                 ! Euclidean distance between cells -> stations
  type dtoS
    real(dp), dimension(:), allocatable        :: S                   ! distance to Station j
  end type dtoS
  type (dtoS), dimension(:), allocatable       :: dS                  ! distance from Station i to all js
  type (dtoS), dimension(:), allocatable       :: dz2S                ! squared diference of Z values
  real(dp)                                     :: xl, xr, yd, yu      ! coordinates of the interpolation block              
  !real(dp)                 , allocatable       :: X(:)  
!
end module kriging

module VarFit
  use mo_kind, only: i4, dp
  implicit none
  ! parameters                                                         
  integer(i4)                                  :: vType               ! variogram type
  integer(i4)                                  :: nParam              ! number of parameters
  integer(i4), dimension(:), allocatable       :: Nh                  ! number of pairs per bin
  integer(i4)                                  :: nbins
  integer(i4)                                  :: nobs                ! number of observations
  real(dp)                                     :: dh                  ! bin size
  real(dp)                                     :: hMax                ! max distance h variogram
  real(dp),  dimension(:,:), allocatable       :: gamma
                                                                     !   variograms
  real(dp)                                     :: gmax(2)
  real(dp)                                     :: m0
  real(dp)                                     :: v0
  real(dp), dimension(:), allocatable          :: beta                ! parameters of the variogram 
  real(dp), dimension(8)                       :: E                   ! efficiency measures
  real(dp), parameter                          :: sRadius = 1.d4      ! searching distance limits
  real(dp), parameter                          :: gradE   = 2.d-1     ! gradient limit = delta E / distance (searching)

end module VarFit

module NetCDFVar
  use mo_kind, only : dp
  implicit none
  ! Variables for NetCDF writing
  character(256)                                   :: fileOut    ! File Name out
  real(dp), dimension(:),   allocatable, target    :: yCoor               ! GK4 (DHDN3-zone 4) easting
  real(dp), dimension(:),   allocatable, target    :: xCoor               ! GK4 (DHDN3-zone 4) northing
  real(dp), dimension(:,:), allocatable, target    :: lons                ! WGS84 lons 
  real(dp), dimension(:,:), allocatable, target    :: lats                ! WGS84 lats
  !
  character(256)                                   :: variable_name ! name of netcdf variable
  character(256)                                   :: variable_unit ! unit of netcdf variable
  character(256)                                   :: variable_long_name ! long name  of netcdf variable
  character(256)                                   :: author_name ! author name of netcdf file
  character(256)                                   :: projection_name ! name of EPSG (EPSG:XXXX)
  character(256)                                   :: variable_standard_name ! standard name of netcdf variable
  character(256)                                   :: variable_calendar_type ! calendar type (time variable attribute)
  logical                                          :: invert_y
  !
  ! netcdf input specifications
  character(256)                                   :: ncIn_variable_name
  character(256)                                   :: ncIn_dem_variable_name
  character(256)                                   :: ncIn_yCoord_name
  character(256)                                   :: ncIn_xCoord_name
  character(256)                                   :: ncOut_dem_variable_name
  character(256)                                   :: ncOut_dem_yCoord_name
  character(256)                                   :: ncOut_dem_xCoord_name
  character(256)                                   :: ncOut_dem_Latitude
  character(256)                                   :: ncOut_dem_Longitude

end module NetCDFVar


