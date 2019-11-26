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
  integer(i4)                                  :: yStart              ! starting year
  integer(i4)                                  :: mStart              ! starting month
  integer(i4)                                  :: dStart              ! starting day
  integer(i4)                                  :: yEnd                ! ending year
  integer(i4)                                  :: mEnd                ! ending month
  integer(i4)                                  :: dEnd                ! ending day
  integer(i4)                                  :: jStart              ! julian day start
  integer(i4)                                  :: jEnd                ! julian day end
  integer(i4)                                  :: tDays               ! number of days
  integer(i4)                                  :: nSta                ! number of stations for a block
  integer(i4)                                  :: nCell               ! number of cells to estimate z
  integer(i4)                                  :: cellFactor          ! > 1 , size grid metereological data
  real(dp)                                     :: DataConvertFactor   ! precipitation & temperature(in 1/10 mm) **** only in NECKAR BASIN *****
  real(dp)                                     :: noDataValue 
  real(dp)                                     :: thresholdDist        ! treshold cellsize  distance          
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
  end type gridGeoRef
  type (gridGeoRef)                            :: grid
  type (gridGeoRef)                            :: gridMeteo           ! reference of the metereological variables
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
  real(dp)                                     :: maxDist             ! max distance [m] search stations
  type CellCoarser
    integer(i4)                                :: nNS                 ! No. Nearest Stations (NS) d <= maxDist
    integer(i4), dimension(:), allocatable     :: listNS              ! list of NS
    real(dp)                                   :: x                   ! x- coordinate
    real(dp)                                   :: y                   ! y- coordinate
    real(dp)                                   :: h                   ! (estimated) elevation [m] (from the nearest cells DEM)
    real(sp)                                   :: z                   ! z values to be interpolated (OUTPUT)
  end type CellCoarser
  type(CellCoarser),  dimension(:), allocatable  :: cell             ! EDK output

  real(dp), dimension(:,:), allocatable        :: dCS                 ! Euclidean distance between cells -> stations
  type dtoS
    real(dp), dimension(:), allocatable        :: S                   ! distance to Station j
  end type dtoS
  type (dtoS), dimension(:), allocatable       :: dS                  ! distance from Station i to all js
  type (dtoS), dimension(:), allocatable       :: dz2S                ! squared diference of Z values
  real(dp)                                     :: xl, xr, yd, yu      ! coordinates of the interpolation block
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
end module NetCDFVar
