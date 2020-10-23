!**********************************************************************************
!  VARIOGRAM: Seting or estimating and fitting
!  PURPOSE:
!          1) Set variagram from DB for a block, or
!          2.1) Estimate an empirical semi variogram for daily met. values
!          2.2) Fitting a teoretical variogram
!          2.3) Set variogram for a block
!  WHERE:
!  UPDATES:
!          Created      Sa  19.02.2004     main structure
!          Last update      12.04.2006
!**********************************************************************************
module mo_setVario

  implicit none

  private

  public :: setVario
  public :: dMatrix
  public :: tVar

contains

subroutine setVario(param)
  use runControl
  use VarFit
  use mainVar
  use mo_EmpVar, only: EmpVar
  use mo_kind, only: i4, dp
  implicit none
  real(dp), intent(out) :: param(3)
  integer(i4) :: jd, y, y0
  !
  ! Estimation
  if (flagVario) then
     ! Initialize
     nobs=0 ; m0=0.0_dp ;  v0=0.0_dp
     y0 = 0
     !
     do jd= jStart, jEnd
        y = floor(float(jd-jStart)/365.)
        if (y > y0 ) then
           y0 = y
           print*, 'VarFit. Processing ', yStart+y0-1 , '...'
        end if
        ! update Variogram bins daily
        call EmpVar(jd, .true.)
     end do
     !
     ! variance
     v0 = v0 / real(nobs, dp)
     ! optimize parameters
     ! open(UNIT = 6,FILE = 'Report_OPT.sol', STATUS='REPLACE')

     print *, "ST: replace old GRG2 opti with something better"
     call opti(param)
  else
    print *, "   ... no variogram estimation, parameters given."
  end if
end subroutine setVario

!
!***********************************************************
!
!  tVAR:: Function to calulate variogram at given distance
!
!***********************************************************
real(8) function tVar(h,c0,c,a)
  use VarFit,  only      : vType
  use mo_kind, only      : dp
  real(dp), intent(in)  :: h                ! distance
  real(dp), intent(in)  :: c0               ! nugget = beta(1) = XU(1)
  real(dp), intent(in)  :: c                ! sill   = beta(2) = XU(2)
  real(dp), intent(in)  :: a                ! range  = beta(3) = XU(3)
  real(dp)              :: r
  !
  select case (vType)
    !
    case (1)
      ! composed:   nugget + spherical + sill
      r = h/a
      if (h == 0.0_dp) then
        tVar = c0 ! 0.0_dp
      elseif ( h <= a) then
        tVar = c0 + c * (1.5_dp * r - 0.5_dp * r**3)
      else
        tVar = c0 + c
      end if
      !
    case (2)
      ! composed:   nugget + exponential + sill
      r = h/a
      tVar = c0 + c * (1.0_dp - dexp(-r))
    end select
  !
end function tVar

!
!*******************************************************
!
!  DMATRIX:: To calculate distance between pairs.......
!
!*******************************************************
subroutine dMatrix
  use mo_kind, only                : i4, dp
  use mainVar
  use kriging
  use runControl

  implicit none
  integer(i4)                     :: i, j, k
  integer(i4)                     :: r, c, ii, jj
  integer(i4)                     :: delta, nTcell
  real(dp)                        :: xc, yc
  integer(i4), allocatable        :: list(:)
  !
  ! Initialize variables
  if ( allocated(cell)) deallocate (cell)
  !
  print*, nCell, "cells", nSta, "stations"
  edk_dist%ncell = nCell
  allocate(edk_dist%cell_pos(nCell, 2))
  allocate ( cell(nCell)     )
  allocate ( list(nSta)      )

  ! cell coordinates and elevation : checked OK
  ! ***************************************
  ! cell numbering convention (1DIM first)
  ! c->1     2                  ncol
  ! r
  ! ---+-----+------...+...-----+
  !    1     nr+1
  !    2
  !    ...             k
  !    nr    2nr                nCell
  ! ***************************************
  !    ii         column in finer grid
  !    yy         row    in finer grid
  !    (xc,yc)    coordinates of meteogrid
  ! ***************************************
  !
  r=1
  c=0
  if (DEMNcFlag /= 1) xc = gridMeteo%xllcorner + dble(gridMeteo%cellsize) * 0.5_dp
  delta = cellFactor / 2
  jj = delta
  do k=1,nCell
     ! advancing the counters
     if (r == 1) then
        c = c + 1
        if (c > 1) then
           jj = jj + cellFactor
        end if
        ii = delta
     else
        ii = ii + cellFactor
     end if

     if (DEMNcFlag == 1) then
        cell(k)%x = gridMeteo%easting(r,c)
        cell(k)%y = gridMeteo%northing(r,c)
     else
        if (r == 1) then
           if (c > 1) then
              xc = xc + dble(gridMeteo%cellsize)
           end if
           yc = gridMeteo%yllcorner + dble(gridMeteo%cellsize) * (dble(gridMeteo%nrows) - 0.5_dp)
        else
           yc = yc - dble(gridMeteo%cellsize)
        end if
        cell(k)%x = xc
        cell(k)%y = yc
      end if
      edk_dist%cell_pos(k,:) = [cell(k)%x, cell(k)%y]

    ! average of only four DEM cells around centre cell (from lower grid scale upto higher grid cell)
    !cell(k)%h = 0.25_dp*(G(ii,jj)%h + G(ii,jj+1)%h + G(ii+1,jj)%h + G(ii+1,jj+1)%h)
    !
    ! average of all DEM cells (from lower grid scale upto higher grid cell)
     if (cellFactor == 1) then
        cell(k)%h = G(ii+1,jj+1)%h
     else
        nTcell =  count(G( (ii-delta+1):(ii+delta) , (jj-delta+1):(jj+delta) )%h  > grid%nodata_value )
        if (nTcell == 0) then
           cell(k)%h = gridMeteo%nodata_value
        else
           cell(k)%h  = sum(G( (ii-delta+1):(ii+delta) , (jj-delta+1):(jj+delta) )%h, &
                G( (ii-delta+1):(ii+delta) , (jj-delta+1):(jj+delta) )%h /= gridMeteo%nodata_value ) / dble(nTcell)
        end if
     end if

     ! advance the counters
    r=r+1
    if (r > gridMeteo%nrows) r = 1
  end do

  ! find the closest stations to cell i (any order): checked  OK
  do i=1,nCell
    list = -9
    do j=1,nSta
      if (edk_dist%getCS(i,j) <= maxDist) list(j) = j
    end do
    cell(i)%nNS = count(list > -9)
    allocate ( cell(i)%listNS( cell(i)%nNS ) )
    cell(i)%listNS = pack(list, MASK = list >-9)
  end do
  !
  deallocate (list)

end subroutine dMatrix

end module mo_setVario
