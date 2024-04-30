!> \file    mo_edk_setvario.f90
!> \brief   \copybrief mo_edk_setvario
!> \details \copydetails mo_edk_setvario

!> \brief   VARIOGRAM: Seting or estimating and fitting
!> \details PURPOSE:
!!          1) Set variagram from DB for a block, or
!!          2.1) Estimate an empirical semi variogram for daily met. values
!!          2.2) Fitting a teoretical variogram
!!          2.3) Set variogram for a block
!> \author  Luis Samaniego
!> \date    19.02.2004
!!          - main structure
!> \date    12.04.2006
!> \copyright Copyright 2005-\today, the CHS Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! EDK is released under the LGPLv3+ license \license_note
module mo_edk_setvario
  use mo_kind, only: i4, dp
  implicit none

  private

  public :: setVario
  public :: dMatrix
  public :: tVar

contains

  !> \brief   find optimal variogram parameters, if wanted
  subroutine setVario(param)
    use runControl
    use VarFit
    use mainVar
    use mo_edk_empvar, only: EmpVar
    implicit none
    real(dp), intent(out) :: param(3) !< parameters (nugget, sill, range)
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

  !> \brief   Function to calulate variogram at given distance
  !> \return  variogram value
  real(dp) function tVar(h,c0,c,a)
    use VarFit,  only      : vType
    use mo_utils, only: eq
    real(dp), intent(in)  :: h                !< distance
    real(dp), intent(in)  :: c0               !< nugget = beta(1) = XU(1)
    real(dp), intent(in)  :: c                !< sill   = beta(2) = XU(2)
    real(dp), intent(in)  :: a                !< range  = beta(3) = XU(3)
    real(dp)              :: r
    !
    if ( a > 0.0_dp ) then
      r = h/a
    else
      r = 0.0_dp
    end if
    select case (vType)
      !
      case (1)
        ! composed:   nugget + spherical + sill
        if ( eq(h, 0.0_dp) ) then
          tVar = c0 ! 0.0_dp
        elseif ( h < a ) then
          tVar = c0 + c * (1.5_dp * r - 0.5_dp * r**3)
        else
          tVar = c0 + c
        end if
        !
      case (2)
        ! composed:   nugget + exponential + sill
        if ( eq(h, 0.0_dp) ) then
          tVar = c0 ! 0.0_dp
        elseif ( a > 0 ) then
          tVar = c0 + c * (1.0_dp - exp(-r))
        else
          tVar = c0 + c
        end if
      end select
    !
  end function tVar

  !> \brief   DMATRIX:: To calculate distance between pairs.
  subroutine dMatrix
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

    do i=1,nSta-1
      ! distance matrix between stations:            checked OK
      do j=i+1, nSta
        if (edk_dist%getSS(i,j) == 0.0_dp) then
          ! check if stations are closer than 5 meter
          print* , '--------------------------------------------------------------------------'
          print* , '!!! Warning: '
          print* , '!!! Stations: ', MetSta(i)%Id, ' and ', MetSta(j)%Id
          print* , '!!! have the same coordinates or are repeated. Check LUT. '
          print* , '!!! Rounded artefacts can be generated when stations have same coordinates'
          print* , '!!! and data at the time.'
          print* , '--------------------------------------------------------------------------'
          !stop
        end if
        if (edk_dist%getSS(i,j) > 0.0_dp .and. edk_dist%getSS(i,j) < 5.0_dp) then
          print* , '--------------------------------------------------------------------------'
          print* , '!!! Warning: '
          print* , '!!! Stations: ', MetSta(i)%Id, ' and ', MetSta(j)%Id
          print* , '!!!  are closer than 5 meter distance. '
          print* , '--------------------------------------------------------------------------'
        end if
      end do
    end do

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

    print*, "   ... find neighborhoods"

    ! find the closest stations to cell i (any order): checked  OK
    do i=1,nCell
      if ( modulo(i,100000) == 0 ) print*, "      ... cells ready", i, "of", nCell
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

  !> \brief   Optimization routine for variograms.
  !> \details Initialization of the Nonlinear Optimization Subroutine GRG2
  !!          The function to be optimized is suplied in subroutine GCOMP
  !> \author  Luis Samaniego
  !> \date    28.05.1999
  subroutine OPTI(pmin)
    use VarFit
    use mo_nelmin, only: nelminrange

    ! parameters for Nelder-Mead algorithm
    real(dp), intent(out) :: pmin(3) !< optimized parameter set
    real(dp) :: pstart(3) ! Starting point for the iteration.
    real(dp) :: prange(3, 2) ! Range of parameters (upper and lower bound).
    real(dp) :: varmin ! the terminating limit for the variance of the function values. varmin>0 is required
    real(dp) :: step(3) ! determines the size and shape of the initial simplex. The relative magnitudes of its elements should reflect the units of the variables. size(step)=size(start)
    integer(i4) :: konvge ! the convergence check is carried out every konvge iterations
    integer(i4) :: maxeval ! the maximum number of function evaluations. default: 1000
    real(dp) :: funcmin
    integer(i4) :: neval ! the number of function evaluations used.
    integer(i4) :: numrestart ! the number of restarts.
    integer(i4) :: ierror ! error indicator.
                          ! 0: no errors detected.
                          ! 1: varmin or konvge have an illegal value.
                          ! 2: iteration terminated because maxeval was exceeded without convergence.
    real(dp), allocatable :: history(:)


    ! ! inputs for GRG2
    ! IMPLICIT  DOUBLE PRECISION(A-H,O-Z), INTEGER(I,J,L,M,N)
    ! INTEGER*4 NCORE,NNVARS,NFUN,MAXBAS,MAXHES,LASTZ
    ! INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
    ! LOGICAL   MAXIM,INPRNT,OTPRNT
    ! DIMENSION BLVAR(100),BUVAR(100),BLCON(10),BUCON(10)
    ! DIMENSION RAMCON(10),RAMVAR(10),INBIND(10),Z(5000)
    ! DIMENSION NONBAS(10),REDGR(10),DEFAUL(20),TITLE(19)
    ! DIMENSION XX(100),FCNS(10),RMULTS(10)
    ! DATA      BIG/1.D31/

    ! Initialization of Nelder-Mead
    pstart = (/0.0, 1., 0.5/) ! Starting point for the iteration.
    prange(:, 1) = (/0., 0., 0./) ! Range of parameters (lower bound).
    prange(:, 2) = (/0.3, 5., 2./) ! Range of parameters (upper bound).
    varmin = 0.001 ! the terminating limit for the variance of the function values. varmin>0 is required
    step = (/0.15, 2.5, 1./) ! determines the size and shape of the initial simplex. The relative magnitudes of its elements should reflect the units of the variables. size(step)=size(start)
    konvge = 100 ! the convergence check is carried out every konvge iterations
    maxeval = 2000 ! the maximum number of function evaluations. default: 1000

    ! Call Nelder-Mead optimizer to reduce GCOMP
    ! pmin = nelmin(obj_func, pstart, varmin, step, konvge, maxeval, &
    !               funcmin, neval, numrestart, ierror, history)
    pmin = nelminrange(obj_func, pstart, prange, varmin, step, konvge, maxeval, &
                      funcmin, neval, numrestart, ierror)

    ! scale up distance h
    where (gamma(:,1) > 0._dp) gamma(:,1) = gamma(:,1) * gmax(1)

    ! scale back parameters (only for range a)
    pmin(3)=pmin(3)*gmax(1)
    beta = pmin

    print *, '=============================='
    print *, ' Results of Nelder-Mead optimization '
    print *, neval, ' of ', maxeval
    print *, "funcmin: ", funcmin
    print *, "p_obj:   ", pmin
    print *, 'error: ', ierror
    print *, 'varmin: ', varmin
    if (allocated(history)) print *, 'history: ', history(1), history(size(history))
    print *, 'gmax: ', gmax(1)


    ! estimate statistics
    call stats
    !
    ! scale up distance h
    where (gamma(:,1) > 0._dp) gamma(:,1) = gamma(:,1) * gmax(1)

  end subroutine OPTI

  !> \brief   Function to be minimised for the nelder mead algorithm
  !> \return  Target function value for given parameter set.
  function obj_func(p)
    use varfit, only      : nbins, gamma, nh
    implicit none
    real(dp), dimension(:), intent(in) :: p
    real(dp) :: obj_func
    real(dp) :: gcal
    integer(i4) :: k
    !
    obj_func = 0.0_dp
    !
    do k=1,nbins
      if (gamma(k,1) > 0.0_dp) then
        ! gcal = tvar (gamma(k,1), xu(1), xu(2), xu(3))
        gcal = tvar(gamma(k,1), p(1), p(2), p(3))
        !
        ! estimator l1 weighted
        if (gamma(k,1) <= 1.0_dp ) obj_func = obj_func + dabs(gamma(k,2) - gcal) * dble(Nh(k))
      end if
    end do
  end function obj_func

  !> \brief   statistics
  subroutine stats
    use varFit, only                     : E, beta, gamma
    implicit none
    integer(i4), parameter               :: incx = 1
    integer(i4)                          :: k
    integer(i4)                          :: ne
    real(dp)                             :: SSE
    real(dp)                             :: zObsMean, zCalMean
    real(dp)                             :: zObsVar,  zCalVar, sumP, NSE_denom
    real(dp), dimension(:), allocatable  :: error, denom
    real(dp), dimension(:), allocatable  :: zCal, zObs
    real(dp), parameter                  :: small = -9.999d3

    !
    !Initialize
    ne = count(gamma(:,1) > 0.0_dp)
    allocate (error(ne), denom(ne), zCal(ne), zObs(ne))
    zObs = gamma(1:ne,2)
    zCal = 0.0_dp
    do k=1,ne
      if (gamma(k,1) > 0._dp)  zCal(k) = tvar(gamma(k,1),beta(1),beta(2),beta(3))
    end do
    !
    error = zObs-zCal
    if ( ne > 1 ) then
      zObsMean  = sum(zObs)/real(ne, dp)
      zCalMean  = sum(zCal)/real(ne, dp)
      sumP      = dot_product(zObs,zCal)
      zObsVar   = dot_product(zObs,zObs) - real(ne, dp) * zObsMean * zObsMean
      zCalVar   = dot_product(zCal,zCal) - real(ne, dp) * zCalMean * zCalMean
      SSE       = dot_product(error,error)
      denom     = zObs - zObsMean
      NSE_denom = dot_product(denom,denom)
    else
    zObsMean = small
    zCalMean = small
    zObsVar  = small
    zCalVar  = small
    end if
    !	****************
    !	Quality measures
    !	****************
    if ( ne > 0 ) then
      !
      ! BIAS
      E(1) = zCalMean - zObsMean
      print *, 'BIAS: ', E(1)
      !
      ! MSE
      E(2) = SSE/real(ne, dp)
      print *, 'MSE:  ', E(2)
      !
      ! RMSE
      if ( E(2) > 0.0_dp ) then
        E(3) = dsqrt(E(2))
      else
        E(3) = small
      end if
      print *, 'RMSE: ', E(3)
      !
      ! RRMSE
      if ( E(3) > 0.0_dp ) then
        E(4)=E(3)/zObsMean
      else
        E(4)= small
      end if
      print *, 'PRMSE:', E(4)
      !
      ! MAE
      E(5)= sum(abs(error))
      E(5)= E(5)/real(ne, dp)
      print *, 'MAE:  ', E(5)
      !
      ! RMAE
      E(6)=E(5)/zObsMean
      print *, 'RMAE: ', E(6)
      !
      ! r
      E(7)= (sumP-real(ne, dp) * zCalMean * zObsMean) / dsqrt(zCalVar * zObsVar)
      print *, 'r:    ', E(7)
      !
      ! NSE
      E(8)= 1.0_dp - (SSE/NSE_denom)
      print *, 'NSE:  ', E(8)
    else
      E = small
    end if
    deallocate (error, denom, zCal, zObs)
  end subroutine stats

end module mo_edk_setvario
