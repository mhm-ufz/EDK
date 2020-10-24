module mo_EDK

  implicit none

  private

  public :: EDK

contains
  !****************************************************************************
  !
  !  SUBROUTINE: External Drift Kriging
  !
  !  PURPOSE:  Perform EDK (daily)
  !            Output: output gridsize = cellFactor * gridsize DEM
  !  UPDATES
  !            Created        L. Samaniego            22.03.2006
  !            Last Update       Sa                   14.04.2006
  !            Last Update       Sa                   14.07.2010   DEM ckeck
  !            Last Update       Zi                   05.05.2011   IWD if No stations < 2
  !            Last Update       Zi                   07.02.2012   correct prec data < 0 --> = 0
  !****************************************************************************
  subroutine EDK(k, jStart, jEnd, edk_dist, MetSta, cell, doOK)
    use mo_kind,    only : i4, dp, sp
    use mainVar,    only : MeteoStation, noDataValue, thresholdDist
    use kriging,    only : CellCoarser
    use runControl, only : correctNeg, distZero
    use mo_edk_types, only: dist_t

    implicit none

    ! input / output variables
    integer(i4),                intent(in)    :: k         ! cell id
    integer(i4),                intent(in)    :: jStart, jEnd
    type(dist_t),               intent(in)    :: edk_dist  ! distances matrix
    type(MeteoStation),         intent(in)    :: MetSta(:) ! MeteoStation input
    type(CellCoarser),          intent(inout) :: cell    ! cell specification
    logical, optional,          intent(in)    :: doOK   ! switch do ordinary kriging

    ! local variables
    logical                         :: doOK_loc, calc_weights
    integer(i4)                     :: jd        ! julian day
    integer(i4)                     :: i, j
    integer(i4)                     :: n_select, n_zero
    integer(i4)                     :: ii, jj
    logical                         :: selectNS(cell%nNS) ! selected neighborhood (no no-data-values) current time-step
    logical                         :: selectNS_old(cell%nNS) ! selected neighborhood (no no-data-values) previous time-step
    integer(i4), allocatable        :: selectNS_ids(:)
    real(dp), allocatable           :: lamda(:)
    real(dp), allocatable           :: weights(:)
    real(dp)                        :: sumLamda

    selectNS_old = .false.
    ! switch ordinary kriging off if not explicitly given
    doOK_loc = .False.
    if (present(doOK)) doOK_loc = doOK
    ! IF NK changed -> re-estimate weights
    timeloop: do jd = jStart, jEnd
      if (jd > jStart) selectNS_old = selectNS
      selectNS = .false. ! reset selected neighborhood
      n_select = 0 ! counter for no-data-values in current Neighborhood
      n_zero = 0 ! counter for zero data values in current Neighborhood
      do i = 1, cell%nNS
         j = cell%listNS(i)
        if ( MetSta(j)%z(jd) /= noDataValue ) then

          !***      zero field     ***********************************
          if (MetSta(j)%z(jd) == 0.0_dp ) n_zero = n_zero + 1
          !**********************************************************
          n_select = n_select + 1
          selectNS(i) = .true.
        end if
      end do
      ! list of IDs of selected neighbors
      if ( allocated(selectNS_ids) ) deallocate(selectNS_ids)
      allocate(selectNS_ids(n_select))
      selectNS_ids = pack(cell%listNS, mask=selectNS)

      if (all(selectNS .eqv. selectNS_old) .and. allocated(weights)) then
        calc_weights = .False. ! same Neighborhood
      else
        calc_weights = .True.  ! new Neighborhood (new stations or old station with missing data)
      end if

      if (.not. ( n_zero == n_select .or.  n_select == 1 .or. n_select == 2 ) ) then
        !>>>>  no value ! avoid indetermination
        ! avoid 0 value calculations n_zero == n_select
        ! avoid calculations where only 1 station is available
        ! avoid numerical instabilities n_select == 2 (may happen that the solver matrix becomes singular)
        if (calc_weights) then
          if ( allocated(weights) ) deallocate(weights)
          call get_kriging_weights(weights, n_select, selectNS_ids, doOK_loc, edk_dist, k, cell%h, MetSta)
        end if
        ! The BLUE of z is then:
        cell%z(jd) = 0.
        do i=1, n_select
          ii = selectNS_ids(i)
          cell%z(jd) = cell%z(jd) + real(weights(i) * MetSta(ii)%z(jd), sp)
        end do

      else if (n_zero /= n_select .AND. n_select == 1) then
        ! only one station available /= 0
        ii = selectNS_ids(1)
        cell%z(jd) = real(MetSta(ii)%z(jd), sp)
        ! TODO: this is actually wrong (should be also calculated by distance)

        ! for precipitation, distant values are set to zero
        if (distZero .and. edk_dist%getCS(k,ii) .gt. thresholdDist) then
          cell%z(jd) = 0.0_sp
        end if

      else if (n_zero /= n_select .AND. n_select == 2) then
        ! avoid numerical instabilities --> IWD insverse weighted squared distance
        ! matrix of solver may become instable if only two or three stations are available
        allocate (lamda(n_select))
        do i=1, n_select
          ii=selectNS_ids(i)
          lamda(i)=1/edk_dist%getCS(k,ii)/edk_dist%getCS(k,ii)
        end do
        sumLamda=sum(lamda)
        lamda=lamda/sum(lamda)
        cell%z = 0.0_dp
        do i=1, n_select
          ii=selectNS_ids(i)
          cell%z(jd) = cell%z(jd) + real(lamda(i) * MetSta(ii)%z(jd), sp)
        end do
        deallocate (lamda)

      else if (n_zero == n_select) then  ! also for empty neighborhood
        ! if all stations have the value 0
        cell%z(jd) = 0.0_sp

      end if
    end do timeloop

    ! correct negative
    if (correctNeg) cell%z = merge(0._sp, cell%z, (cell%z .gt. -9999._sp) .and. (cell%z .lt. 0.))
    !
  end subroutine EDK

  subroutine get_kriging_weights(X, nNmax, Nk, doOK_loc, edk_dist, k, cell_h, MetSta)

    use mo_kind,     only : dp, i4
    use mainVar,     only : MeteoStation
    use varfit,      only : beta
    use mo_setVario, only : tVar
    use mo_edk_types, only : dist_t

    implicit none

    real(dp), allocatable, intent(out) :: X(:)
    integer(i4),           intent(in)  :: nNmax
    integer(i4),           intent(in)  :: Nk(:)
    logical,               intent(in)  :: doOK_loc
    type(dist_t),          intent(in)  :: edk_dist  ! distances
    integer(i4),           intent(in)  :: k         ! cell index
    real(dp),              intent(in)  :: cell_h    ! cell elevation
    type(MeteoStation),    intent(in)  :: MetSta(:) ! MeteoStation input

    ! local variables
    integer(i4)              :: i, j, ii, jj, info
    integer(i4), allocatable :: ipvt(:)
    real(dp),    allocatable :: A (:,:), B(:), C(:,:)

    !
    ! initialize matrices
    if (doOK_loc) then
      allocate (A(nNmax+1,nNmax+1), B(nNmax+1), X(nNmax+1), ipvt(nNmax + 1))
    else
      allocate (A(nNmax+2,nNmax+2), B(nNmax+2), X(nNmax+2), ipvt(nNmax + 2))
    end if
    !
    ! assembling the system of equations: OK
    A=0.0_dp

    ! loop over available stations nNmax in neighborhood
    do i = 1, nNmax
      ii = Nk(i)

      do j = i + 1, nNmax
        jj = Nk(j)

        ! available only the upper triangular matrix
        if (jj > ii) then
          A(i,j)=tVar(edk_dist%getSS(ii,jj),beta(1),beta(2),beta(3))
          A(j,i)=A(i,j)
        else
          A(i,j)=tVar(edk_dist%getSS(jj,ii),beta(1),beta(2),beta(3))
          A(j,i)=A(i,j)
        end if
      end do
      A(i,i) = tVar(0._dp, beta(1),beta(2),beta(3))

      A(i,nNmax+1) = 1.0_dp
      A(nNmax+1,i) = 1.0_dp
      if (.not. doOK_loc) A(i,nNmax+2) = MetSta(ii)%h
      if (.not. doOK_loc) A(nNmax+2,i) = MetSta(ii)%h

      B(i)=tVar(edk_dist%getCS(k,ii), beta(1),beta(2),beta(3))
    end do

    !
    B(nNmax+1) = 1.0_dp
    if (.not. doOK_loc) B(nNmax+2) = cell_h

    ! NOTE: only the upper triangular matrix is needed!
    ! call D_LSASF (A, B, X)
    ! print *, '<<<<<<<<<<<<<<'
    ! print *, 'rhs = ', B(:10)
    C = A
    X = B
    ! print *, 'A = ', C(:10, 1)
    ! print *, '<<<<<<<<<<<<<<'
    call dgesv(size(A, 1), 1, A, size(A, 1), ipvt, B, size(A, 1), info)
    ! print *, '<<<<<<<<<<<<<<'
    ! print *, 'B = ', B(:10)
    ! print *, 'A = ', A(:10, 1)
    ! print *, '<<<<<<<<<<<<<<'
    if (maxval(abs(matmul(C, B) - X)) .gt. 1e-10) print *, 'maximum error: ', maxval(abs(matmul(C, B) - X))
    X = B
    if (info .ne. 0_i4) then
      print *, '***WARNING: calculation of weights failed'
    end if
    ! print *, 'shape of A: ', shape(A)
    ! print *, 'dem of stations: ', A(:, 10)
    ! print *, 'maximum of dem at stations: ', maxval(MetSta(:)%h)
    ! print *, 'easting: ', cell%x
    ! print *, 'northing: ', cell%y
    ! print *, 'number of neighbors: ', nNmax
    ! print *, ''
    ! ! print *, 'ipvt: ', ipvt
    ! print *, 'info: ', info
    if (abs(sum(X(:nNmax)) - 1._dp) .gt. 1.e-4) then
       print *, '***WARNING: sum of weights is not 1, calculation of weights failed'
       print *, 'sum of weights: ', sum(X(:nNmax))
       ! stop 'testing'
    end if
    deallocate (A, B, C, ipvt)
  end subroutine get_kriging_weights

end module mo_EDK
