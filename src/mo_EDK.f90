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
  subroutine EDK(k, jStart, jEnd, dCS, MetSta, dS, cell, X, Nk_old, doOK)
    use mo_kind,    only : i4, dp, sp
    use mainVar,    only : MeteoStation, noDataValue, nSta, thresholdDist
    use kriging,    only : CellCoarser, dtoS
    use runControl, only : correctNeg, distZero

    implicit none

    ! input / output variables
    integer(i4),                intent(in)    :: k         ! cell id
    integer(i4),                intent(in)    :: jStart, jEnd
    integer(i4),                intent(inout) :: Nk_old(nSta) ! added Nk_old(nSta)
    real(dp),                   intent(in)    :: dCS(:, :) ! distances matrix
    real(dp), allocatable,      intent(inout) :: X(:)      ! added X(:) 
    type(MeteoStation),         intent(in)    :: MetSta(:) ! MeteoStation input
    type(dtoS),                 intent(in)    :: dS(:)     ! distance among stations
    type(CellCoarser),          intent(inout) :: cell(:)    ! cell specification
    logical, optional,          intent(in)    :: doOK   ! switch do ordinary kriging

    ! local variables
    logical                         :: doOK_loc, calc_weights
    integer(i4)                     :: jd        ! julian day
    integer(i4)                     :: i, j, l, ll
    integer(i4)                     :: ii, jj
    integer(i4)                     :: Nk(nSta),nNmax !Nk_old !Nk(nSta), Nk_old(nSta), nNmax ! deleted: Nk_old(nSta)   
    real(dp), allocatable           :: A (:,:), B(:), C(:,:) ! deleted: X(:)  
    real(dp), allocatable           :: lamda(:)
    real(dp)                        :: sumLamda

    !
    ! Check stations have valid data ( =/ -9, 999 etc. ) store them in Nk(:)
    l      = 0
    ll     = 0
    Nk     = 0
    ! switch ordinary kriging off if not explicitly given
    doOK_loc = .False.
    if (present(doOK)) doOK_loc = doOK
    ! IF NK changed -> re-estimate weights
    timeloop: do jd = jStart, jEnd
      if (jd > jStart) Nk_old = Nk
      Nk = 0_i4
      l = 0
      ll = 0
      do i = 1, cell(k)%nNS
         j = cell(k)%listNS(i)
        if ( MetSta(j)%z(jd) /= noDataValue ) then

          !***      zero field     ***********************************
          if (MetSta(j)%z(jd) == 0.0_dp ) ll = ll + 1
          !**********************************************************
          l = l + 1
          Nk(l) = j
        end if
      end do
      nNmax = l
      if (all(Nk == Nk_old) .and. allocated(X)) then
        calc_weights = .False.
      else
        calc_weights = .True.
      end if
      !>>>>  no value ! avoid indetermination
      ! avoid 0 value calculations ll == nNmax
      ! avoid calculations where only 1 station is available
      ! avoid numerical instabilities nNmax == 2 (may happen that the solver matrix becomes singular)
      if (.not. ( ll == nNmax .or.  nNmax == 1 .or. nNmax == 2 ) ) then     

        if (calc_weights) then
          ! print *, 'cell:       ', k
           !print *, '#neighbors: ', nNmax
           !print *, 'Nk:     ', Nk
           !print *, 'Nk_old: ', Nk_old
           !write(*,*),"calc_weights inside = ",calc_weights
          call get_kriging_weights(X, nNmax, Nk, doOK_loc, dCS(k, :), dS, cell(k), MetSta)
          !write(*,*),"X after kriging weights: ",X
        end if
        ! The BLUE of z is then:
        cell(k)%z(jd) = 0.
        do i=1,nNmax
          ii=Nk(i)
          cell(k)%z(jd) = cell(k)%z(jd) + X(i) * MetSta(ii)%z(jd)
        end do
        !
        ! only one station available /= 0
      else if (ll /= nNmax .AND. nNmax == 1) then
        ii = Nk(1)
        cell(k)%z(jd) = MetSta(ii)%z(jd)
        ! 
        ! for precipitation, distant values are set to zero
        if (distZero .and. dCS(k,ii) .gt. thresholdDist) then
          cell(k)%z(jd) = 0.0_dp
        end if
        !  
        ! if all stations have the value 0
      else if (ll == nNmax ) then 
        cell(k)%z(jd) = 0.0_dp
        ! avoid numerical instabilities --> IWD insverse weighted squared distance
        ! matrix of solver may become instable if only two or three stations are available 
      else if (ll /= nNmax .AND. nNmax == 2) then
        allocate (lamda(nNmax))
        do i=1,nNmax
          ii=Nk(i)
          lamda(i)=1/dCS(k,ii)/dCS(k,ii)
        end do
        sumLamda=sum(lamda)
        lamda=lamda/sum(lamda)
        cell(k)%z = 0.0_dp
        do i=1,nNmax
          ii=Nk(i)
          cell(k)%z(jd) = cell(k)%z(jd) + lamda(i) * MetSta(ii)%z(jd)
        end do
        deallocate (lamda)
      end if
      ! stop 'TESTING'
    end do timeloop
    !if (allocated(X)) deallocate (X)
    !
    ! correct negative
    if (correctNeg) cell(k)%z = merge(0._sp, cell(k)%z, (cell(k)%z .gt. -9999._sp) .and. (cell(k)%z .lt. 0.))
    !
  end subroutine EDK

  subroutine get_kriging_weights(X, nNmax, Nk, doOK_loc, dCS, dS, cell, MetSta)

    use mo_kind,     only : dp, i4
    use mainVar,     only : MeteoStation, noDataValue, nSta, thresholdDist
    use kriging,     only : CellCoarser, dtoS
    use varfit,      only : beta
    use mo_setVario, only : tVar

    implicit none

    real(dp), allocatable, intent(out) :: X(:)
    integer(i4),           intent(in)  :: nNmax
    integer(i4),           intent(in)  :: Nk(nSta)!Nk(nSta)
    logical,               intent(in)  :: doOK_loc
    real(dp),              intent(in)  :: dCS(:)    ! distances matrix
    type(dtoS),            intent(in)  :: dS(:)     ! distance among stations
    type(CellCoarser),     intent(in)  :: cell      ! cell specification
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

    ! loop over available stations nNmax
    do i = 1, nNmax
      ii = Nk(i)

      do j = i + 1, nNmax
        jj = Nk(j)

        ! available only the upper triangular matrix
        if (jj > ii) then
          A(i,j)=tVar(dS(ii)%S(jj),beta(1),beta(2),beta(3))
          A(j,i)=A(i,j)
          ! print *, 'A: ', A(i,j), dS(ii)%S(jj),beta(1),beta(2),beta(3)
        else
          A(i,j)=tVar(dS(jj)%S(ii),beta(1),beta(2),beta(3))
          A(j,i)=A(i,j)
        end if
      end do
      A(i,i) = tVar(0._dp, beta(1),beta(2),beta(3))

      A(i,nNmax+1) = 1.0_dp
      A(nNmax+1,i) = 1.0_dp
      if (.not. doOK_loc) A(i,nNmax+2) = MetSta(ii)%h
      if (.not. doOK_loc) A(nNmax+2,i) = MetSta(ii)%h

      B(i)=tVar(dCS(ii), beta(1),beta(2),beta(3))
      ! print *, 'B: ', B(i), dCS(k,ii), beta(1),beta(2),beta(3)
    end do

    !
    B(nNmax+1) = 1.0_dp
    if (.not. doOK_loc) B(nNmax+2) = cell%h

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
    if (abs(sum(X(:nNmax)) - 1._dp) .gt. 1.e-4) then
      print *, '***WARNING: sum of weights is not 1, calculation of weights failed'
      print *, 'sum of weights: ', sum(X(:nNmax))
    end if
    ! print *, 'easting: ', cell%x
    ! print *, 'northing: ', cell%y
    ! print *, 'number of neighbors: ', nNmax
    ! print *, ''
    ! ! print *, 'ipvt: ', ipvt
    ! print *, 'info: ', info
    ! stop 'testing'
    deallocate (A, B, C, ipvt)
  end subroutine get_kriging_weights

end module mo_EDK
