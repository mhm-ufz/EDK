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
subroutine EDK(jd, k, dCS, MetSta, dS, cell)
  use mo_kind, only                : i4, dp
  use mainVar, only                : MeteoStation, noDataValue, nSta, thresholdDist
  use kriging, only                : CellCoarser, dtoS
  use runControl, only             : flagVarTyp
  use varfit, only                 : beta
  use mo_setVario, only            : tVar
  implicit none

  ! input / output variables
  integer(i4), intent(in)          :: jd        ! julian day
  integer(i4), intent(in)          :: k         ! cell id
  real(dp), intent(in)             :: dCS(:, :) ! distances matrix
  type(MeteoStation), intent(in)   :: MetSta(:) ! MeteoStation input
  type(dtoS), intent(in)           :: dS(:)     ! distance among stations
  type(CellCoarser), intent(inout) :: cell(:)   ! cell specification

  ! local variables
  integer(i4)                     :: i, j, l, ll, info
  integer(i4)                     :: ii, jj
  integer(i4)                     :: Nk(nSta), nNmax
  integer(i4), allocatable        :: ipvt(:)
  real(dp), allocatable           :: A (:,:), B(:), X(:) !, C(:,:)
  real(dp), allocatable           :: lamda(:)
  real(dp)                        :: sumLamda
  !
  ! Check which stations have valid data ( =/ -9, 999 etc. ) store them in Nk(:)
  l  = 0
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

  !>>>>  no value ! avoid indetermination
  ! avoid 0 value calculations ll == nNmax
  ! avoid calculations where only 1 station is available
  ! avoid numerical instabilities nNmax == 2 (may happen that the solver matrix becomes singular)
  if (.not. ( ll == nNmax .or.  nNmax == 1 .or. nNmax == 2 ) ) then     
    !
    ! initialize matrices  
    allocate (A(nNmax+2,nNmax+2), B(nNmax+2), X(nNmax+2), ipvt(nNmax + 2))
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
        else
          A(i,j)=tVar(dS(jj)%S(ii),beta(1),beta(2),beta(3))
          A(j,i)=A(i,j)
        end if
      end do
      A(i,i) = tVar(0._dp,beta(1),beta(2),beta(3))

      A(i,nNmax+1) = 1.0_dp
      A(nNmax+1,i) = 1.0_dp
      A(i,nNmax+2) = MetSta(ii)%h
      A(nNmax+2,i) = MetSta(ii)%h

      B(i)=tVar(dCS(k,ii),beta(1),beta(2),beta(3))
   end do
   
    !
    B(nNmax+1) = 1.0_dp
    B(nNmax+2) = cell(k)%h

    ! NOTE: only the upper triangular matrix is needed!
    ! call D_LSASF (A, B, X)
    ! print *, '<<<<<<<<<<<<<<'
    ! print *, 'rhs = ', B
    ! C = A
    ! print *, 'A = ', C(:, 1)
    ! print *, '<<<<<<<<<<<<<<'
    call dgesv(size(A, 1), 1, A, size(A, 1), ipvt, B, size(A, 1), info)
    ! print *, '<<<<<<<<<<<<<<'
    ! print *, 'B = ', B
    ! print *, 'A = ', A(:, 1)
    ! print *, '<<<<<<<<<<<<<<'
    ! print *, 'result: ', sum(C(:, 1) * B)
    X = B
    if (info .ne. 0_i4) then
      print *, '***WARNING: calculation of weights failed'
    end if
    if (abs(sum(X(:nNmax)) - 1._dp) .gt. 1.e-4) then
      print *, '***WARNING: sum of weights is not 1, calculation of weights failed'
      print *, 'sum of weights: ', sum(X(:nNmax))
    end if
    ! print *, 'number of neighbors: ', nNmax
    ! print *, 
    ! print *, 'ipvt: ', ipvt
    ! print *, 'info: ', info
    ! stop 'testing'
    !
    ! The BLUE of z is then:
    cell(k)%z = 0.
    do i=1,nNmax
      ii=Nk(i)
      cell(k)%z = cell(k)%z + X(i) * MetSta(ii)%z(jd)
    end do
    ! correct in case of negative values
    ! sick system of equations caused by lack of 
    ! precipitation values in nearby stations
    if (cell(k)%z < 0.0_dp .and. flagVarTyp == 1) then
      cell(k)%z = 0.0_dp
    end if
    deallocate (A,B,X,ipvt)
    !
    ! only one precipiation station available /= 0
  else if (ll /= nNmax .AND. nNmax == 1 .AND. flagVarTyp == 1) then
     ii = Nk(1)
    ! problably convective rain at the only station
     if (dCS(k,ii) <=  thresholdDist) then
        cell(k)%z = MetSta(ii)%z(jd)
     else
        cell(k)%z = 0.0_dp
     end if

  ! only one station available, which has values /= 0
  else if (ll /= nNmax .AND. nNmax == 1 .AND. flagVarTyp /= 1) then
     ii = Nk(1)
     cell(k)%z = MetSta(ii)%z(jd)

 ! if all stations have the value 0
  else if (ll == nNmax ) then 
     cell(k)%z = 0.0_dp

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
        cell(k)%z = cell(k)%z + lamda(i) * MetSta(ii)%z(jd)
     end do
     deallocate (lamda)
  end if
  !
end subroutine EDK

end module mo_EDK
