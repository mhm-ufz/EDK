!****************************************************************************
!
!  SUBROUTINE: Ordinary Kriging
!
!  PURPOSE:  Perform OK (daily)
!            Output: output gridsize = cellFactor * gridsize DEM
!
!            Created        L. Samaniego            14.04.2006
!            Last Update       Sa                   14.04.2006
!
!****************************************************************************
subroutine OK(jd,k)
  use mainVar
  use kriging
  use runControl
  use varfit, only      : beta
  use mo_kind, only     : i4, dp
  use mo_setVario, only : tVar
  ! use LSASF_INT
  ! use lapack95, only: gesv
  implicit none
  integer(i4), intent(in)         :: jd                   ! day
  integer(i4), intent(in)         :: k                    ! cell id
  integer(i4)                     :: i, j, l, ll
  integer(i4)                     :: ii, jj
  integer(i4)                     :: Nk(nSta), nNmax 
  real(dp), allocatable           :: A (:,:), B(:), X(:)
  real(dp), allocatable           :: lamda(:)
  real(dp)                        :: sumLamda
  !
  ! check DEM
  if (cell(k)%h == grid%nodata_value ) then
     cell(k)%z = gridMeteo%noData_value
     return
  end if 
  !
  ! Check which stations have valid data ( =/ -9, 999 etc. ) store them in Nk(:)
  l  = 0
  do i= 1, cell(k)%nNS
    j = cell(k)%listNS(i)
    if ( MetSta(j)%z(jd) /= noDataValue ) then
      l = l + 1
      Nk(l) = j
    end if
  end do
  nNmax = l
  !
  !>>>>  no value ! avoid indetermination
  ! avoid 0 value calculations ll == nNmax
  ! avoid calculations where only 1 station is available
  ! avoid numerical instabilities nNmax == 2 (may happen that the solver matrix becomes singular) 
  if (.not. ( ll == nNmax .or.  nNmax == 1 .or. nNmax == 2 ) ) then 
     ! initialize matrices  
     allocate (A(nNmax+1,nNmax+1), B(nNmax+1), X(nNmax+1))
     !
     ! assembling the system of equations: OK
     A=0.0_dp
     do i=1,nNmax-1
        ii=Nk(i)
        do j=i+1,nNmax
           jj=Nk(j)
           ! available only the upper triangular matrix
           if (jj > ii) then
              A(i,j)=tVar(dS(ii)%S(jj),beta(1),beta(2),beta(3))
           else
              A(i,j)=tVar(dS(jj)%S(ii),beta(1),beta(2),beta(3))
           end if
        end do
        A(i,nNmax+1) = 1.0_dp
        B(i)=tVar(dCS(k,ii),beta(1),beta(2),beta(3))
     end do
     A(nNmax,nNmax+1) = 1.0_dp
     ii=Nk(nNmax)
     B(nNmax)=tVar(dCS(k,ii),beta(1),beta(2),beta(3))
     !
     B(nNmax+1) = 1.0_dp
     !NOTE: only the upper triangular matrix is needed!
     ! call D_LSASF (A, B, X)
     call dgesv(A, B)
     X = B
     !
     ! The BLUE of z is then:
     cell(k)%z = 0.
     do i=1,nNmax
        ii=Nk(i)
        cell(k)%z = cell(k)%z + real(X(i) * MetSta(ii)%z(jd))
     end do
    !
    ! only one station available /= 0
  else if (ll /= nNmax .AND. nNmax == 1) then
    ii = Nk(1)
    cell(k)%z = MetSta(ii)%z(jd)
    ! 
    ! for precipitation, distant values are set to zero
    if (distZero .and. dCS(k,ii) .gt. thresholdDist) then
      cell(k)%z = 0.0_dp
    end if
    !  
    ! if all stations have the value 0

 ! if all stations have the value 0
  else if (ll == nNmax ) then 
     cell(k)%z = 0.0_dp

  ! avoid numerical instabilities --> IWD insverse weighted squared distance
  ! matrix of DLSASF may become instable if only two or three stations are available 
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
     
end subroutine OK


