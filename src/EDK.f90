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
subroutine EDK(jd,k)
  use mo_kind, only                : i4, dp
  use mainVar
  use kriging
  use runControl
  use varfit, only                 : beta
  ! use LSASF_INT
  ! use mkl95_lapack, only: gesv
  ! use lapack95, only: gesv

  implicit none
  integer(i4), intent(in)         :: jd                  ! julian day
  integer(i4), intent(in)         :: k                   ! cell id
  integer(i4)                     :: i, j, l, ll, info
  integer(i4)                     :: ii, jj
  integer(i4)                     :: Nk(nSta), nNmax
  integer(i4), allocatable        :: ipvt(:)
  real(dp), allocatable           :: A (:,:), B(:), X(:)
  real(dp)                        :: tVar
  real(dp), allocatable           :: lamda(:)
  real(dp)                        :: sumLamda
  !
  ! check DEM
  print *, 'ha...'
  if (nint(cell(k)%h) == grid%nodata_value ) then
     cell(k)%z = gridMeteo%nodata_value
     return
  end if 
  !
  ! Check which stations have valid data ( =/ -9, 999 etc. ) store them in Nk(:)
  l  = 0
  ll = 0

  do i= 1, cell(k)%nNS
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
    do i = 1, nNmax-1
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
   
    A(nNmax,nNmax+1) = 1.0_dp
    ii=Nk(nNmax)
    A(nNmax,nNmax+2) = MetSta(ii)%h


    B(nNmax)=tVar(dCS(k,ii),beta(1),beta(2),beta(3))
    !
    B(nNmax+1) = 1.0_dp
    B(nNmax+2) = cell(k)%h

    print *, 'hu...'
    ! NOTE: only the upper triangular matrix is needed!
    ! call D_LSASF (A, B, X)
    ! call gesv(A, B) ! EVE CentOS 7
    ! call dgesv(A, B) ! MacOS Sierra: Accelerate Framework
    call dgesv(size(A, 1), 1, A, size(A, 1), ipvt, B, size(A, 1), info)
    X = B
    print *, 'hu...'
    print *, 'info: ', info
    stop 'testing'
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
    deallocate (A,B,X)
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
  !
end subroutine EDK
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
        tVar = 0.0_dp
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
  integer(i4)                     :: NoCellsFiner
  real(dp)                        :: xc, yc
  integer(i4), allocatable        :: list(:)
  !
  ! Initialize variables
  if ( allocated(dCS) ) deallocate (dCS)
  if ( allocated(dS)  ) deallocate (dS)
  if ( allocated(cell)) deallocate (cell)
  if ( allocated(dz2S)) deallocate (dz2S)

  !
  allocate ( dz2S(nSta-1) )
  allocate ( dCS(nCell,nSta) )
  allocate ( dS(nSta-1)      )
  allocate ( cell(nCell)     )
  allocate ( list(nSta)      )
  
  !
  do i=1,nSta-1
    allocate ( dS(i)%S(i+1:nSta) )
    allocate ( dz2S(i)%S(i+1:nSta) )
    ! distance matrix between stations:            checked OK
    do j=i+1, nSta
      dS(i)%S(j) = dsqrt( ( MetSta(i)%x - MetSta(j)%x )**2 +  ( MetSta(i)%y - MetSta(j)%y )**2 )
      if (dS(i)%S(j) == 0.0_dp) then
        print* , 'Stations: ', MetSta(i)%Id, MetSta(j)%Id, ' have the same coordinates, or are repeated. Check LUT.'
        !stop
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
  NoCellsFiner = (gridMeteo%cellsize / grid%cellsize) **2 
  xc = gridMeteo%xllcorner + dble(gridMeteo%cellsize) * 0.5_dp
  delta = cellFactor / 2
  jj = delta
  do k=1,nCell 
 
    if (r == 1) then
       c = c + 1
       if (c > 1) then 
         xc = xc + dble(gridMeteo%cellsize)
         jj = jj + cellFactor
       end if
       yc = gridMeteo%yllcorner + dble(gridMeteo%cellsize) * (dble(gridMeteo%nrows) - 0.5_dp)
       ii = delta 
    else
       yc = yc - dble(gridMeteo%cellsize)
       ii = ii + cellFactor
    end if
    cell(k)%x = xc
    cell(k)%y = yc

    ! average of only four DEM cells around centre cell (from lower grid scale upto higher grid cell)
    !cell(k)%h = 0.25_dp*(G(ii,jj)%h + G(ii,jj+1)%h + G(ii+1,jj)%h + G(ii+1,jj+1)%h)
    !
    ! average of all DEM cells (from lower grid scale upto higher grid cell)
    nTcell =  count(G( (ii-delta+1):(ii+delta) , (jj-delta+1):(jj+delta) )%h  > grid%nodata_value )
    if (nTcell == 0) then
      cell(k)%h = gridMeteo%nodata_value
    else
      cell(k)%h  = sum(G( (ii-delta+1):(ii+delta) , (jj-delta+1):(jj+delta) )%h, &
                       G( (ii-delta+1):(ii+delta) , (jj-delta+1):(jj+delta) )%h /= gridMeteo%nodata_value ) / dble(nTcell)
      if ( NoCellsFiner < NoCellsFiner) then
         print*, 'Cells with non matching DEM on the finer scale (nodata values in finer grid)'
      end if
    end if
    r=r+1
    if (r > gridMeteo%nrows) r = 1
    ! MZMZMZMZ - delete !!!!!!!
    if (cellFactor == 1) then
       !print*, ii, jj
       cell(k)%h = G(ii+1,jj+1)%h
       cycle
    end if
    ! MZMZMZMZ - delete !!!!!!!
  end do

  ! distance matrix cell to stations: checked OK
  do j=1, nSta
    do i=1,nCell
      dCS(i,j) = dsqrt( ( cell(i)%x - MetSta(j)%x )**2 + ( cell(i)%y - MetSta(j)%y )** 2)
    end do
  end do

  ! find the closest stations to cell i (any order): checked  OK
  do i=1,nCell
    list = -9
    do j=1,nSta
      if (dCS(i,j) <= maxDist) list(j) = j
    end do
    cell(i)%nNS = count(list > -9)
    allocate ( cell(i)%listNS( cell(i)%nNS ) )
    cell(i)%listNS = pack(list, MASK = list >-9)
  end do
  !
  deallocate (list)

end subroutine dMatrix

