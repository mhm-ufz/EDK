!> \file    mo_edk_empvar.f90
!> \copydoc mo_edk_empvar

!> \brief   This program calculates the empirical variogram
!> \details gamma(:,1) : distance
!!          gamme(:,2) : semivarance
!!          botm are normalized with gmax(:)
!> \author  Luis Samaniego
!> \date    22.11.2001
!> \date    19.02.2004
module mo_edk_empvar

  implicit none

  private

  public :: EmpVar

contains

  !> \brief   calculate the empirical variogram
  subroutine EmpVar(jd, flagMax)
    use mainVar
    use mo_kind   , only       : i4, dp
    use kriging
    use VarFit
    implicit none
    integer(i4), intent(in)   :: jd
    logical, intent(in)       :: flagMax
    integer(i4)               :: i, j, k, nhk, ni
    real(dp)                  :: gk, delta

    ! variance h=0 (one pass algorithm)
    do i = 1,nSta
      if ( MetSta(i)%z(jd) ==  noDataValue ) cycle
      nobs  = nobs + 1
      delta = MetSta(i)%z(jd) - m0
      m0    = m0 + delta / real(nobs, dp)
      v0    = v0 + delta * (MetSta(i)%z(jd) - m0)
      !write (999,'(4i6, 3f15.4)') y,d,i, nobs, MetSta(i)%z(jd), m0, v0
    end do
    !
    !
    ! calculate the empiric variogram
    !
    ! selecting pair - bins
    if (jd == jStart) then
      nbins = ceiling(hMax/dh)
      if (.not. allocated (nH))   allocate (Nh(nbins),      gamma(nbins,2))
      Nh = 0
      gamma = 0.0_dp
    end if
    !
    print *, '***WARNING: Removal of outliers in the estimation of the variogram is activated'
    do i=1,nSta-1
      do j=i+1,nSta
        if (edk_dist%getZ2(i,j,jd) /=  noDataValue ) then
          ! take values up to max distance
          if (edk_dist%getSS(i,j) > hMax ) cycle
          ! ! remove outliers for the estimation of the variogram
          ! if (flagVarTyp == 2 ) then
            if (  dabs( MetSta(i)%h - MetSta(j)%h ) / edk_dist%getSS(i,j)  > gradE ) then
              ! write(999,*), 'pair removed', MetSta(i)%id, MetSta(j)%id
                cycle
            end if
          ! end if
          !
          k=max(1, ceiling(edk_dist%getSS(i,j)/dh))
          Nh(k)=Nh(k)+1
          gamma(k,2)=gamma(k,2) + edk_dist%getZ2(i,j,jd)
        end if
      end do
    end do
    !
    ! only at the end
    if (jd == jEnd) then
      !
      ! consolidate bins N(h)>=30
      !forward
      do k=1,nbins-1
        nhk=Nh(k)
        gk=gamma(k,2)
        if ((nhk>0) .and. (nhk<30)) then
          Nh(k)=-9
          gamma(k,2)=0.0_dp
          if (Nh(k+1)>0) then
            Nh(k+1)=Nh(k+1)+nhk
            gamma(k+1,2)=gamma(k+1,2)+gk
          else
            Nh(k+1)=nhk
            gamma(k+1,2)=gk
          end if
        end if
      end do
      !backward
      do k=nbins,2,-1
        nhk=Nh(k)
        gk=gamma(k,2)
        if ((nhk>0) .and. (nhk<30)) then
          Nh(k)=-9
          gamma(k,2)=0_dp
          if (Nh(k-1)>0) then
            Nh(k-1)=Nh(k-1)+nhk
            gamma(k-1,2)=gamma(k-1,2)+gk
          else
            Nh(k-1)=nhk
            gamma(k-1,2)=gk
          end if
        end if
      end do
      !
      ni=0
      ! estimate semi-variance gamma(i,1)=h, gamma(i,2)=g(h)
      do k=1,nbins
        if (Nh(k) > 0) then
          ni=ni+1
          ! Classsical
          gamma(k,2)=gamma(k,2)/2.0_dp/real(Nh(k), dp)
          !
          ! Cressi-Hawkins: adjust bias
          !gamma(k,2)=0.5_dp/(0.457_dp+0.494_dp/dfloat(Nh(k)))*(gamma(k,2)/dfloat(Nh(k)))**4
          !
          gamma(k,1)=real(k, dp)*dh-real(ni, dp)*dh/2.0_dp
          ni=0
        else
          ni=ni+1
          gamma(k,1)=-9.0_dp
        end if
      end do
      !
      ! scaling
      if (flagMax) gmax=maxval(gamma, dim=1)

      do k=1,nbins
        if (gamma(k,1) > 0.0_dp) then
          gamma(k,1) = gamma(k,1)/gmax(1)
          gamma(k,2) = gamma(k,2)/gmax(2)
        end if
      end do
      !
      !keep variogram
    end if

  end subroutine EmpVar

end module mo_edk_empvar
