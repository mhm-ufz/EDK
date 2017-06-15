!******************************************************************************
! Function                                 To be optimised:
! VARIOGRAM TYPES:                         sum of errors=G(1)=> min!
!                 pure nugget + potential
!     
!******************************************************************************
subroutine GCOMP(G,XU)
  use mo_kind, only  : i4, dp
  use VarFit, only   : vType, nbins, gamma, Nh
  real(dp)              G(10),XU(100)
  real(dp)              gCal, tvar
  integer(i4)           k 
  !
  G(1) = 0.0_dp    
  !
  do k=1,nbins
    if (gamma(k,1) > 0.0_dp) then
      gCal = tvar (gamma(k,1), XU(1), XU(2), XU(3))
      !
      ! Estimator L1 weighted
      if (gamma(k,1) <= 1.0_dp ) G(1)=G(1)+dabs(gamma(k,2) - gCal)*dble(Nh(k))
    end if
  end do
end
