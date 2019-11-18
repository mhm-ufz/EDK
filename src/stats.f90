!	**************************************************************************
!
! SUBROUTINE      statistics
!
!	**************************************************************************
subroutine stats
  use varFit, only                     : E, beta, gamma, nbins
  use mo_kind, only                    : i4, dp
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
  real(dp)                             :: tVar

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

