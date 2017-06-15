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
subroutine setVario
  use runControl
  use VarFit
  use mainVar
  use mo_kind, only: i4, dp
  implicit none
  integer(i4)               :: jd, y, y0
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
     v0 = v0 / dfloat(nobs)
     ! optimize parameters
     open(UNIT = 6,FILE = 'Report_OPT.sol', STATUS='REPLACE')
     call opti

  end if
end subroutine setVario
