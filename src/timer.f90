!**********************************************************************************
!  TIMER
!
!  PURPOSE:  To estimate the running time of the simulation. It calls the system's
!            function DTIME().
!  UPDATES:
!            Created      Sa   24.01.2003 
!            Last update  Sa   24.01.2003
!**********************************************************************************
subroutine Timer
  use RunControl
  use mo_kind, only : i4
  implicit none
  real(i4)          :: DTIME
  if (RunTime <= 0.0) then
      RunTime = DTIME(TA)
  else
    open (unit=500, file='timer.out', status='unknown')
      RunTime = DTIME(TA)
      write (500, 100) RunTime, TA(1), TA(2)
     close (500, status='keep')
  end if


  ! formats
  100  format (  '-------------------------------------',&
                /'        End of the Simulation        ',&
                /'-------------------------------------',&
                /' Elapsed CPU time:    ', F10.1, ' s'  ,&
                /' Elapsed user time:   ', F10.1, ' s'  ,&
                /' Elapsed system time: ', F10.1, ' s'  ,&
                /'-------------------------------------')
end subroutine Timer

