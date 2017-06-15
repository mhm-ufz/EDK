!****************************************************************************
!
!  SUBROUTINE: clean
!
!  PURPOSE:  free allocated space of all arrays
!  UPDATES
!            Created        L. Samaniego            10.08.2010
!            Last Update
!****************************************************************************
subroutine clean
  use mainVar
  use mo_kind, only: i4
  use kriging
  use runControl
  integer(i4)     :: i
  !
  ! DEM will be reused...
  !
  
  ! Stations
  do i = 1, nSta
    if ( allocated( MetSta(i)%z ) ) deallocate( MetSta(i)%z )
  end do
  if ( allocated(MetSta)  ) deallocate (MetSta)

  do i=1,nSta-1
    if ( allocated( dS(i)%S ) ) deallocate( dS(i)%S )
  end do  
  if ( allocated(dS)  ) deallocate (dS)

  do i=1,nCell
    if ( allocated( cell(i)%listNS ) )  deallocate ( cell(i)%listNS )
  end do
  if ( allocated(cell)) deallocate (cell)
  if ( allocated(dCS) ) deallocate (dCS)
  if ( allocated(dz2S)) deallocate (dz2S)

end subroutine clean



