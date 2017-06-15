!*************************************************************************
!    PURPOSE    WRITE METEREOLOGIC VARIABLES
!    FORMAT     xxx_yyyy.bin
!                            
!               xxx   : | pre   flag = 1
!                       | tem          2
!                       | pet          3
!               yyyy  :   year
!               rec j :   grid   lenght (ncol x nrow x 4) bytes 
!                   j :   julian day
!
!               note      direct access, unformatted file
!                         with records sized for the whole grid
!
!    AUTHOR:    Luis E. Samaniego-Eguiguren, UFZ
!    UPDATES
!               Created        Sa   21.03.2006
!               Last Update    Sa   
!**************************************************************************
subroutine WriteDataMeteo(y,d,wFlag) 
  use mo_kind, only         : i4, sp, dp      
  use mainVar
  use runControl
  use kriging
  use VarFit
  !
  implicit none
  integer(i4), intent (in)  :: y, d, wFlag
  integer(i4)               :: i, j, k
  integer(i4)               :: leap             ! leap day either 0 or 1 
  character(256)            :: dummy
  character(256)            :: fileName
  logical                   :: wasOpened
  real(dp)                  :: tVar

  !
  select case (wFlag)
  case (1)
    !---------------------
    ! write Binary files
    !---------------------
    write (dummy, 110) y
    fileName = trim(dataPathOut)//trim(dummy)

    ! open file
    inquire (file=trim(fileName), opened = wasOpened)
    if (.not. wasOpened) then
       open (unit=100, file=fileName, form='unformatted', access='direct', status='unknown', recl=4*nCell)                     ! 4 bytes / cell (real 4)
    end if

    ! write data
    write (100,rec=d) cell(:)%z
        
    ! is leapyear ?
    if (  ( (mod(y,4) .EQ. 0) .AND. (mod(y,100) .NE. 100) )  .OR. (mod(y,400) .EQ. 0)  ) then
       leap = 1
    else
       leap = 0
    end if

    ! close file on last day
    !if ( (d .EQ. 365+leap) .OR. ( (y == yEnd) .AND. (d == ) ) ) close (100)
    if ( d .EQ. 365+leap ) then
       print*, 'saving ', trim(fileName)
       close (100)
    end if
    !
    !
  case (2)
    !---------------------------
    ! write variogram parameters
    !---------------------------
   
    ! print varfit cross-val depending on variogram estimation (true or false)
    if ( flagVario == .TRUE. ) then    !
      fileName=trim(dataPathOut)//'varFit.txt'
      inquire(201, OPENED = wasOpened)
      if (.not.wasOpened) then
         open (201, file=filename, status='unknown', action='write')
         write (201, 200) 'nugget', 'sill', 'range', 'Type', 'Easting', 'Northing','BIAS', 'RMSE', 'r'
      end if
      write (201, 201) (beta(i),i=1,nParam), vType, (xl+xr)*0.5_dp, (yd+yu)*0.5_dp, E(1), E(2), E(7)
      close(201)

      fileName = trim(dataPathOut)//trim('variogram.dat')
      open (21, file=trim(fileName), status='unknown', action='write')
      write(21,203)
      do k=1,nbins
         if (nh(k) >0)   write(21,204) (gamma(k,j), j=1,2),  tVar(gamma(k,1),beta(1),beta(2),beta(3)),  nh(k)
      end do
      close(21)
    end if 

  end select
  !
  ! formats
  110 format (i4,'.bin')
  200 format (3a12  , a6, 2a11  , 3a9)
  201 format (3e12.4, i6, 2f11.1, 3f9.4)
  203 format (19x,'h',12x,'gamma(h)',12x,'g_cal(h)', 6x, 'N(h)')
  204 format (3es20.5,i10)

  !
end subroutine WriteDataMeteo


!  *************************************************************************
!
!  SUBROUTINE   WRITE HEADER GRIDs
!
!  *************************************************************************
subroutine writeHeader(ic, fName)
  use mainVar
  use mo_kind, only         : i4
  implicit none
  integer(i4), intent(in)        :: ic                  ! input channel 
  character(256), intent(in)    :: fName
  !
  open (unit=ic, file=fName, status='unknown')
  write (ic, 1)  'ncols       ',    gridMeteo%ncols
  write (ic, 1)  'nrows       ',    gridMeteo%nrows
  write (ic, 2)  'xllcorner   ',    gridMeteo%xllcorner
  write (ic, 2)  'yllcorner   ',    gridMeteo%yllcorner
  write (ic, 1)  'cellsize    ',    gridMeteo%cellsize
  write (ic, 1)  'NODATA_value',    gridMeteo%nodata_value
  close (ic)
  ! formats
  1 format (a12, 2x, i10)
  2 format (a12, 2x, f10.1)
  !
end subroutine writeHeader


