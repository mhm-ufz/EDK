!******************************************************************************
!  MODULE for EDK distance Type
!  AUTHOR  Sebastian Mueller (UFZ)
!  DESCRIPTION
!          This file contains memory efficient types for EDK
!  UPDATES
!          Created        21.10.2020
!******************************************************************************
module mo_edk_types
  use mo_kind, only: i4, dp
  implicit none
  ! parameters
  type dist_t
    integer(i4)                                :: ncell       ! no of cells
    integer(i4)                                :: nstat       ! no of stations
    real(dp)                                   :: noDataValue ! masking value
    ! positions: (i, j) - i: cell/station index; j: 1 - x, 2 - y
    real(dp), dimension(:,:), allocatable      :: cell_pos    ! (x,y) pos of cells
    real(dp), dimension(:,:), allocatable      :: stat_pos    ! (x,y) pos of stations
    ! z values (i,j) - i: station; j - timeID
    real(dp), dimension(:,:), allocatable      :: stat_z
  contains
    procedure :: clean => clean_dist
    procedure :: getCS => get_dist_cell_station
    procedure :: getSS => get_dist_station_station
    procedure :: getZ2 => get_squared_z_diff
  end type dist_t

contains

  subroutine clean_dist(self)
    implicit none
    class(dist_t), intent(inout) :: self
    if ( allocated(self%cell_pos)) deallocate (self%cell_pos)
    if ( allocated(self%stat_pos)) deallocate (self%stat_pos)
    if ( allocated(self%stat_z)) deallocate (self%stat_z)
  end subroutine clean_dist

  real(dp) function get_dist_cell_station(self, i, j) result(distCS)
    implicit none
    class(dist_t), intent(in) :: self
    integer, intent(in) :: i, j
    distCS = dsqrt( ( self%cell_pos(i,1) - self%stat_pos(j,1) )**2 + ( self%cell_pos(i,2) - self%stat_pos(j,2) )** 2)
  end function get_dist_cell_station

  real(dp) function get_dist_station_station(self, i, j) result(distSS)
    implicit none
    class(dist_t), intent(in) :: self
    integer, intent(in) :: i, j
    distSS = dsqrt( ( self%stat_pos(i,1) - self%stat_pos(j,1) )**2 + ( self%stat_pos(i,2) - self%stat_pos(j,2) )** 2)
  end function get_dist_station_station

  real(dp) function get_squared_z_diff(self, i, j, jd) result(distZ2)
    implicit none
    class(dist_t), intent(in) :: self
    integer, intent(in) :: i, j, jd
    if ( (self%stat_z(i, jd) /=  self%noDataValue   ) .and. &
         (self%stat_z(j, jd) /=  self%noDataValue   )        )   then
      distZ2 = ( self%stat_z(i, jd) - self%stat_z(j, jd) )**2
    else
      distZ2 = self%noDataValue
    end if
  end function get_squared_z_diff

end module mo_edk_types
