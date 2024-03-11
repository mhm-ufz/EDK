!> \file    mo_edk_types.f90
!> \brief   \copybrief mo_edk_types
!> \details \copydetails mo_edk_types

!> \brief   Module for EDK distance Type.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    21.10.2020
!> \details This file contains memory efficient types for EDK.
!> \copyright Copyright 2005-\today, the CHS Developers, Sabine Attinger: All rights reserved.
!! EDK is released under the LGPLv3+ license \license_note
module mo_edk_types
  use mo_kind, only: i4, dp
  implicit none
  !> \class dist_t
  !> \brief Cell distances on demand.
  type dist_t
    integer(i4)                                :: ncell       !< no of cells
    integer(i4)                                :: nstat       !< no of stations
    real(dp)                                   :: noDataValue !< masking value
    real(dp), dimension(:,:), allocatable      :: cell_pos    !< (x,y) pos of cells (i, j) - i: cell index; j: 1 - x, 2 - y
    real(dp), dimension(:,:), allocatable      :: stat_pos    !< (x,y) pos of stations (i, j) - i: station index; j: 1 - x, 2 - y
    real(dp), dimension(:,:), allocatable      :: stat_z      !< z values (i,j) - i: station; j - timeID
  contains
    !> \copydoc mo_edk_types::clean_dist
    procedure :: clean => clean_dist !< \see mo_edk_types::clean_dist
    !> \copydoc mo_edk_types::get_dist_cell_station
    procedure :: getCS => get_dist_cell_station !< \see mo_edk_types::get_dist_cell_station
    !> \copydoc mo_edk_types::get_dist_station_station
    procedure :: getSS => get_dist_station_station !< \see mo_edk_types::get_dist_station_station
    !> \copydoc mo_edk_types::get_squared_z_diff
    procedure :: getZ2 => get_squared_z_diff !< \see mo_edk_types::get_squared_z_diff
  end type dist_t

contains

  !> \brief   Clean up.
  subroutine clean_dist(self)
    implicit none
    class(dist_t), intent(inout) :: self
    if ( allocated(self%cell_pos)) deallocate (self%cell_pos)
    if ( allocated(self%stat_pos)) deallocate (self%stat_pos)
    if ( allocated(self%stat_z)) deallocate (self%stat_z)
  end subroutine clean_dist

  !> \brief   Get distance between cell i and station j.
  !> \return  Distance.
  real(dp) function get_dist_cell_station(self, i, j) result(distCS)
    implicit none
    class(dist_t), intent(in) :: self
    integer, intent(in) :: i, j
    distCS = dsqrt( ( self%cell_pos(i,1) - self%stat_pos(j,1) )**2 + ( self%cell_pos(i,2) - self%stat_pos(j,2) )** 2)
  end function get_dist_cell_station

  !> \brief   Get distance between station i and station j.
  !> \return  Distance.
  real(dp) function get_dist_station_station(self, i, j) result(distSS)
    implicit none
    class(dist_t), intent(in) :: self
    integer, intent(in) :: i, j
    distSS = dsqrt( ( self%stat_pos(i,1) - self%stat_pos(j,1) )**2 + ( self%stat_pos(i,2) - self%stat_pos(j,2) )** 2)
  end function get_dist_station_station

  !> \brief   Get squared elevation difference between station i and station j.
  !> \return  Elevation difference squared.
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
