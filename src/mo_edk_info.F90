!> \file    mo_edk_info.f90
!> \copydoc mo_edk_info

#ifndef EDKVERSION
#define EDKVERSION "0.0.0-dev0"
#endif
#ifndef EDKDATE
#define EDKDATE __DATE__
#endif
#define set_version(x) character(len = *), parameter :: version = x
#define set_date(x) character(len = *), parameter :: version_date = x

!> \brief   module with edk program information
!> \details Provides information about the edk program as parameters.
!!          The \p version parameter will be set during compilation
!!          to content of \a version.txt.
!!          The \p version_date parameter will be set during compilation
!!          to content of \a version_date.txt,
!!          if it is a release version, otherwise it will be the current date.
!> \authors Sebastian Mueller
!> \date    May 2022
!> \copyright Copyright 2005-\today, the CHS Developers, Sabine Attinger: All rights reserved.
!! EDK is released under the LGPLv3+ license \license_note
module mo_edk_info

  implicit none

  set_version(EDKVERSION)
  !< Current edk version

  set_date(EDKDATE)
  !< Time of current edk version release

  !> Driver file
  character(len = *), parameter :: file_main = 'EDK_driver.f90'
  !> Namelist file name
  character(:), allocatable :: file_namelist ! = 'edk.nml'

end module mo_edk_info
