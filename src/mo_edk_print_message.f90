!> \file    mo_edk_print_message.f90
!> \brief   \copybrief mo_edk_print_message
!> \details \copydetails mo_edk_print_message

!> \brief   Module containing messages for the EDK program.
!> \copyright Copyright 2005-\today, the CHS Developers, Sabine Attinger: All rights reserved.
!! EDK is released under the LGPLv3+ license \license_note
module mo_edk_print_message

  implicit none

  private

  public :: print_start_message
  public :: print_end_message

contains

  !> \brief   Start up messages for the EDK program.
  subroutine print_start_message()

    use mo_edk_info, only: version
    use mo_message, only: message

    implicit none

    call message('')
    call message('==============================================')
    call message('!!                                          !!')
    call message('!!          THE KRIGING PROGRAM             !!')
    call message('!!              VERSION '// trim(version) //'          !!')
    call message('!!                                          !!')
    call message('==============================================')
    call message('')

  end subroutine print_start_message

  !> \brief   Ending messages for the EDK program.
  subroutine print_end_message

    use mo_message, only: message

    implicit none

    call message('')
    call message('==============================================')
    call message('!!                                          !!')
    call message('!!          THE KRIGING PROGRAM             !!')
    call message('!!              IS FINISHED!                !!')
    call message('!!                                          !!')
    call message('==============================================')
    call message('')

  end subroutine print_end_message

end module mo_edk_print_message
