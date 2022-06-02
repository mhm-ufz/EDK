module mo_edk_print_message

  implicit none

  private

  public :: print_start_message
  public :: print_end_message

contains

  subroutine print_start_message()

    use mo_edk_info, only: version
    use mo_message, only: message

    implicit none

    call message('')
    call message('==============================================')
    call message('!!                                          !!')
    call message('!!          THE KRIGING PROGRAM             !!')
    call message('!!              VERSION '// trim(version) //'                 !!')
    call message('!!                                          !!')
    call message('==============================================')
    call message('')

  end subroutine print_start_message

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
