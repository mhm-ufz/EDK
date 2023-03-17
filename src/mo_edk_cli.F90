!> \file    mo_edk_cli.f90
!> \copydoc mo_edk_cli

!> \brief   Module to parse command line arguments of edk.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Jun 2022
!> \details A simple parser for command line arguments for edk.
!!          You can pass the path to the config namelist (main.dat by default)
!!
!!          You can also pass the CWD as plain last argument and get a help or version text.
!!
!!          \code{.sh}
!!          edk -h
!!          edk -v
!!          \endcode
!> \copyright Copyright 2005-\today, the CHS Developers, Sabine Attinger: All rights reserved.
!! EDK is released under the LGPLv3+ license \license_note
module mo_edk_cli

#ifdef NAG
  USE f90_unix_dir, ONLY : CHDIR
#endif

  implicit none

  private

  public :: parse_command_line

contains

  !> \brief parse the given command line arguments.
  subroutine parse_command_line()
    use mo_cli, only: cli_parser
    use mo_edk_info, only: version, file_namelist

    implicit none

    type(cli_parser) :: parser

    parser = cli_parser( &
      description="The External Drift Kriging - EDK program", &
      add_version_option=.true., &
      version=version)

    call parser%add_option( &
      name="cwd", &
      blank=.true., &
      help="The desired working directory (optional).")

    call parser%add_option( &
      name="nml", &
      s_name="n", &
      has_value=.true., &
      value_name="path", &
      default="edk.nml", &
      help="The edk configuration namelist.")

    ! parse given command line arguments
    call parser%parse()

    ! change working directory first
    if (parser%option_was_read("cwd")) call chdir(parser%option_value("cwd"))

    ! set nml file path
    file_namelist = parser%option_value("nml")

  end subroutine parse_command_line

end module mo_edk_cli
