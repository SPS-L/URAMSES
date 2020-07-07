!> \mainpage Welcome to ramses
!! \tableofcontents
!! \section intro_sec Introduction
!!
!! Please check the code for comments and the user guide.
!!
!! @file
!! @brief Main file

!> @brief Main routine executed at startup
!> @details .
program ramses_init
   use c_interface_mod, only: get_named_prm
   implicit none

   interface ramses_interface
      integer(c_int) function ramses(ccmdname, coutputfile) bind(C, name="ramses")
         use, intrinsic :: ISO_C_BINDING
         character(c_char), dimension(*), optional, intent(in) :: ccmdname, coutputfile
      end function ramses
   end interface ramses_interface

   call exit(ramses()) 

end program ramses_init
