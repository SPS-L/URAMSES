!> @file
!> @brief Dispatch table mapping discrete controller model names to procedure pointers
!! @details Contains the assoc_dctl_ptr subroutine, which resolves a discrete
!!          controller (DCTL) model name to the corresponding procedure pointer.
!!          Discrete controllers implement event-driven logic such as protection
!!          relays and switching devices.  The only built-in model currently
!!          registered is dctl_line_prot (line protection relay).  The prefix
!!          "dctl_" is added automatically when absent.

!> @brief Map a discrete controller model name to its procedure pointer
!! @details Prepends "dctl_" to modelname if not already present, then resolves it
!!          against the built-in select-case table.
!> @param[inout] modelname  Model name string (may be modified to add "dctl_" prefix)
!> @param[out]   dctl_ptr   On exit, points to the discrete-controller subroutine for the named model
subroutine assoc_dctl_ptr(modelname,dctl_ptr)

   use MODELING

   implicit none

   character(len=20), intent(inout):: modelname    !< model name; "dctl_" prefix added if absent
   character(len=25):: modelname5                  !< internal name with guaranteed "dctl_" prefix
   procedure(dctl), pointer, intent(out) :: dctl_ptr  !< procedure pointer set to the matched discrete controller
   external dctl_line_prot, dctl_injprot
   !<<STEPSS-GUI:EXTERNALS>>

   if(modelname(1:5)=='dctl_')then
      modelname5=modelname
   else
      modelname5='dctl_'//modelname
   endif

   select case (modelname5)

      case('dctl_injprot')
         dctl_ptr=>dctl_injprot

      case('dctl_line_prot')
         dctl_ptr=>dctl_line_prot

   !<<STEPSS-GUI:CASES>>
   end select

end subroutine assoc_dctl_ptr

