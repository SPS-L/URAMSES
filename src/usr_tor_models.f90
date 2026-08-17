!> @file
!> @brief Dispatch table mapping torque/governor model names to procedure pointers
!! @details Contains the assoc_torque_ptr subroutine, which resolves a governor or
!!          mechanical-torque model name to the corresponding procedure pointer.
!!          Built-in models include simplified ENTSO-E governor (ENTSOE_simp),
!!          hydro governor (HYGOV), gas turbine (GAST), steam governor (TGOV1D /
!!          TGOV1), and the DEGOV1 diesel governor.  The prefix "tor_" is added
!!          automatically when absent.

!> @brief Map a torque/governor model name to its procedure pointer
!! @details Prepends "tor_" to modelname if not already present, then resolves it
!!          against the built-in select-case table.
!> @param[inout] modelname  Model name string (may be modified to add "tor_" prefix)
!> @param[out]   tor_ptr    On exit, points to the torque subroutine for the named model
subroutine assoc_torque_ptr(modelname,tor_ptr)

   use MODELING

   implicit none

   character(len=20), intent(inout):: modelname    !< model name; "tor_" prefix added if absent
   character(len=24):: modelname4                  !< internal name with guaranteed "tor_" prefix
   procedure(torque), pointer, intent(out) :: tor_ptr  !< procedure pointer set to the matched torque model
   external :: tor_ENTSOE_simp, tor_DEGOV1
   external :: tor_hygov, tor_GAST, tor_TGOV1D
   !<<STEPSS-GUI:EXTERNALS>>

   if(modelname(1:4)=='tor_')then
      modelname4=modelname
   else
      modelname4='tor_'//modelname
   endif

   select case (modelname4)
   case('tor_ENTSOE_simp')
         tor_ptr=>tor_ENTSOE_simp

   case('tor_HYGOV')
         tor_ptr=>tor_hygov
   case('tor_GAST')
         tor_ptr=>tor_GAST
   case('tor_TGOV1')
         tor_ptr=>tor_TGOV1D
   case('tor_DEGOV1')
         tor_ptr=>tor_DEGOV1

   !<<STEPSS-GUI:CASES>>
   end select

end subroutine assoc_torque_ptr
