!> @file
!> @brief associates the name of the torque with the actual subroutine

!> @brief associates the name of the torque with the actual subroutine
!> @details .
subroutine assoc_torque_ptr(modelname,tor_ptr)

   use MODELING

   implicit none

   character(len=20), intent(in):: modelname
   procedure(torque), pointer, intent(out) :: tor_ptr
   !external tor_DEGOV1

   select case (modelname)


      !case('DEGOV1')
      !   tor_ptr => tor_DEGOV1
         
   end select

end subroutine assoc_torque_ptr
