!> @file
!> @brief associates the name of the torque with the actual subroutine

!> @brief associates the name of the torque with the actual subroutine
!> @details .
subroutine assoc_torque_ptr(modelname,tor_ptr)

   use MODELING

   implicit none

   character(len=20), intent(in):: modelname
   procedure(torque), pointer, intent(out) :: tor_ptr
   !external tor_ENTSOE_simp

   tor_ptr => null()

   select case (modelname)

      !case('tor_ENTSOE_simp')
      !   tor_ptr => tor_ENTSOE_simp
         
   end select

end subroutine assoc_torque_ptr
