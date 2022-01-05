!> @brief associates the name of the DCTL with the actual subroutine
!> @details .
subroutine assoc_dctl_ptr(modelname,dctl_ptr)

   use MODELING

   implicit none

   character(len=20), intent(inout):: modelname
   procedure(dctl), pointer, intent(out) :: dctl_ptr
   ! external dctl_line_prot

   select case (modelname)
       
   !    case('dctl_line_prot')
   !      dctl_ptr=>dctl_line_prot

   case default

   end select

end subroutine assoc_dctl_ptr

