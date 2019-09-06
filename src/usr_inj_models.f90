!> @file
!> @brief associates the name of the injector with the actual subroutine

!> @brief associates the name of the injector with the actual subroutine
!> @details .
subroutine assoc_inj_ptr(modelname,inj_ptr)

   use MODELING

   implicit none

   character(len=20), intent(in):: modelname
   procedure(injector), pointer, intent(out) :: inj_ptr
   external inj_vfault


   inj_ptr=>null()
   
   select case (modelname)

      case('VFAULT')
         inj_ptr=>inj_vfault
         
      
   end select

end subroutine assoc_inj_ptr
