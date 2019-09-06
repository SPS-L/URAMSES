!> @file
!> @brief associates data to subroutines

!> @brief associates the name of the injector with the actual subroutine
!> @details .
subroutine assoc_twop_ptr(modelname,twop_ptr)

   use MODELING

   implicit none

   character(len=20), intent(in):: modelname
   procedure(twop_injector), pointer, intent(out) :: twop_ptr
   ! external twop_HQSVC

   twop_ptr=>null()
   
   select case (modelname)
   !   case('twop_HQSVC')
   !      twop_ptr => twop_HQSVC
         
   end select

end subroutine assoc_twop_ptr
