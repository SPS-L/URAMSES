subroutine assoc_exciter_ptr(modelname,exc_ptr)

   use MODELING

   implicit none

   character(len=20), intent(in):: modelname
   procedure(exciter), pointer, intent(out) :: exc_ptr
   ! external exc_ENTSOE_simp


   select case (modelname)

   !   case('exc_ENTSOE_simp')
   !      exc_ptr => exc_ENTSOE_simp


   end select

end subroutine assoc_exciter_ptr
