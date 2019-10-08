!> @file
!> @brief associates the name of the DCTL with the actual subroutine

!> @brief associates the name of the DCTL with the actual subroutine
!> @details .
subroutine assoc_dctl_ptr(modelname,dctl_ptr)

   use MODELING

   implicit none

   character(len=20), intent(inout):: modelname
   character(len=25):: modelname5
   procedure(dctl), pointer, intent(out) :: dctl_ptr
   ! external dctl_OLTC2
   integer(C_INTPTR_T) :: p_USERFUNC
   integer(BOOL) :: free_status
   integer i, ret

   if(modelname(1:5)=='dctl_')then
      modelname5=modelname
   else
      modelname5='dctl_'//modelname
   endif

#if defined __INTEL_COMPILER && (defined _WIN64 || defined _WIN32)
   do i=1,dll_handleno
      if (dll_handle(i) .ne. NULL) then
         p_USERFUNC = GetProcAddress (hModule=dll_handle(i), lpProcName=trim(modelname5)//C_NULL_CHAR)
         if (p_USERFUNC .ne. NULL) then
            call C_F_PROCPOINTER (TRANSFER(p_USERFUNC, C_NULL_FUNPTR), dctl_ptr)
            return
         endif
      endif
   enddo
#endif

   select case (modelname5)

   case default

   end select

end subroutine assoc_dctl_ptr

