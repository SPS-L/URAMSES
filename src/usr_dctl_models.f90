!> @file
!> @brief Dispatch table mapping discrete controller model names to procedure pointers
!! @details Contains the assoc_dctl_ptr subroutine, which resolves a discrete
!!          controller (DCTL) model name to the corresponding procedure pointer.
!!          Discrete controllers implement event-driven logic such as protection
!!          relays and switching devices.  The only built-in model currently
!!          registered is dctl_line_prot (line protection relay).  DLL-supplied
!!          models are checked first on Windows/Intel builds.  The prefix "dctl_"
!!          is added automatically when absent.

!> @brief Map a discrete controller model name to its procedure pointer
!! @details Prepends "dctl_" to modelname if not already present, then queries
!!          loaded DLL modules (Windows/Intel) before the built-in select-case table.
!> @param[inout] modelname  Model name string (may be modified to add "dctl_" prefix)
!> @param[out]   dctl_ptr   On exit, points to the discrete-controller subroutine for the named model
subroutine assoc_dctl_ptr(modelname,dctl_ptr)

   use MODELING

   implicit none

   character(len=20), intent(inout):: modelname    !< model name; "dctl_" prefix added if absent
   character(len=25):: modelname5                  !< internal name with guaranteed "dctl_" prefix
   procedure(dctl), pointer, intent(out) :: dctl_ptr  !< procedure pointer set to the matched discrete controller
   external dctl_line_prot
   integer(C_INTPTR_T) :: p_USERFUNC  !< address returned by GetProcAddress for DLL lookup
   integer i, ret
#if defined __INTEL_COMPILER && (defined _WIN64 || defined _WIN32)
   integer(BOOL) :: free_status
#endif

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
       
       case('dctl_line_prot')
         dctl_ptr=>dctl_line_prot

   case default

   end select

end subroutine assoc_dctl_ptr

