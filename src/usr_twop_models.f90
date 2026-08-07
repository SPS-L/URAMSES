!> @file
!> @brief Dispatch table mapping two-port model names to procedure pointers
!! @details Contains the assoc_twop_ptr subroutine, which resolves a two-port
!!          device model name to the corresponding procedure pointer.  Two-port
!!          elements connect two AC buses and include HVDC links (LCC and VSC
!!          variants) and DC-line / weak-coupling-link (DCL_WCL) models.  DLL
!!          models are checked first on Windows/Intel builds; built-in models are
!!          resolved through a select-case block.  The prefix "twop_" is added
!!          automatically when absent.

!> @brief Map a two-port model name to its procedure pointer
!! @details Prepends "twop_" to modelname if not already present, then queries
!!          loaded DLL modules (Windows/Intel) before the built-in select-case table.
!> @param[inout] modelname  Model name string (may be modified to add "twop_" prefix)
!> @param[out]   twop_ptr   On exit, points to the two-port subroutine for the named model
subroutine assoc_twop_ptr(modelname,twop_ptr)

   use MODELING

   implicit none

   character(len=20), intent(inout):: modelname    !< model name; "twop_" prefix added if absent
   character(len=25):: modelname5                  !< internal name with guaranteed "twop_" prefix
   procedure(twop_injector), pointer, intent(out) :: twop_ptr  !< procedure pointer set to the matched two-port model
   external twop_HVDC_LCC,  twop_HVDC_VSC, twop_HVDC_VSC_SC, twop_DCL_WCL
   !<<STEPSS-GUI:EXTERNALS>>
   integer(C_INTPTR_T) :: p_USERFUNC  !< address returned by GetProcAddress for DLL lookup
   integer i, ret
#if defined __INTEL_COMPILER && (defined _WIN64 || defined _WIN32)
   integer(BOOL) :: free_status
#endif

   if(modelname(1:5)=='twop_')then
      modelname5=modelname
   else
      modelname5='twop_'//modelname
   endif

#if defined __INTEL_COMPILER && (defined _WIN64 || defined _WIN32)
   do i=1,dll_handleno
      if (dll_handle(i) .ne. NULL) then
         p_USERFUNC = GetProcAddress (hModule=dll_handle(i), lpProcName=trim(modelname5)//C_NULL_CHAR)
         if (p_USERFUNC .ne. NULL) then
            call C_F_PROCPOINTER (TRANSFER(p_USERFUNC, C_NULL_FUNPTR), twop_ptr)
            return
         endif
      endif
   enddo
#endif

   select case (modelname5)

   case('twop_HVDC_LCC')
         twop_ptr=>twop_HVDC_LCC

   case('twop_HVDC_VSC')
         twop_ptr => twop_HVDC_VSC

   case('twop_DCL_WCL')
         twop_ptr => twop_DCL_WCL

   case('twop_HVDC_VSC_SC')
         twop_ptr => twop_HVDC_VSC_SC

   !<<STEPSS-GUI:CASES>>
   end select

end subroutine assoc_twop_ptr
