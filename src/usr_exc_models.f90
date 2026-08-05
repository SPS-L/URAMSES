!> @file
!> @brief Dispatch table mapping excitation controller model names to procedure pointers
!! @details Contains the assoc_exciter_ptr subroutine, which resolves a model-name
!!          string to the corresponding exciter procedure pointer.  Resolution order:
!!          1. On Windows/Intel builds, loaded DLL modules are queried first via
!!             GetProcAddress so that user-supplied models take precedence.
!!          2. Built-in models are resolved through a select-case block covering
!!             all excitation controllers shipped with RAMSES (e.g. ST1A, SEXS,
!!             AC1A, AC4A, EXPIC1, and generic IEEE models).
!!          The prefix "exc_" is automatically prepended to the model name if not
!!          already present.

!> @brief Map an excitation controller model name to its procedure pointer
!! @details Prepends "exc_" to modelname if not already present, then attempts
!!          to resolve the resulting string against loaded DLL modules (Windows/Intel
!!          only) and subsequently against the built-in model select-case table.
!!          The caller must pass a valid procedure pointer variable; if the model
!!          name is not found the pointer is left unchanged and no error is raised.
!> @param[inout] modelname  Model name string (may be modified to add "exc_" prefix)
!> @param[inout] exc_ptr    On exit, points to the exciter subroutine for the named model
subroutine assoc_exciter_ptr(modelname,exc_ptr)

   use MODELING

   implicit none

   character(len=20), intent(inout):: modelname    !< model name; "exc_" prefix added if absent
   character(len=24):: modelname4                  !< internal name with guaranteed "exc_" prefix
   procedure(exciter), pointer, intent(inout) :: exc_ptr  !< procedure pointer set to the matched exciter
   external :: exc_KUNDUR, exc_GENERIC3, exc_GENERIC4
   external :: exc_ENTSOE_simp, exc_ST1A, exc_AC1A, exc_AC4A, exc_IEEET5
   external :: exc_ST1A_IEEEST, exc_ST1A_PSS4B, exc_ST1A_PSS2B
   external :: exc_EXPIC1_PSS2B, exc_SEXS, exc_SEXS_IEEEST
   ! Models built from custom_models/ in this repository. Declare your own
   ! here alongside the pre-compiled ones above.
   external :: exc_ENTSOE_lim
   integer(C_INTPTR_T) :: p_USERFUNC  !< address returned by GetProcAddress for DLL lookup
   integer i, ret
#if defined __INTEL_COMPILER && (defined _WIN64 || defined _WIN32)
   integer(BOOL) :: free_status
#endif

   if(modelname(1:4)=='exc_')then
      modelname4=modelname
   else
      modelname4='exc_'//modelname
   endif


#if defined __INTEL_COMPILER && (defined _WIN64 || defined _WIN32)
   do i=1,dll_handleno
      if (dll_handle(i) .ne. NULL) then
         p_USERFUNC = GetProcAddress (hModule=dll_handle(i), lpProcName=trim(modelname4)//C_NULL_CHAR)
         if (p_USERFUNC .ne. NULL) then
            call C_F_PROCPOINTER (TRANSFER(p_USERFUNC, C_NULL_FUNPTR), exc_ptr)
            return
         endif
      endif
   enddo
#endif

   select case (modelname4)
      case('exc_kundur')
         exc_ptr=>exc_KUNDUR
      case('exc_ENTSOE_simp')
         exc_ptr=>exc_ENTSOE_simp
      case('exc_ST1A')
         exc_ptr=>exc_ST1A
      case('exc_SEXS')
         exc_ptr=>exc_SEXS
      case('exc_SEXS_IEEEST')
         exc_ptr=>exc_SEXS_IEEEST
      case('exc_GENERIC3')
         exc_ptr=>exc_GENERIC3
      case('exc_GENERIC4')
         exc_ptr=>exc_GENERIC4
      case('exc_AC1A')
         exc_ptr=>exc_AC1A
      case('exc_AC4A')
         exc_ptr=>exc_AC4A
      case('exc_IEEET5')
         exc_ptr=>exc_IEEET5
      case('exc_ST1A_IEEEST')
         exc_ptr=>exc_ST1A_IEEEST
      case('exc_ST1A_PSS4B')
         exc_ptr=>exc_ST1A_PSS4B
      case('exc_ST1A_PSS2B')
         exc_ptr=>exc_ST1A_PSS2B
      case('exc_EXPIC1_PSS2B')
         exc_ptr=>exc_EXPIC1_PSS2B

      ! Models compiled from custom_models/ in this repository. Add your own
      ! below; the name is what a .dat file's EXC record refers to, with or
      ! without the "exc_" prefix.
      case('exc_ENTSOE_lim')
         exc_ptr=>exc_ENTSOE_lim

      case default
   end select

end subroutine assoc_exciter_ptr
