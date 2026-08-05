!> @file
!> @brief Dispatch table mapping injector model names to procedure pointers
!! @details Contains the assoc_inj_ptr subroutine, which resolves a model-name
!!          string to the corresponding injector procedure pointer.  Resolution
!!          follows the same DLL-first, built-in-second strategy used by all
!!          usr_*_models dispatch files.  Built-in injector models include
!!          constant-power loads (PQ), voltage-fault injectors (VFAULT),
!!          voltage-frequency-dependent loads (vfd_load), inverter-based generators
!!          (IBG, WT3, WT4, BESS), and grid-forming/following converters (GFOR, GFOL).
!!          The prefix "inj_" is automatically prepended when absent.

!> @brief Map an injector model name to its procedure pointer
!! @details Prepends "inj_" to modelname if not already present, then queries
!!          loaded DLL modules (Windows/Intel) and the built-in select-case table.
!> @param[inout] modelname  Model name string (may be modified to add "inj_" prefix)
!> @param[out]   inj_ptr    On exit, points to the injector subroutine for the named model
subroutine assoc_inj_ptr(modelname,inj_ptr)

   use MODELING

   implicit none

   character(len=20), intent(inout):: modelname    !< model name; "inj_" prefix added if absent
   character(len=24):: modelname4                  !< internal name with guaranteed "inj_" prefix
   procedure(injector), pointer, intent(out) :: inj_ptr  !< procedure pointer set to the matched injector
   external inj_vfault ! Do not remove
   external inj_vfd_load
   external inj_PQ, inj_IBG, inj_WT3, inj_WT4, inj_BESS, inj_GFOR, inj_GFOL
   external inj_PMU
   integer(C_INTPTR_T) :: p_USERFUNC  !< address returned by GetProcAddress for DLL lookup
   integer i, ret
#if defined __INTEL_COMPILER && (defined _WIN64 || defined _WIN32)
   integer(BOOL) :: free_status
#endif

   if(modelname(1:4)=='inj_')then
      modelname4=modelname
   else
      modelname4='inj_'//modelname
   endif

#if defined __INTEL_COMPILER && (defined _WIN64 || defined _WIN32)
   do i=1,dll_handleno
      if (dll_handle(i) .ne. NULL) then
         p_USERFUNC = GetProcAddress (hModule=dll_handle(i), lpProcName=trim(modelname4)//C_NULL_CHAR)
         if (p_USERFUNC .ne. NULL) then
            call C_F_PROCPOINTER (TRANSFER(p_USERFUNC, C_NULL_FUNPTR), inj_ptr)
            return
         endif
      endif
   enddo
#endif

   select case (modelname4)

      case('inj_VFAULT') ! Do not remove
         inj_ptr=>inj_vfault

      case('inj_vfd_load')
         inj_ptr=>inj_vfd_load

      case('inj_PQ')
         inj_ptr=>inj_PQ

      case('inj_IBG')
         inj_ptr=>inj_IBG

      case('inj_WT3')
         inj_ptr=>inj_WT3

      case('inj_WT4')
         inj_ptr=>inj_WT4

      case('inj_BESS')
         inj_ptr=>inj_BESS

      case('inj_GFOL')
         inj_ptr=>inj_GFOL

      case('inj_GFOR')
         inj_ptr=>inj_GFOR

      case('inj_PMU')
         inj_ptr=>inj_PMU


   end select

end subroutine assoc_inj_ptr
