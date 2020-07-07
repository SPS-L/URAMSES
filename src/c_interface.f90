!> @file
!> @brief List of subroutines for C interfacing

module c_interface_mod

   use, intrinsic :: iso_c_binding
   use sets
   use DIMENSIONS, only: mxzon
   use readline_utility
   implicit none

   !> This value is computed once to give the maximum number of parameters used in any model
   integer(c_int) :: mxprm = 0

   !> This is a vector with pointers to subsystems
   type(setType), dimension(mxzon) :: subsystem


   interface ramses_interface
      integer(c_int) function ramses(ccmdname) bind(C, name="ramses")
         use, intrinsic :: ISO_C_BINDING
         character(c_char), dimension(*), optional, intent(in) :: ccmdname
      end function ramses
   end interface ramses_interface

contains


   
!> @brief This is a C wrapper for getting voltage magnitude
!> @details .
integer(c_int) function get_volt_mag(bus_name, bus_volt) bind(C, name="get_volt_mag")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_volt_mag
#endif
use VOLTAGE, only: vx_h, vy_h
use search_mod, only: searn
implicit none
doubleprecision :: vx, vy
character(c_char), dimension(*), intent(in) :: bus_name
real(c_double), intent(inout) :: bus_volt
character(len=20) :: fname
integer :: busn_num

fname=c_to_f_string(bus_name)
!write(log,"(a)") fname
call searn(fname(1:18),busn_num)
if (busn_num==0) then
   bus_volt=0.d0
   get_volt_mag=1
   return
endif

!!$omp atomic read
vx=vx_h(busn_num)

!!$omp atomic read
vy=vy_h(busn_num)

bus_volt=hypot(vx, vy)
get_volt_mag=0
end function get_volt_mag

!> @brief This is a C wrapper for getting voltage phase
!> @details .
integer(c_int) function get_volt_pha(bus_name, bus_pha) bind(C, name="get_volt_pha")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_volt_pha
#endif
use VOLTAGE, only: vx_h, vy_h
use search_mod, only: searn
use SETTINGS, only: rad
implicit none
doubleprecision :: vx, vy
character(c_char), dimension(*), intent(in) :: bus_name
real(c_double), intent(inout) :: bus_pha
character(len=20) :: fname
integer :: busn_num

fname=c_to_f_string(bus_name)
!write(log,"(a)") fname
call searn(fname(1:18),busn_num)
if (busn_num==0) then
   bus_pha=0.d0
   get_volt_pha=1
   return
endif

!!$omp atomic read
vx=vx_h(busn_num)

!!$omp atomic read
vy=vy_h(busn_num)

bus_pha=datan2(vy,vx)
get_volt_pha=0
end function get_volt_pha

!> @brief This is a C wrapper for getting line power
!> @details .
integer(c_int) function get_line_pow(line_name, p_orig, q_orig, p_extr, q_extr) bind(C, name="get_line_pow")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_line_pow
#endif
use observ_mod, only: pqbra
use search_mod, only: searb
implicit none
character(c_char), dimension(*), intent(in) :: line_name
real(c_double), intent(inout) :: p_orig, q_orig, p_extr, q_extr
double precision :: p_orig_t, q_orig_t, p_extr_t, q_extr_t
character(len=20) :: fname
integer :: bra_num

bra_num = 0
fname=c_to_f_string(line_name)
!write(log,"(a)") fname
call searb(fname(1:20),bra_num)
!write(log,"(i)") bra_num
if (bra_num==0) then
   get_line_pow=1
   return
endif

call pqbra(bra_num, p_orig_t, q_orig_t, p_extr_t, q_extr_t, .true.)

p_orig = p_orig_t
q_orig = q_orig_t
p_extr = p_extr_t
q_extr = q_extr_t

get_line_pow=0
end function get_line_pow

!> @brief This is a C wrapper for getting line current
!> @details .
integer(c_int) function get_line_cur(line_name, ix_orig, iy_orig, ix_extr, iy_extr) bind(C, name="get_line_cur")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_line_cur
#endif
use observ_mod, only: pqbra
use search_mod, only: searb
implicit none
character(c_char), dimension(*), intent(in) :: line_name
real(c_double), intent(inout) :: ix_orig, iy_orig, ix_extr, iy_extr
double precision :: ix_orig_t, iy_orig_t, ix_extr_t, iy_extr_t
character(len=20) :: fname
integer :: bra_num

bra_num = 0
fname=c_to_f_string(line_name)
!write(log,"(a)") fname
call searb(fname(1:20),bra_num)
!write(log,"(i)") bra_num
if (bra_num==0) then
   get_line_cur=1
   return
endif

call pqbra(bra_num, ix_orig_t, iy_orig_t, ix_extr_t, iy_extr_t, .false.)

ix_orig = ix_orig_t
iy_orig = iy_orig_t
ix_extr = ix_extr_t
iy_extr = iy_extr_t

get_line_cur=0
end function get_line_cur


!> @brief Get Jacobian matrix
!> @details .
integer(c_int) function get_Jac() bind(C, name="get_Jac")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_Jac
#endif
use DISTURB, only: jacfile
use simul_decomp_mod, only: dumpjac, forcejac
use SimTIME, only: pause_time
use SETTINGS,only: error_flag
implicit none
jacfile="py"
dumpjac=.true.
forcejac=.true.
pause_time = pause_time+0.001d0
error_flag=.false.
get_Jac=ramses()
end function get_Jac

!> @brief Get current simulation time
!> @details .
real(c_double) function get_sim_time() bind(C, name="get_sim_time")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_sim_time
#endif
use SimTIME, only: t_h
implicit none
!!$omp atomic read
get_sim_time=t_h(0)
end function get_sim_time

!> @brief Get largest float in Fortran
!> @details .
real(c_double) function get_huge_double() bind(C, name="get_huge_double")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_huge_double
#endif
implicit none
get_huge_double=huge(1.d0)
end function get_huge_double

!> @brief Returns the value of a named observable
!> @details .
integer(c_int) function get_named_obs(comp_type, comp_name, obs_name, obs_value) bind(C, name="get_named_obs")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_named_obs
#endif
use simul_decomp_mod
use search_mod
use observ_mod, only: pqsync
use SYNC, only: adxsync, xsync_h
implicit none
character(c_char), dimension(*), intent(in) :: comp_name, comp_type, obs_name
real(c_double), intent(inout) :: obs_value
character(len=20) :: fcomp_name
character(len=10) :: fobs_name
character(len=10) :: fcomp_type ! 'EXC','TOR','INJ','TWOP','DCTL','SYN'
character(len=200) :: msg
integer :: j, i, nbobs, k, ni, nj, shift
character(len=10) :: names(mxobsuser)
double precision :: tmpObsBuf(mxobsuser), p, q, if, v

fcomp_name=c_to_f_string(comp_name)
fobs_name=c_to_f_string(obs_name)
fcomp_type=c_to_f_string(comp_type)

select case(fcomp_type)
   case ('SYN')
      call searm(fcomp_name,j)
      if(j == 0)then
         write(msg,"('t = ',f9.4,' unknown sync machine ',a)")t_h(0),trim(fcomp_name)
         call write_msg('get_named_obs',trim(msg))
         get_named_obs=1
         return
      endif
      call pqsync(j,p,q)
      select case(fobs_name)
      case ('P')
         obs_value=p*sbases(bus_inj(j))
         get_named_obs=0
         return
      case ('Q')
         obs_value=q*sbases(bus_inj(j))
         get_named_obs=0
         return
      case ('Omega')
         obs_value=xsync_h(adxsync(j)+9,0)
         get_named_obs=0
         return
      case ('S') !SYNC MVA
         obs_value = hypot(p, q)*sbases(bus_inj(j))
         get_named_obs=0
         return
      case ('SNOM')
         obs_value = snom_sync(j)*sbases(bus_inj(j))
         get_named_obs=0
         return
      case ('PNOM')
         obs_value = pnom_sync(j)*sbases(bus_inj(j))
         get_named_obs=0
         return
      end select
      write(msg,"('t = ',f9.4,' unknown observable name ',a)")t_h(0),trim(fobs_name)
      call write_msg('get_named_obs',trim(msg))
      get_named_obs=1
      return

   case ('EXC')
      call searm(fcomp_name,j)
      if(j == 0)then
         write(msg,"('t = ',f9.4,' unknown sync machine ',a)")t_h(0),trim(fcomp_name)
         call write_msg('get_named_obs',trim(msg))
         get_named_obs=1
         return
      endif
      call def_obs_exc_model(j,exc_model(j),nbobs,names)
      ObsNameSearchLoopExc: &
      do k=1,nbobs
         if(names(k)==fobs_name)then
            tmpObsBuf=0.d0
            if(injbr(j)==1)then
               i=bus_inj(j)
               call pqsync(j,p,q)
               shift=adxsync(j)
               v=dsqrt(vx_h(i)**2+vy_h(i)**2)
               if=(xsync_h(shift+4,0)-xsync_h(shift+2,0))/llf(i) !9) SFC field current (pu machine base)
               if=if*rf(j)/puf(j)
               call eval_obs_exc_model(j,exc_model(j),syncname(j),t_h(0),v,p/snom_sync(j),q/snom_sync(j), &
                   xsync_h(shift+9,0),if,prmexc(adprmexc(j)),xsync_h(shift+nbparkeq,0),zexc(adzexc(j)),tmpObsBuf)
            else
               exit ObsNameSearchLoopExc
            endif
            obs_value=tmpObsBuf(k)
            get_named_obs=0
            return
         endif
      enddo ObsNameSearchLoopExc
      write(msg,"('t = ',f9.4,' unknown observable name ',a)")t_h(0),trim(fobs_name)
      call write_msg('get_named_obs',trim(msg))
      get_named_obs=1
      return

   case ('TOR')
      call searm(fcomp_name,j)
      if(j == 0)then
         write(msg,"('t = ',f9.4,' unknown sync machine ',a)")t_h(0),trim(fcomp_name)
         call write_msg('get_named_obs',trim(msg))
         get_named_obs=1
         return
      endif
      call def_obs_tor_model(j,tor_model(j),nbobs,names)
      ObsNameSearchLoopTor: &
      do k=1,nbobs
         if(names(k)==fobs_name)then
            tmpObsBuf=0.d0
            if(injbr(j)==1)then
               i=bus_inj(j)
               call pqsync(j,p,q)
               shift=adxsync(j)
               call eval_obs_tor_model(j,tor_model(j),syncname(j),t_h(0),p/pnom_sync(j),xsync_h(shift+9,0), &
                   prmtor(adprmtor(j)),xsync_h(shift+nbparkeq+nbxexc(j),0),ztor(adztor(j)),tmpObsBuf)
            else
               exit ObsNameSearchLoopTor
            endif
            obs_value=tmpObsBuf(k)
            get_named_obs=0
            return
         endif
      enddo ObsNameSearchLoopTor
      write(msg,"('t = ',f9.4,' unknown observable name ',a)")t_h(0),trim(fobs_name)
      call write_msg('get_named_obs',trim(msg))
      get_named_obs=1
      return

   case ('INJ')
      call seari(fcomp_name,j)
      if(j == 0)then
         write(msg,"('t = ',f9.4,' unknown injector ',a)")t_h(0),trim(fcomp_name)
         call write_msg('get_named_obs',trim(msg))
         get_named_obs=1
         return
      endif
      call def_obs_inj_model(j,inj_model(j),nbobs,names)
      ObsNameSearchLoopInj: &
      do k=1,nbobs
         if(names(k)==fobs_name)then
            tmpObsBuf=0.d0
            if(injbr(nbsync+j)==1)then
               i=bus_inj(nbsync+j)
               call eval_obs_inj_model(j,inj_model(j),injname(j),t_h(0),vx_h(i),vy_h(i),omegacoi(isl(i),0), &
                  prminj(adprminj(j)),xinj_h(adxinj(j),0),zinj(adzinj(j)),tmpObsBuf,sbases(bussubnet(j)))
            else
               exit ObsNameSearchLoopInj
            endif
            obs_value=tmpObsBuf(k)
            get_named_obs=0
            return
         endif
      enddo ObsNameSearchLoopInj
      write(msg,"('t = ',f9.4,' unknown observable name ',a)")t_h(0),trim(fobs_name)
      call write_msg('get_named_obs',trim(msg))
      get_named_obs=1
      return

   case ('TWOP')
      call seart(fcomp_name,j)
      if(j == 0)then
         write(msg,"('t = ',f9.4,' unknown twoport ',a)")t_h(0),trim(fcomp_name)
         call write_msg('get_named_obs',trim(msg))
         get_named_obs=1
         return
      endif
      call def_obs_twop_model(j,twop_model(j),nbobs,names)
      ObsNameSearchLoopTwop: &
      do k=1,nbobs
         if(names(k)==fobs_name)then
            tmpObsBuf=0.d0
            if(twopbr(j)==1)then
               ni=twop_orig(j)
               nj=twop_extr(j)
               call eval_obs_twop_model(j,twop_model(j),twopname(j),t_h(0),vx_h(ni),vy_h(ni),vx_h(nj),vy_h(nj), &
                  omegacoi(isl(ni),0),omegacoi(isl(nj),0),prmtwop(adprmtwop(j)),xtwop_h(adxtwop(j),0),ztwop(adztwop(j)), &
                  tmpObsBuf,sbases(bussubnet(ni)),sbases(bussubnet(nj)))
            else
               exit ObsNameSearchLoopTwop
            endif
            obs_value=tmpObsBuf(k)
            get_named_obs=0
            return
         endif
      enddo ObsNameSearchLoopTwop
      write(msg,"('t = ',f9.4,' unknown observable name ',a)")t_h(0),trim(fobs_name)
      call write_msg('get_named_obs',trim(msg))
      get_named_obs=1
      return

   case ('DCTL')
      call seard(fcomp_name,j)
      if(j == 0)then
         write(msg,"('t = ',f9.4,' unknown DCTL ',a)")t_h(0),trim(fcomp_name)
         call write_msg('get_named_obs',trim(msg))
         get_named_obs=1
         return
      endif
      call def_obs_dctl_model(j,dctl_model(j),nbobs,names)
      ObsNameSearchLoopDCTL: &
      do k=1,nbobs
         if(names(k)==fobs_name)then
            tmpObsBuf=0.d0
            call eval_obs_dctl_model(j,dctl_model(j),t_h(0),wdctl(adwdctl(j)),tmpObsBuf)
            obs_value=tmpObsBuf(k)
            get_named_obs=0
            return
         endif
      enddo ObsNameSearchLoopDCTL
      write(msg,"('t = ',f9.4,' unknown observable name ',a)")t_h(0),trim(fobs_name)
      call write_msg('get_named_obs',trim(msg))
      get_named_obs=1
      return
end select

get_named_obs=1

end function get_named_obs

!===================================================================================================================
integer function get_named_obs_inj(j, fobs_name, obs_value, flg_nom, flg_mva)
   use simul_decomp_mod
   implicit none
   integer, intent(in) :: j
   double precision, intent(inout) :: obs_value
   character(len=10), intent(in) :: fobs_name
   integer :: i, k, nbobs
   character(len=10) :: names(mxobsuser)
   double precision :: tmpObsBuf(mxobsuser)
   logical :: flg_nom, flg_mva
   character(len=200) :: msg

   call def_obs_inj_model(j,inj_model(j),nbobs,names)
   ObsNameSearchLoopInj: &
   do k=1,nbobs
      if(names(k)==fobs_name)then
         tmpObsBuf=0.d0
         if(injbr(nbsync+j)==1)then
            i=bus_inj(nbsync+j)
            call eval_obs_inj_model(j,inj_model(j),injname(j),t_h(0),vx_h(i),vy_h(i),omegacoi(isl(i),0), &
               prminj(adprminj(j)),xinj_h(adxinj(j),0),zinj(adzinj(j)),tmpObsBuf,sbases(bussubnet(j)))
         else
            exit ObsNameSearchLoopInj
         endif
         obs_value=tmpObsBuf(k)
         if (flg_nom) obs_value=obs_value/(vx_h(i)**2 + vy_h(i)**2)
         if (flg_mva) obs_value=obs_value*sbases(bussubnet(j))
         get_named_obs_inj=0
         return
      endif
   enddo ObsNameSearchLoopInj
   write(msg,"('t = ',f9.4,' unknown observable name ',a)")t_h(0),trim(fobs_name)
   call write_msg('get_named_obs',trim(msg))
   get_named_obs_inj=1
end function
!===================================================================================================================

!> @brief Returns the value of a named parameter
!> @details .
integer(c_int) function get_named_prm(comp_type, comp_name, prm_name, prm_value) bind(C, name="get_named_prm")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_named_prm
#endif
use search_mod
use SETTINGS, only: write_msg
use SYNC, only: nameprmexc, nameprmtor, prmexc, prmtor, adprmexc,adprmtor
use UDIM, only: nameprminj, prminj, adprminj
use TWOP, only: nameprmtwop, prmtwop, adprmtwop
use DCTL, only: adwdctl, wdctl, namewdctl
use SimTIME, only: t_h
implicit none
character(c_char), dimension(*), intent(in) :: comp_name, comp_type, prm_name
real(c_double), intent(inout) :: prm_value
character(len=20) :: fcomp_name
character(len=10) :: fprm_name
character(len=10) :: fcomp_type ! 'EXC','TOR','INJ','TWOP','DCTL'
character(len=200) :: msg
integer :: j, i

fcomp_name=c_to_f_string(comp_name)
fprm_name=c_to_f_string(prm_name)
fcomp_type=c_to_f_string(comp_type)

select case(fcomp_type)
   case ('EXC')
      call searm(fcomp_name,j)
      if(j == 0)then
         write(msg,"('t = ',f9.4,'  unknown sync machine ',a)")t_h(0),trim(fcomp_name)
         call write_msg('get_named_prm',trim(msg))
         get_named_prm=1
         return
      endif
      do i=adprmexc(j),adprmexc(j+1)-1
         if(nameprmexc(i)==fprm_name)then
            !!$omp atomic read
            prm_value=prmexc(i)
            get_named_prm=0
            return
         endif
      enddo

   case ('TOR')
      call searm(fcomp_name,j)
      if(j == 0)then
         write(msg,"('t = ',f9.4,'  unknown sync machine ',a)")t_h(0),trim(fcomp_name)
         call write_msg('get_named_prm',trim(msg))
         get_named_prm=1
         return
      endif
      do i=adprmtor(j),adprmtor(j+1)-1
         if(nameprmtor(i)==fprm_name)then
            !!$omp atomic read
            prm_value=prmtor(i)
            get_named_prm=0
            return
         endif
      enddo

   case ('INJ')
      call seari(fcomp_name,j)
      if(j == 0)then
         write(msg,"('t = ',f9.4,'  unknown injector ',a)")t_h(0),trim(fcomp_name)
         call write_msg('get_named_prm',trim(msg))
         get_named_prm=1
         return
      endif
      do i=adprminj(j),adprminj(j+1)-1
         if(nameprminj(i)==fprm_name)then
            !!$omp atomic read
            prm_value=prminj(i)
            get_named_prm=0
            return
         endif
      enddo

   case ('TWOP')
      call seart(fcomp_name,j)
      if(j == 0)then
         write(msg,"('t = ',f9.4,'  unknown twoport ',a)")t_h(0),trim(fcomp_name)
         call write_msg('get_named_prm',trim(msg))
         get_named_prm=1
         return
      endif
      do i=adprmtwop(j),adprmtwop(j+1)-1
         if(nameprmtwop(i)==fprm_name)then
            !!$omp atomic read
            prm_value=prmtwop(i)
            get_named_prm=0
            return
         endif
      enddo

   case ('DCTL')
      call seard(fcomp_name,j)
      if(j == 0)then
         write(msg,"('t = ',f9.4,'  unknown DCTL ',a)")t_h(0),trim(fcomp_name)
         call write_msg('get_named_prm',trim(msg))
         get_named_prm=1
         return
      endif
      do i=adwdctl(j),adwdctl(j+1)-1
         if(namewdctl(i)==fprm_name)then
            !!$omp atomic read
            prm_value=wdctl(i)
            get_named_prm=0
            return
         endif
      enddo
end select

get_named_prm=1

end function get_named_prm

!> @brief Return nbbus
!> @details .
integer(c_int) function get_nbbus() bind(C, name="get_nbbus")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_nbbus
#endif
use BUS, only: nbbus
implicit none
get_nbbus=nbbus
end function get_nbbus

!-----------------------------------------------------------------------------------------------------------------
!> @brief Return bus name
!> @details .
integer(c_int) function get_last_err_log(cerr_msg) bind(C, name="get_last_err_log")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_last_err_log
#endif
use SETTINGS, only: lasterrormsglog
implicit none
character(c_char), dimension(1024), intent(inout) :: cerr_msg
character(len=1024) :: ferr_msg
integer :: n, j

get_last_err_log = 0

cerr_msg(1)=c_null_char
ferr_msg=trim(lasterrormsglog)
n = len_trim(ferr_msg)
if(n>=0)then
   do j = 1, n
      cerr_msg(j) = ferr_msg(j:j)
   end do
   cerr_msg(n + 1) = c_null_char
   return
endif
get_last_err_log = 1
end function get_last_err_log

!> @brief Return bus name
!> @details .
integer(c_int) function get_bus_name(i,name) bind(C, name="get_bus_name")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_bus_name
#endif
use BUS, only: busname, nbbus
implicit none
integer(c_int), INTENT(IN), value :: i
character(c_char), dimension(21), intent(inout) :: name
character(len=21) :: f_name
integer :: n, j
name(1)=c_null_char
if(i<=nbbus .and. i>0)then
   f_name=trim(busname(i))
   n = len_trim(f_name)
   do j = 1, n
      name(j) = f_name(j:j)
   end do
   name(n + 1) = c_null_char
   get_bus_name = 0
   return
endif
get_bus_name = 1
end function get_bus_name

!> @brief Return nbsync
!> @details .
integer(c_int) function get_nbsync() bind(C, name="get_nbsync")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_nbsync
#endif
use SYNC, only: nbsync
implicit none
get_nbsync=nbsync
end function get_nbsync

!> @brief Return sync name
!> @details .
integer(c_int) function get_sync_name(i,name) bind(C, name="get_sync_name")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_sync_name
#endif
use SYNC, only: syncname, nbsync
implicit none
integer(c_int), INTENT(IN), value :: i
character(c_char), dimension(21), intent(inout) :: name
character(len=21) :: f_name
integer :: n, j
name(1)=c_null_char
if(i<=nbsync .and. i>0)then
   f_name=trim(syncname(i))
   n = len_trim(f_name)
   do j = 1, n
      name(j) = f_name(j:j)
   end do
   name(n + 1) = c_null_char
   get_sync_name = 0
   return
endif
get_sync_name = 1
end function get_sync_name


!> @brief Return nbinj
!> @details .
integer(c_int) function get_nbinj() bind(C, name="get_nbinj")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_nbinj
#endif
use UDIM, only: nbinj
implicit none
get_nbinj=nbinj
end function get_nbinj

!> @brief Return inj name
!> @details .
integer(c_int) function get_inj_name(i,name) bind(C, name="get_inj_name")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_inj_name
#endif
use UDIM, only: injname, nbinj
implicit none
integer(c_int), INTENT(IN), value :: i
character(c_char), dimension(21), intent(inout) :: name
character(len=21) :: f_name
integer :: n, j
name(1)=c_null_char
if(i<=nbinj .and. i>0)then
   f_name=trim(injname(i))
   n = len_trim(f_name)
   do j = 1, n
      name(j) = f_name(j:j)
   end do
   name(n + 1) = c_null_char
   get_inj_name = 0
   return
endif
get_inj_name = 1
end function get_inj_name

!> @brief Return nbdctl
!> @details .
integer(c_int) function get_nbdctl() bind(C, name="get_nbdctl")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_nbdctl
#endif
use DCTL, only: nbdctl
implicit none
get_nbdctl=nbdctl
end function get_nbdctl

!> @brief Return dctl name
!> @details .
integer(c_int) function get_dctl_name(i,name) bind(C, name="get_dctl_name")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_dctl_name
#endif
use DCTL, only: dctlname, nbdctl
implicit none
integer(c_int), INTENT(IN), value :: i
character(c_char), dimension(21), intent(inout) :: name
character(len=21) :: f_name
integer :: n, j
name(1)=c_null_char
if(i<=nbdctl .and. i>0)then
   f_name=trim(dctlname(i))
   n = len_trim(f_name)
   do j = 1, n
      name(j) = f_name(j:j)
   end do
   name(n + 1) = c_null_char
   get_dctl_name = 0
   return
endif
get_dctl_name = 1
end function get_dctl_name

!> @brief Return nbbra
!> @details .
integer(c_int) function get_nbbra() bind(C, name="get_nbbra")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_nbbra
#endif
use BRANCH, only: nbbra
implicit none
get_nbbra=nbbra
end function get_nbbra

!> @brief Return branch name
!> @details .
integer(c_int) function get_branch_name(i,name) bind(C, name="get_branch_name")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_branch_name
#endif
use BRANCH, only: braname, nbbra
implicit none
integer(c_int), INTENT(IN), value :: i
character(c_char), dimension(21), intent(inout) :: name
character(len=21) :: f_name
integer :: n, j
name(1)=c_null_char
if(i<=nbbra .and. i>0)then
   f_name=trim(braname(i))
   n = len_trim(f_name)
   do j = 1, n
      name(j) = f_name(j:j)
   end do
   name(n + 1) = c_null_char
   get_branch_name = 0
   return
endif
get_branch_name = 1
end function get_branch_name

!> @brief Return nbtwop
!> @details .
integer(c_int) function get_nbtwop() bind(C, name="get_nbtwop")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_nbtwop
#endif
use TWOP, only: nbtwop
implicit none
get_nbtwop=nbtwop
end function get_nbtwop

!> @brief Return twoport name
!> @details .
integer(c_int) function get_twop_name(i,name) bind(C, name="get_twop_name")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_twop_name
#endif
use TWOP, only: nbtwop, twopname
implicit none
integer(c_int), INTENT(IN), value :: i
character(c_char), dimension(21), intent(inout) :: name
character(len=21) :: f_name
integer :: n, j
name(1)=c_null_char
if(i<=nbtwop .and. i>0)then
   f_name=trim(twopname(i))
   n = len_trim(f_name)
   do j = 1, n
      name(j) = f_name(j:j)
   end do
   name(n + 1) = c_null_char
   get_twop_name = 0
   return
endif
get_twop_name = 1
end function get_twop_name

!> @brief Return nbshunt
!> @details .
integer(c_int) function get_nbshunt() bind(C, name="get_nbshunt")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_nbshunt
#endif
use SHUNT, only: nbshunt
implicit none
get_nbshunt=nbshunt
end function get_nbshunt

!> @brief Return mxprm
!> @details Gets the maximum number of parameters in the any type of model
integer(c_int) function get_mxprm() bind(C, name="get_mxprm")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_mxprm
#endif
use SYNC, only: nbsync, adprmexc, adprmtor
use UDIM, only: nbinj, adprminj
use TWOP, only: nbtwop, adprmtwop
implicit none
integer :: i

if (mxprm == 0) then
   do i=1,nbsync
      mxprm=max(mxprm, adprmexc(i+1)-adprmexc(i), adprmtor(i+1)-adprmtor(i))
   enddo
   do i=1,nbinj
      mxprm=max(mxprm, adprminj(i+1)-adprminj(i))
   enddo
   do i=1,nbtwop
      mxprm=max(mxprm, adprmtwop(i+1)-adprmtwop(i))
   enddo
endif

get_mxprm=mxprm
end function get_mxprm

!> @brief Return list of parameter names of component
!> @details .
integer(c_int) function get_comp_prm_names(comp_type,comp_name,mxprm,names) bind(C, name="get_comp_prm_names")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_comp_prm_names
#endif
use search_mod
use SETTINGS, only: write_msg
use SYNC, only: nameprmexc, nameprmtor, prmexc, prmtor, adprmexc,adprmtor
use UDIM, only: nameprminj, prminj, adprminj
use TWOP, only: nameprmtwop, prmtwop, adprmtwop
use DCTL, only: adwdctl, wdctl, namewdctl
use SimTIME, only: t_h
implicit none

integer(c_int), INTENT(in), value :: mxprm
character(c_char), dimension(11*mxprm), intent(inout) :: names
character(c_char), dimension(*), intent(inout) :: comp_name, comp_type
character(len=10) :: f_name

character(len=20) :: fcomp_name
character(len=10) :: fprm_name
character(len=10) :: fcomp_type ! 'EXC','TOR','INJ','TWOP','DCTL'
character(len=200) :: msg
integer :: j, i, k, n

fcomp_name=c_to_f_string(comp_name)
fcomp_type=c_to_f_string(comp_type)

names(1)=c_null_char
k=0

select case(fcomp_type)
   case ('EXC')
      call searm(fcomp_name,j)
      if(j == 0)then
         write(msg,"('t = ',f9.4,'  unknown sync machine ',a)")t_h(0),trim(fcomp_name)
         call write_msg('get_comp_prm_names',trim(msg))
         get_comp_prm_names=1
         return
      endif
      do i=adprmexc(j),adprmexc(j+1)-1
         if(nameprmexc(i).ne.'')then
            k=k+1
            f_name=nameprmexc(i)
            do n=1,10
               names(11*(k-1)+n)=f_name(n:n)
            enddo
            names(11*k)=' '
         endif
      enddo

   case ('TOR')
      call searm(fcomp_name,j)
      if(j == 0)then
         write(msg,"('t = ',f9.4,'  unknown sync machine ',a)")t_h(0),trim(fcomp_name)
         call write_msg('get_comp_prm_names',trim(msg))
         get_comp_prm_names=1
         return
      endif
      do i=adprmtor(j),adprmtor(j+1)-1
         if(nameprmtor(i).ne.'')then
            k=k+1
            f_name=nameprmtor(i)
            do n=1,10
               names(11*(k-1)+n)=f_name(n:n)
            enddo
            names(11*k)=' '
         endif
      enddo

   case ('INJ')
      call seari(fcomp_name,j)
      if(j == 0)then
         write(msg,"('t = ',f9.4,'  unknown injector ',a)")t_h(0),trim(fcomp_name)
         call write_msg('get_comp_prm_names',trim(msg))
         get_comp_prm_names=1
         return
      endif
      do i=adprminj(j),adprminj(j+1)-1
         if(nameprminj(i).ne.'')then
            k=k+1
            f_name=nameprminj(i)
            do n=1,10
               names(11*(k-1)+n)=f_name(n:n)
            enddo
            names(11*k)=' '
         endif
      enddo

   case ('TWOP')
      call seart(fcomp_name,j)
      if(j == 0)then
         write(msg,"('t = ',f9.4,'  unknown twoport ',a)")t_h(0),trim(fcomp_name)
         call write_msg('get_comp_prm_names',trim(msg))
         get_comp_prm_names=1
         return
      endif
      do i=adprmtwop(j),adprmtwop(j+1)-1
         if(nameprmtwop(i).ne.'')then
            k=k+1
            f_name=nameprmtwop(i)
            do n=1,10
               names(11*(k-1)+n)=f_name(n:n)
            enddo
            names(11*k)=' '
         endif
      enddo

   case ('DCTL')
      call seard(fcomp_name,j)
      if(j == 0)then
         write(msg,"('t = ',f9.4,'  unknown DCTL ',a)")t_h(0),trim(fcomp_name)
         call write_msg('get_comp_prm_names',trim(msg))
         get_comp_prm_names=1
         return
      endif
      do i=adwdctl(j),adwdctl(j+1)-1
         if(namewdctl(i).ne.'')then
            k=k+1
            f_name=namewdctl(i)
            do n=1,10
               names(11*(k-1)+n)=f_name(n:n)
            enddo
            names(11*k)=' '
         endif
      enddo
end select
if (k>0) names(11*k)=c_null_char
get_comp_prm_names = 0
return

end function get_comp_prm_names

!> @brief Return shunt name
!> @details .
integer(c_int) function get_shunt_name(i,name) bind(C, name="get_shunt_name")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_shunt_name
#endif
use SHUNT, only: nbshunt, shuname
implicit none
integer(c_int), INTENT(IN), value :: i
character(c_char), dimension(21), intent(inout) :: name
character(len=21) :: f_name
integer :: n, j
name(1)=C_NULL_CHAR
if(i<=nbshunt .and. i>0)then
   f_name=trim(shuname(i))
   n = len_trim(f_name)
   do j = 1, n
      name(j) = f_name(j:j)
   end do
   name(n + 1) = c_null_char
   get_shunt_name = 0
   return
endif
get_shunt_name = 1
end function get_shunt_name

!> @brief Return nbload
!> @details .
integer(c_int) function get_nbload() bind(C, name="get_nbload")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_nbload
#endif
use LOAD, only: nbload
implicit none
get_nbload=nbload
end function get_nbload

!> @brief Return load name
!> @details .
integer(c_int) function get_load_name(i,name) bind(C, name="get_load_name")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_load_name
#endif
use LOAD, only: nbload, loadname
implicit none
integer(c_int), INTENT(IN), value :: i
character(c_char), dimension(21), intent(inout) :: name
character(len=21) :: f_name
integer :: n, j
name(1)=c_null_char
if(i<=nbload .and. i>0)then
   f_name=trim(loadname(i))
   n = len_trim(f_name)
   do j = 1, n
      name(j) = f_name(j:j)
   end do
   name(n + 1) = c_null_char
   get_load_name = 0
   return
endif
get_load_name = 1
end function get_load_name

!> @brief Set pause time
!> @details .
integer(c_int) function set_pause_time(c_pause_time) bind(C, name="set_pause_time")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: set_pause_time
#endif
   use SimTIME, only: pause_time
   implicit none
   real(c_double), intent(in), value :: c_pause_time
   !$omp atomic write
   pause_time=c_pause_time
   set_pause_time=0
end function set_pause_time

!> @brief Continue simulation
!> @details .
integer(c_int) function continue_simul() bind(C, name="continue_simul")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: continue_simul
#endif
   use SETTINGS,only: error_flag
   implicit none
   error_flag=.false.
   continue_simul=ramses()
end function continue_simul

!> @brief Activate backup
!> @details .
integer(c_int) function set_to_backup() bind(C, name="set_to_backup")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: set_to_backup
#endif
   use backup_data, only: to_back_up
   implicit none
   to_back_up=.true.
   set_to_backup=0
end function set_to_backup

!> @brief End the simulation
!> @details .
integer(c_int) function set_end_simul() bind(C, name="set_end_simul")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: set_end_simul
#endif
use SETTINGS, only: end_simul
implicit none
!$omp atomic write
end_simul=.true.
set_end_simul=0
end function set_end_simul

!> @brief Check if simulation has ended
!> @details .
integer(c_int) function get_end_simul() bind(C, name="get_end_simul")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_end_simul
#endif
use SETTINGS, only: end_simul
implicit none

if(end_simul .eq. .true.)then
!$omp atomic write
   get_end_simul=1
else
!$omp atomic write
   get_end_simul=0
endif
end function get_end_simul

!> @brief This is a C wrapper for adding a disturbance
!> @details .
integer(c_int) function add_disturb(disturbtime,disturbdesc) bind(C, name="add_disturb")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: add_disturb
#endif
use DISTURB
use SETTINGS, only: write_msg
use SimTIME, only: t_h
use simul_decomp_mod, only: t_horizon
implicit none
real(c_double), intent(in), value :: disturbtime
character(c_char), dimension(*), intent(in) :: disturbdesc
character(len=200) :: msg
character(len=256) :: fdisturbdesc
double precision :: t_dist_new

integer :: i,j

fdisturbdesc=c_to_f_string(disturbdesc)
t_dist_new=disturbtime

write(msg,"('PyRAMSES: Adding disturbance ', f15.6,' ', a)")t_dist_new, trim(fdisturbdesc)
call write_msg('add_disturb',msg)

if(nbdist == mxdist)then
   write(msg,"('more than ',i5,' disturbances')")mxdist
   call write_msg('add_disturb',msg)
   add_disturb=1
   return
endif

if(t_dist_new < t_h(0) .or. t_dist_new < 0.d0)then
   write(msg,"('you are trying to add a disturbance in the past')")
   call write_msg('add_disturb',msg)
   add_disturb=2
   return
endif

if (trim(fdisturbdesc) .eq. 'STOP')then
   t_dist(nbdist)=t_dist_new
   add_disturb=0
   return
endif

if(t_dist_new > t_dist(nbdist))then
   write(msg,"('you are trying to add a disturbance after the STOP time')")
   call write_msg('add_disturb',msg)
   add_disturb=3
   return
endif

distloop: &
do i=1, nbdist
   if(t_dist(i) > t_dist_new)then
      do j=nbdist+1,i+1 ! shift disturbances
         t_dist(j)=t_dist(j-1)
         desc_dist(j)=desc_dist(j-1)
      enddo
      t_dist(i)=t_dist_new
      desc_dist(i)=trim(fdisturbdesc)
      exit distloop
   endif
enddo distloop

if (t_horizon > t_dist_new)t_horizon=t_dist_new

nbdist=nbdist+1

add_disturb=0
end function add_disturb


!> @brief Define a subsystem as a list of buses
!> @details .
integer(c_int) function define_SS(ssID,filter1,num1,filter2,num2,filter3,num3) bind(C, name="define_SS")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: define_SS
#endif
use search_mod
use SETTINGS, only: write_msg
use BUS, only: nbbus, vnom, busname
use ZONES
use tokenize
implicit none

integer(c_int), INTENT(in), value :: ssID, num1, num2, num3
character(c_char), dimension(*), intent(in) :: filter1,filter2,filter3
character(len=:), allocatable :: ffilter1,ffilter2,ffilter3
character(len=200) :: msg
integer :: j, i, k, n, length, flag
double precision :: kv
type(tokenizer)   :: token
character(len=20) :: tempStr
type(setType) :: tmpSet1, tmpSet2, tmpSet3, tmpSet

if (ssID<=0 .or. ssID > mxzon) then
   define_SS = -1
   return
end if

call init_set(subsystem(ssID), nbbus) ! initialize the set to the proper size
call init_set(tmpSet1, nbbus) ! initialize the set to the proper size
call init_set(tmpSet2, nbbus) ! initialize the set to the proper size
call init_set(tmpSet3, nbbus) ! initialize the set to the proper size
call init_set(tmpSet, nbbus) ! initialize the set to the proper size

call set_tokenizer( token, token_whitespace//token_tsv, token_empty, token_quotes ) ! setup the tokenizer for spaces and tabs and quotes

! Voltage level filter
if (num1>0) then
   ffilter1=c_to_f_string(filter1)
   tempStr=first_token( token, ffilter1, length )

   if(length<=0) then
      call write_msg('define_SS','Error with filter 1')
      define_SS = -1
      return
   end if

   do while ( length .ge. 0 )
      if(length .gt. 20)then
         call write_msg('define_SS','An item in filter 1 has a name bigger than 20')
         define_SS = -1
         return
      endif
      call str2dp(tempStr, kv, flag)
      if ( flag == 0 ) then
         do i=1, nbbus
            if (vnom(i) == kv) call set_member( tmpSet1, i )
!            print *, i, busname(i), kv, vnom(i)
         end do
      else
         call write_msg('define_SS','An item in filter 1 cannot be converted to integer')
         define_SS = -1
         return
      endif
      tempStr = next_token( token, ffilter1, length )
   enddo
else
   call set_all_elem( tmpSet1 )
end if

! Zone
if (num2>0) then
   ffilter2=c_to_f_string(filter2)
   tempStr=first_token( token, ffilter2, length )
   if(length<=0) then
      call write_msg('define_SS','Error with filter 2')
      define_SS = -1
      return
   end if
   do while ( length .ge. 0 )
      if(length .gt. 20)then
         call write_msg('define_SS','An item in filter 2 has a name bigger than 20')
         define_SS = -1
         return
      endif
      call searzon(tempStr,n)
      if (n<=0) then
         call write_msg('define_SS','The zone in Filter 2 was not found')
         define_SS = -1
         return
      end if
      do i=adzonbus(n), adzonbus(n+1)-1
         call set_member(tmpSet2, zonbus(i) )
      end do
      tempStr = next_token( token, ffilter2, length )
   enddo
else
   call set_all_elem( tmpSet2 )
end if

! Bus name filter
if (num3>0) then
   ffilter3=c_to_f_string(filter3)
   tempStr=first_token( token, ffilter3, length )
   if(length<=0) then
      call write_msg('define_SS','Error with filter 3')
      define_SS = -1
      return
   end if
   do while ( length .ge. 0 )
      if(length .gt. 20)then
         call write_msg('define_SS','An item in filter 3 has a name bigger than 20')
         define_SS = -1
         return
      endif
      call searn(tempStr(1:18),n)
      if (n<=0) then
         call write_msg('define_SS','A bus in Filter 3 was not found')
         define_SS = -1
         return
      end if
      call set_member(tmpSet3,n)
      tempStr = next_token( token, ffilter3, length )
   enddo
else
   call set_all_elem( tmpSet3 )
end if

call set_intersection( tmpSet1, tmpSet2, subsystem(ssID) )
tmpSet1 = subsystem(ssID)
call set_intersection( tmpSet1, tmpSet3, subsystem(ssID) )

define_SS = 0
end function define_SS

!> @brief Retrieve a subsystem as a list of buses
!> @details .
integer(c_int) function get_SS(ssID,mxreclen,busNames) bind(C, name="get_SS")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_SS
#endif
use search_mod
use SETTINGS, only: write_msg
use SYNC, only: syncname
use UDIM, only: injname
use TWOP, only: twopname
use DCTL, only: dctlname
use BUS, only: busname,nbbus
use SimTIME, only: t_h
use tokenize
implicit none

integer(c_int), INTENT(in), value :: ssID, mxreclen
character(c_char), dimension(19*mxreclen), intent(inout) :: busNames
integer :: i, n, k
character(len=200) :: msg
character(len=20) :: f_name

if (ssID<=0 .or. ssID > mxzon .or. .not. set_is_valid( subsystem(ssID) )) then
   call write_msg('get_SS','Error with the ssID, either does not exist or not initialized')
   get_SS = -1
   return
end if

busNames(1) = c_null_char

k=0
do i = 1, nbbus
   if (set_has_member( subsystem(ssID), i ))then
      k = k+1
      f_name=''
      f_name=trim(busname(i))
      do n=1,18
         busNames(18*(k-1)+n)=f_name(n:n)
      enddo
      busNames(19*k)=' '
   end if
end do

if (k>0) busNames(19*k)=c_null_char

get_SS = 0
end function get_SS

!> @brief Retrieve list of transformers in SS
!> @details .
integer(c_int) function get_transformer_in_SS(ssID, location, in_service, rettype, mxreclen, retStr, retdp, retint, elem) bind(C, name="get_transformer_in_SS")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_transformer_in_SS
#endif
use SETTINGS, only: write_msg
use BRANCH, only: nbbra,braname,brabr_extr,brabr_orig,bratype, origin, extrem
use BUS, only: busname, vnom
use sets
use observ_mod, only: pqbra
use NET_TOPO, only: sbases, bussubnet
implicit none

integer(c_int), INTENT(in), value :: ssID, mxreclen, location, in_service
character(c_char), dimension(21*mxreclen), intent(inout) :: retStr
real(c_double), dimension(mxreclen), intent(inout) :: retdp
integer(c_int), dimension(mxreclen), intent(inout) :: retint
integer(c_int), intent(inout) :: elem
integer :: i, n, k
double precision :: p_orig, q_orig, p_extr, q_extr
character(len=200) :: msg
character(len=20) :: f_name
character(c_char), dimension(*), intent(in) :: rettype
character(len=10) :: frettype
type(setType) :: trfoSet

call init_set(trfoSet, nbbra) ! initialize the set to the proper size

frettype=c_to_f_string(rettype)

if (ssID<=0 .or. ssID > mxzon .or. .not. set_is_valid( subsystem(ssID) )) then
   call write_msg('get_transformer_in_SS','Error with the ssID, either does not exist or not initialized')
   get_transformer_in_SS = -1
   return
end if

retStr(1) = c_null_char
retdp(1:mxreclen) = 0.d0
retint(1:mxreclen) = 0
elem = 0

do i=1, nbbra
   if ( bratype(i) == 'trfo' ) then

      if ( location == 1 .or. location == 3) then
         if (set_has_member( subsystem(ssID), origin(i) ) .and. set_has_member( subsystem(ssID), extrem(i) )) &
          call set_member ( trfoSet, i )
      end if
      if ( location == 2 .or. location == 3) then
         if ( ( set_has_member( subsystem(ssID), origin(i) ) .and. .not. set_has_member( subsystem(ssID), extrem(i) )) &
               .or. ( .not. set_has_member( subsystem(ssID), origin(i) ) .and. set_has_member( subsystem(ssID), extrem(i) ))) &
          call set_member ( trfoSet, i )
      end if

      if ( In_service == 1 .and. (brabr_orig(i) == 0 .or. brabr_extr(i) == 0) ) call set_clear_member ( trfoSet, i )

   end if
end do

k=0
do i=1, nbbra
   if ( set_has_member ( trfoSet, i ) )then
      k = k+1
      if (k >= mxreclen) then
         call write_msg('get_transformer_in_SS','The number of returning values is higher than mxreclen')
         get_transformer_in_SS = -1
         return
      end if
      f_name=''
      select case(frettype)
         case ('NAME')
            f_name=trim(braname(i))
            do n=1,20
               retStr(20*(k-1)+n)=f_name(n:n)
            enddo
            retStr(21*k)=' '
         case ('From')
            f_name=trim(busname(origin(i)))
            do n=1,20
               retStr(20*(k-1)+n)=f_name(n:n)
            enddo
            retStr(21*k)=' '
         case ('To')
            f_name=trim(busname(extrem(i)))
            do n=1,20
               retStr(20*(k-1)+n)=f_name(n:n)
            enddo
            retStr(21*k)=' '
         case ('Status')
            if (brabr_orig(i)==1 .and. brabr_extr(i)==1) retint(k) = 1
         case ('Tap')
!           To be discussed. We don't save the tap number!
         case ('Currentf')
            call pqbra(i,p_orig,q_orig,p_extr,q_extr,.false.)
            retdp(k) = hypot(p_orig, q_orig) * sbases(bussubnet(origin(i)))/(vnom(origin(i))*dsqrt(3.0d0)) ! in kA
         case ('Currentt')
            call pqbra(i,p_orig,q_orig,p_extr,q_extr,.false.)
            retdp(k) = hypot(p_extr, q_extr) * sbases(bussubnet(extrem(i)))/(vnom(extrem(i))*dsqrt(3.0d0)) ! in kA
         case ('Pf')
            call pqbra(i,p_orig,q_orig,p_extr,q_extr,.true.)
            retdp(k) = p_orig * sbases(bussubnet(origin(i))) ! in MW
         case ('Qf')
            call pqbra(i,p_orig,q_orig,p_extr,q_extr,.true.)
            retdp(k) = q_orig * sbases(bussubnet(origin(i))) ! in MVA
         case ('Pt')
            call pqbra(i,p_orig,q_orig,p_extr,q_extr,.true.)
            retdp(k) = p_extr * sbases(bussubnet(extrem(i))) ! in MW
         case ('Qt')
            call pqbra(i,p_orig,q_orig,p_extr,q_extr,.true.)
            retdp(k) = q_extr * sbases(bussubnet(extrem(i))) ! in MVA
      end select
   end if
end do
if (k>0) retStr(21*k)=c_null_char
elem = k
get_transformer_in_SS = 0
end function get_transformer_in_SS

!> @brief Initialize observable selection (structure and trajectory file)
!> @details .
integer(c_int) function initObserv(traj_filenm) bind(C, name="initObserv")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: initObserv
#endif
use SETTINGS, only: error_flag
use observ_mod, only: observ_init
implicit none

character(c_char), dimension(*), intent(in) :: traj_filenm
character(len=256) :: ftraj_filenm

ftraj_filenm=c_to_f_string(traj_filenm)
call observ_init(ftraj_filenm)
if(error_flag)then
   initObserv=1
   return
end if
initObserv=0
end function initObserv

!> @brief Add an element to be observed
!> @details .
integer(c_int) function addObserv(string) bind(C, name="addObserv")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: addObserv
#endif
use SETTINGS, only: warn_flag
use observ_mod, only: add_observ
implicit none

character(c_char), dimension(*), intent(in) :: string
character(len=256) :: fstring

addObserv=0
fstring=c_to_f_string(string)
warn_flag=.false.
call add_observ(fstring)
if(warn_flag)then
   addObserv=1
   return
end if
end function addObserv

!> @brief Finalize observable selection, allocate buffer, and write header of trajectory file
!> @details .
integer(c_int) function finalObserv() bind(C, name="finalObserv")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: finalObserv
#endif
use SETTINGS, only: error_flag
use observ_mod, only: observ_fin
implicit none

call observ_fin()
if(error_flag)then
   finalObserv=1
   return
end if
finalObserv=0
end function finalObserv

!=================================================================================================================
!Single Element Data Retrieval
!=================================================================================================================
!> @brief Single Element Data Retrieval of bus real values (To complete with int and complex value)
!> @details .
integer(c_int) function get_bus_data(bus_name, rettype, bus_value) bind(C, name="get_bus_data")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_bus_data
#endif
!   bus_name: bus bame
!   rettype   : bus quantity to get
!   bus_value : returned real vector
!
!   returns   : errorflag
!
use BUS, only: nbbus, vnom, busname
use VOLTAGE, only: vx_h, vy_h
use SETTINGS, only: rad
use search_mod, only: searn
implicit none
doubleprecision :: vx, vy
character(c_char), dimension(*), intent(in) :: bus_name
character(c_char), dimension(*), intent(in) :: rettype
real(c_double), intent(inout) :: bus_value
character(len=20) :: fname
character(len=:), allocatable :: frettype
integer :: bus_number

fname=c_to_f_string(bus_name)
call searn(fname(1:18),bus_number)
if (bus_number==0) then
   bus_value=0.d0
   get_bus_data=1
   return
endif

frettype=c_to_f_string(rettype)

select case(frettype)
    case ('BASE')
        bus_value=vnom(bus_number)
    case ('PU')
        bus_value=hypot(vx_h(bus_number), vy_h(bus_number))
    case ('KV')
        bus_value=hypot(vx_h(bus_number), vy_h(bus_number))*vnom(bus_number)
    case ('ANGLE')
        bus_value=datan2(vy_h(bus_number), vx_h(bus_number))
    case ('ANGLED')
        bus_value=datan2(vy_h(bus_number), vx_h(bus_number))*rad
        case DEFAULT
        get_bus_data = 2
        return
end select

get_bus_data=0
end function get_bus_data

!=================================================================================================================
!> @brief Single Element Data Retrieval of branch real data
!> @details .
integer(c_int) function get_branch_real_data(branch_name, rettype, retdp) bind(C, name="get_branch_real_data")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_branch_real_data
#endif
use SETTINGS, only: write_msg, rad
use BRANCH, only: nbbra,braname,brabr_extr,brabr_orig, origin, extrem, phan, bratype
use BUS, only: busname, vnom
use search_mod, only: searb
use observ_mod, only: pqbra
use NET_TOPO, only: sbases, bussubnet
implicit none

real(c_double), intent(inout) :: retdp
integer :: i, ii, j, k, n_str, stat
double precision :: p_orig, q_orig, p_extr, q_extr
character(len=20) :: f_name
character(c_char), dimension(*), intent(in) :: branch_name
character(c_char), dimension(*), intent(in) :: rettype
type(setType) :: trfoSet
character(len=:), allocatable :: frettype

character(len=20) :: fname
integer :: branch_num

branch_num = 0
fname=c_to_f_string(branch_name)
!write(log,"(a)") fname
call searb(fname(1:20),branch_num)
!write(log,"(i)") bra_num
if (branch_num==0) then
   get_branch_real_data=1
   return
endif

if ( bratype(branch_num) == 'trfo') then
    get_branch_real_data=1
    return
endif


frettype=c_to_f_string(rettype)

select case(frettype)
    case ('AMPS') !FROM kA
        ii=branch_num
        call pqbra(ii,p_orig,q_orig,p_extr,q_extr,.false.)  !returns current .false.
        retdp = hypot(p_orig, q_orig) * sbases(bussubnet(origin(ii)))/(vnom(origin(ii))*dsqrt(3.0d0)) ! in kA
    case ('PUCUR') !FROM pu
        ii=branch_num
        call pqbra(ii,p_orig,q_orig,p_extr,q_extr,.false.)  !returns current .false.
        retdp = hypot(p_orig, q_orig) ! in pu
    case ('P') !FROM BUS MW
        ii=branch_num
        call pqbra(ii,p_orig,q_orig,p_extr,q_extr,.true.)  !returns power .true.
        retdp = p_orig * sbases(bussubnet(origin(ii))) ! in MW
    case ('Q') !FROM BUS MVAR
        ii=branch_num
        call pqbra(ii,p_orig,q_orig,p_extr,q_extr,.true.)  !returns power .true.
        retdp = q_orig * sbases(bussubnet(origin(ii))) ! in MVA
    case ('MVA') !FROM BUS MVA
        ii=branch_num
        call pqbra(ii,p_orig,q_orig,p_extr,q_extr,.true.)  !returns power .true.
        retdp = hypot(p_orig, q_orig) * sbases(bussubnet(origin(ii))) ! in MVA
        case DEFAULT
        get_branch_real_data = 5
        return
    end select

get_branch_real_data = 0
end function get_branch_real_data

!=================================================================================================================
!> @brief Single Element Data Retrieval of branch integer data
!> @details .
integer(c_int) function get_branch_int_data(branch_name, rettype, retint) bind(C, name="get_branch_int_data")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_branch_int_data
#endif
use SETTINGS, only: write_msg, rad
use BRANCH, only: nbbra,braname,brabr_extr,brabr_orig, origin, extrem, phan, bratype
use BUS, only: busname, vnom
use search_mod, only: searb
use observ_mod, only: pqbra
use NET_TOPO, only: sbases, bussubnet
implicit none

integer(c_int), intent(inout) :: retint
integer :: i, ii, j, k, n_str, stat
double precision :: p_orig, q_orig, p_extr, q_extr
character(len=20) :: f_name
character(c_char), dimension(*), intent(in) :: branch_name
character(c_char), dimension(*), intent(in) :: rettype
character(len=:), allocatable :: frettype

character(len=20) :: fname
integer :: branch_num

branch_num = 0
fname=c_to_f_string(branch_name)
!write(log,"(a)") fname
call searb(fname(1:20),branch_num)
!write(log,"(i)") bra_num
if (branch_num==0) then
    get_branch_int_data=2
    return
endif

frettype=c_to_f_string(rettype)

select case(frettype)
    case ('STATUS')
        ii=branch_num
        if (brabr_orig(ii)==1 .and. brabr_extr(ii)==1)then
            retint = 1
        else
            retint = 0
        endif
    case DEFAULT
        get_branch_int_data = 3
        return
end select

get_branch_int_data = 0
end function get_branch_int_data

!=================================================================================================================
!> @brief Single Element Data Retrieval of xfo real data
!> @details .
integer(c_int) function get_xfo_real_data(branch_name, rettype, retdp) bind(C, name="get_xfo_real_data")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_xfo_real_data
#endif
use SETTINGS, only: write_msg, rad
use BRANCH, only: nbbra,braname,brabr_extr,brabr_orig, origin, extrem, phan, bratype
use BUS, only: busname, vnom
use search_mod, only: searb
use observ_mod, only: pqbra
use NET_TOPO, only: sbases, bussubnet
implicit none

real(c_double), intent(inout) :: retdp
integer :: i, ii, j, k, n_str, stat
double precision :: p_orig, q_orig, p_extr, q_extr
character(len=20) :: f_name
character(c_char), dimension(*), intent(in) :: branch_name
character(c_char), dimension(*), intent(in) :: rettype
type(setType) :: trfoSet
character(len=:), allocatable :: frettype

character(len=20) :: fname
integer :: branch_num

branch_num = 0
fname=c_to_f_string(branch_name)
!write(log,"(a)") fname
call searb(fname(1:20),branch_num)
!write(log,"(i)") bra_num
if (branch_num==0) then
   get_xfo_real_data=1
   return
endif

if (bratype(branch_num) /= 'trfo') then
    get_xfo_real_data=1
    return
endif

frettype=c_to_f_string(rettype)

select case(frettype)
    case ('ANGLE')
            ii=branch_num
            retdp = phan(ii)*rad !ideal transformer phase angle in degrees
        case DEFAULT
        get_xfo_real_data = 5
        return
end select

get_xfo_real_data = 0
end function get_xfo_real_data

!=================================================================================================================
!> @brief Single Element Data Retrieval of sync integer data
!> @details .
integer(c_int) function get_sync_int_data(sync_name, rettype, retint) bind(C, name="get_sync_int_data")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_sync_int_data
#endif
use SETTINGS, only: write_msg, rad
use SYNC, only: nbsync,syncname
use INJ, only: injbr
use search_mod, only: searm
implicit none

integer(c_int), intent(inout) :: retint
integer :: i, ii, j, k, n_str, stat
character(len=20) :: f_name
character(c_char), dimension(*), intent(in) :: sync_name
character(c_char), dimension(*), intent(in) :: rettype
character(len=:), allocatable :: frettype

character(len=20) :: fname
integer :: sync_num

sync_num = 0
fname=c_to_f_string(sync_name)
call searm(fname(1:20),sync_num)
if (sync_num==0) then
   get_sync_int_data=2
   return
endif


frettype=c_to_f_string(rettype)

select case(frettype)
    case ('STATUS')
        retint=injbr(sync_num)
    case DEFAULT
        get_sync_int_data = 5
        return
end select

get_sync_int_data = 0
end function get_sync_int_data

!=================================================================================================================
!> @brief Single Element Data Retrieval of sync real data
!> @details .
integer(c_int) function get_sync_real_data(sync_name, rettype, retdp) bind(C, name="get_sync_real_data")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_sync_real_data
#endif
use SETTINGS, only: write_msg, rad
use SYNC, only: nbsync,syncname,snom_sync,adxsync,xsync_h
use INJ, only: injbr
use search_mod, only: searm
use observ_mod, only: pqsync
use NET_TOPO, only: sbases, bussubnet
implicit none

real(c_double), intent(inout) :: retdp
integer :: i, ii, j, k, n_str, stat
double precision :: p_sync, q_sync
character(len=20) :: f_name
character(c_char), dimension(*), intent(in) :: sync_name
character(c_char), dimension(*), intent(in) :: rettype
character(len=:), allocatable :: frettype

character(len=20) :: fname
integer :: sync_num

sync_num = 0
fname=c_to_f_string(sync_name)
!write(log,"(a)") fname
call searm(fname(1:20),sync_num)
!write(log,"(i)") bra_num
if (sync_num==0) then
   get_sync_real_data=1
   return
endif

frettype=c_to_f_string(rettype)

select case(frettype)
    case ('MBASE')
        ii=sync_num
        retdp = snom_sync(ii)
    case ('P') !SYNC MW
        ii=sync_num
        call pqsync(ii,p_sync,q_sync)  !returns p & q
        retdp = p_sync * snom_sync(ii) ! in MW
    case ('Q') !SYNC MVAR
        ii=sync_num
        call pqsync(ii,p_sync,q_sync)  !returns p & q
        retdp = q_sync * snom_sync(ii) ! in MVAR
    case ('MVA') !SYNC MVA
        ii=sync_num
        call pqsync(ii,p_sync,q_sync)  !returns p & q
        retdp = hypot(p_sync, q_sync) * snom_sync(ii) ! in MVA
    case DEFAULT
        get_sync_real_data = 5
        return
end select

get_sync_real_data = 0
end function get_sync_real_data

!=================================================================================================================
!> @brief Single Element Data Retrieval of load complex data
!> @details .
integer(c_int) function get_load_cplx_data(load_name, rettype1, rettype2, retdp) bind(C, name="get_load_cplx_data")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_load_cplx_data
#endif
!   ssID      : Id of subsystem
!   flag      : Load    Bus   (INS = inservice)
!          =1    INS    INS
!          =2    ALL    INS
!          =3    INS    ALL
!          =4    ALL    ALL
!   rettype   : tokenized string - list of values to get
!   mxreclen  : vector length for returned variables retdp & retint
!   retdp     : returned real vector
!   nb_load_ss : number of elements returned for each asked value in rettype
!
!   returns   : errorflag
!
use BUS, only: nbbus, busname
use LOAD
use VOLTAGE, only: vx_h, vy_h
use DIMENSIONS, only :mxinj
use INJ, only: injbr
use search_mod, only: seari
use observ_mod, only: pqsync

implicit none
integer :: i, ii, j, n, k, n_str, ierr
character(c_char), dimension(*), intent(in) :: load_name,rettype1, rettype2
real(c_double), dimension(2), intent(out) :: retdp
character(len=10) :: fobs_name
character(len=:),allocatable :: rettype
double precision :: obs_value
logical :: flg_act

character(len=20) :: fname, frettype1, frettype2
integer :: load_num

load_num = 0
fname=c_to_f_string(load_name)
!write(log,"(a)") fname
call seari(fname(1:20),load_num)
!write(log,"(i)") bra_num
if (load_num==0) then
   get_load_cplx_data=2
   return
endif

frettype1 = c_to_f_string(rettype1)
frettype2 = c_to_f_string(rettype2)

rettype = trim(frettype1) // trim(frettype2)

get_load_cplx_data = 0

select case(rettype)
    case ('TOTALACT') !complex value - doubled in retdp !???? same name as real
        ii=load_num
        fobs_name='P'
        ierr = get_named_obs_inj(ii, fobs_name, retdp(1), .true., .true.)
        fobs_name='Q'
        ierr = get_named_obs_inj(ii, fobs_name, retdp(2), .true., .true.)
    case DEFAULT
        get_load_cplx_data = 5
    return
    end select

get_load_cplx_data = 0
end function get_load_cplx_data

!=====================================================================================
!> @brief Get Apparent Power Base of System
!> @details .
real(c_double) function get_mva_base() bind(C, name="get_mva_base")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: get_mva_base
#endif
use SETTINGS, only: sbasetransm
implicit none
!!$omp atomic read
get_mva_base=sbasetransm
end function get_mva_base


!=====================================================================================
!> @brief write message from python to output file
!> @details .
integer(c_int) function send_msg_to_output(string_msg) bind(C, name="send_msg_to_output")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: send_msg_to_output
#endif
   use SETTINGS, only: write_msg
   implicit none

   character(c_char), dimension(*), intent(in) :: string_msg
   character(len=:), allocatable :: f_str_msg

   f_str_msg=c_to_f_string(string_msg)

   call write_msg('send_msg_to_output',f_str_msg)
   send_msg_to_output=0
   return
end function send_msg_to_output

!=====================================================================================
!> @brief load user model DLL
!> @details .
integer(c_int) function c_load_MDL(string_msg) bind(C, name="c_load_MDL")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: c_load_MDL
#endif
   use MODELING, only: load_MDL
   implicit none

   character(c_char), dimension(*), intent(in) :: string_msg
   character(len=:), allocatable :: f_str_msg

   f_str_msg=c_to_f_string(string_msg)

   call load_MDL(f_str_msg)
   c_load_MDL=0
   return
end function c_load_MDL

!=====================================================================================
!> @brief unload user model DLL
!> @details .
integer(c_int) function c_unload_MDL(string_msg) bind(C, name="c_unload_MDL")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: c_unload_MDL
#endif
   use MODELING, only: unload_MDL
   implicit none

   character(c_char), dimension(*), intent(in) :: string_msg
   character(len=:), allocatable :: f_str_msg

   f_str_msg=c_to_f_string(string_msg)

   call unload_MDL(f_str_msg)
   c_unload_MDL=0
end function c_unload_MDL

!=====================================================================================
!> @brief returns number of user model DLL loaded
!> @details .
integer(c_int) function c_get_MDL_no() bind(C, name="c_get_MDL_no")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: c_get_MDL_no
#endif
   use MODELING, only: dll_handleno,RAMSESMDLfile
   implicit none
   integer i,n

   n=0
   do i=1,dll_handleno
      if(RAMSESMDLfile(i)/="")n=n+1
   enddo

   c_get_MDL_no=n
end function c_get_MDL_no

!=====================================================================================
!> @brief returns list of user model DLL loaded
!> @details .
integer(c_int) function c_get_MDL_names(mxreclen,DLL_Names) bind(C, name="c_get_MDL_names")
#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT :: c_get_MDL_names
#endif
   use MODELING, only: dll_handleno, RAMSESMDLfile
   implicit none
   integer(c_int), INTENT(in), value :: mxreclen
   character(c_char), dimension(mxreclen), intent(inout) :: DLL_Names
   character(len=mxreclen):: DLL_name
   integer ::i,j,k,n

   DLL_Names(1) = c_null_char

   n=0
   k=1
   do i=1,dll_handleno
      if(RAMSESMDLfile(i)/="")then
         if(n>0)then
            DLL_Names(k:k)="|"
            k=k+1
         endif
         n=n+1
         do j = 1, len_trim(RAMSESMDLfile(i))
            DLL_Names(k) = RAMSESMDLfile(i)(j:j)
            k=k+1
         end do
      endif
   enddo

   if (n>0) then
      DLL_names(k)=c_null_char
   endif
   c_get_MDL_names=n
end function c_get_MDL_names

end module c_interface_mod
