!> @file
!> @brief Implements the vfault subroutines
!! LOAD inj_name bus_name FP FQ P Q  DP A1 alpha1 A2 alpha2 alpha3 DQ B1 beta1 B2 beta2 beta3 ;
!! prm=( DP A1 alpha1 A2 alpha2 alpha3 DQ B1 beta1 B2 beta2 beta3 )
!! prm=(  1  2   3    4    5        6   7  8   9   10   11    12  )
subroutine inj_vfault(nb,name,mode,nbxvar,nbzvar,nbdata,nbaddpar,prm,parname,nbobs, &
      obsname,adix,adiy,eqtyp,eqtyp_tc,t,omega,sbase,busnum,vx,vy,ix,iy,x,z,f,obs)

   use MODELING
   use DISTURB, only: vfault_val, vfault_inj_num, vfault_v_pre, vfault_bfault_val, vfault_bfault_val_pre
   use INJ, only: injbr, bus_inj
   use BUS, only: gfault, bfault, busfault
   use SYNC, only: nbsync
   use SETTINGS, only: write_msg_and_stop

   implicit none
   double precision, intent(in):: t,vx,vy,omega,sbase,ix,iy
   double precision, intent(out):: f(*)
   double precision :: obs(*)
   double precision, intent(inout):: x(*),prm(*),eqtyp_tc(*)
   integer, intent(in):: nb,mode,busnum
   integer, intent(inout):: nbxvar,nbzvar,nbdata,nbaddpar,nbobs,eqtyp(*),z(*),adix,adiy
   character(len=20), intent(in):: name
   character(len=10) :: parname(*),obsname(*)

   double precision :: tmp, alpha, beta

   select case (mode)
      case (define_var_and_par)
         nbxvar=2
         nbzvar=1
         nbdata=0
         nbaddpar=0
         adiy=1
         adix=2
         vfault_inj_num=nb+nbsync
         return
      case (define_obs)
         nbobs=1
         obsname(1)='Xfault'
         return
      case (initialize)
         z(1)=0
         vfault_bfault_val_pre=0.d0
         vfault_bfault_val=0.d0

         eqtyp(1)=0
         eqtyp(2)=0
         injbr(vfault_inj_num)=0

      case (evaluate_eqs)
         f(1)=vfault_bfault_val*vx-x(1)
         f(2)=-vfault_bfault_val*vy-x(2)

      case (update_disc)
         if(abs(vfault_val-hypot(vx,vy))<=0.005d0 .and. z(1)<=0)then
            z(1)=1
!            print *, 't=', t, 'Final B=', vfault_bfault_val,'DeltaV=',abs(vfault_val-hypot(vx,vy))
         elseif(z(1)<=-10)then
            call write_msg_and_stop('VFAULT: ','Tried more than 10 times to reach voltage set-point. Exiting...')
            return
         elseif(z(1)<=0)then
            z(1)=z(1)-1
            alpha=(hypot(vx,vy)-vfault_v_pre)/(vfault_bfault_val-vfault_bfault_val_pre)
            beta=vfault_v_pre-alpha*vfault_bfault_val_pre
            tmp=(vfault_val-beta)/alpha
            vfault_bfault_val_pre=vfault_bfault_val
            vfault_bfault_val=tmp
            vfault_v_pre=hypot(vx,vy)
!            print *, 't=', t, 'New B=', vfault_bfault_val,'DeltaV=',abs(vfault_val-hypot(vx,vy))
         endif

      case (evaluate_obs)
         obs(1)=vfault_bfault_val
         return
   end select

end subroutine inj_vfault

