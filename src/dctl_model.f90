!> @file
!! @brief discrete controller subroutines

!> @brief subroutine to define discrete controler
!> @details .
!> @param[in] modelname name of the model
!> @param[in] name name of the controller
!> @param[in] field field of record
!> @param[in] nbfields
!> @param[in] nbwvar
!> @param[in] parname
!> @param[in] w working variables
subroutine def_eq_dctl_model(nb,modelname,name,field,w,nbfields,nbwvar,parname)

   use UNITS
   use SETTINGS, only: write_msg_and_stop

   implicit none

   double precision, intent(out):: w(*)
   integer, intent(out):: nbwvar
   character(len=20), intent(in):: modelname,name,field(*)
   character(len=10), intent(out):: parname(*)
   integer, intent(in) :: nb,nbfields

   select case (modelname)

      case('PST')
         call def_eq_dctl_pst(nb,name,field,w,nbfields,nbwvar,parname)

      case('LTC')
         call def_eq_dctl_ltc(nb,name,field,w,nbfields,nbwvar,parname)

      case('LTC2')
         call def_eq_dctl_ltc2(nb,name,field,w,nbfields,nbwvar,parname)

      case('LTCINV')
         call def_eq_dctl_ltcinv(nb,name,field,w,nbfields,nbwvar,parname)
         
      case('MAIS')
         call def_eq_dctl_mais(nb,name,field,w,nbfields,nbwvar,parname)

      case('UVLS')
         call def_eq_dctl_uvls(nb,name,field,w,nbfields,nbwvar,parname)

      case('RT')
         call def_eq_dctl_rt(nb,name,field,w,nbfields,nbwvar,parname)

      case('UVPROT')
         call def_eq_dctl_uvprot(nb,name,field,w,nbfields,nbwvar,parname)

      case('FRT')
         call def_eq_dctl_FRT(nb,name,field,w,nbfields,nbwvar,parname)

      case('SIM_MINMAXVOLT')
         call def_eq_dctl_sim_minmaxvolt(nb,name,field,w,nbfields,nbwvar,parname)

      case('SIM_MINMAXSPEED')
         call def_eq_dctl_sim_minmaxspeed(nb,name,field,w,nbfields,nbwvar,parname)

      case default
         call write_msg_and_stop('def_eq_dctl_model','')
         write(log,"('a DCTL record involves an unknown discrete controller model : ',a20)")modelname
         return

   end select
end subroutine def_eq_dctl_model
!-----------------------------------------------------------------------
!> @brief define injector model observables
!> @details .
subroutine def_obs_dctl_model(nb,modelname,nbobs,obsname)

   implicit none

   integer, intent(out):: nbobs
   character(len=20), intent(in):: modelname
   character(len=10), intent(out):: obsname(*)
   integer, intent(in) :: nb

   select case (modelname)

      case('RT')
         call def_obs_dctl_rt(nb,nbobs,obsname)

      case('FRT')
         call def_obs_dctl_FRT(nb,nbobs,obsname)

      case('VOLT_VAR')
         call def_obs_dctl_volt_var(nb,nbobs,obsname)

      case('SIM_MINMAXVOLT')
         call def_obs_dctl_sim_minmaxvolt(nb,nbobs,obsname)

      case('SIM_MINMAXSPEED')
         call def_obs_dctl_sim_minmaxspeed(nb,nbobs,obsname)

      case default

   end select
end subroutine def_obs_dctl_model
!-----------------------------------------------------------------------
!> @brief  subroutine to initialize discrete controler
!> @details .
!> @param[in] modelname name of the model
!> @param[in] w working variables
subroutine ini_stat_dctl_model(nb,modelname,w)

   implicit none

   double precision, intent(inout):: w(*)
   character(len=20), intent(in):: modelname
   integer, intent(in) :: nb

   select case (modelname)

      case('PST')
         call ini_stat_dctl_pst(nb,w)

      case('LTC')
         call ini_stat_dctl_ltc(nb,w)

      case('LTC2')
         call ini_stat_dctl_ltc2(nb,w)

      case('LTCINV')
         call ini_stat_dctl_ltcinv(nb,w)
         
      case('MAIS')
         call ini_stat_dctl_mais(nb,w)

      case('UVLS')
         call ini_stat_dctl_uvls(nb,w)

      case('RT')
         call ini_stat_dctl_rt(nb,w)

      case('UVPROT')
         call ini_stat_dctl_uvprot(nb,w)

      case('FRT')
         call ini_stat_dctl_FRT(nb,w)

      case('VOLT_VAR')
         call ini_stat_dctl_volt_var(nb,w)

      case('SIM_MINMAXVOLT')
         call ini_stat_dctl_sim_minmaxvolt(nb,w)

      case('SIM_MINMAXSPEED')
         call ini_stat_dctl_sim_minmaxspeed(nb,w)

   end select
end subroutine ini_stat_dctl_model
!-----------------------------------------------------------------------
!> @brief  subroutine to update discrete controler working variables
!> @details .
!> @param[in] modelname name of the model
!> @param[out] w working variables
subroutine upd_w_dctl_model(nb,modelname,w)

   implicit none

   double precision, intent(inout):: w(*)
   character(len=20), intent(in):: modelname
   integer, intent(in) :: nb

   select case (modelname)

      case('PST')
         call upd_w_dctl_pst(nb,w)

      case('LTC')
         call upd_w_dctl_ltc(nb,w)

      case('LTC2')
         call upd_w_dctl_ltc2(nb,w)

      case('LTCINV')
         call upd_w_dctl_ltcinv(nb,w)
         
      case('MAIS')
         call upd_w_dctl_mais(nb,w)

      case('UVLS')
         call upd_w_dctl_uvls(nb,w)

      case('RT')
         call upd_w_dctl_rt(nb,w)

      case('UVPROT')
         call upd_w_dctl_uvprot(nb,w)

      case('FRT')
         call upd_w_dctl_FRT(nb,w)

      case('VOLT_VAR')
         call upd_w_dctl_volt_var(nb,w)

      case('SIM_MINMAXVOLT')
         call upd_w_dctl_sim_minmaxvolt(nb,w)

      case('SIM_MINMAXSPEED')
         call upd_w_dctl_sim_minmaxspeed(nb,w)

   end select
end subroutine upd_w_dctl_model
!-----------------------------------------------------------------------
!> @brief export dctl model observables
!> @details .
subroutine eval_obs_dctl_model(nb,modelname,t,w,obs)

   implicit none

   double precision, intent(in):: t, w(*)
   character(len=20), intent(in):: modelname
   double precision, intent(out):: obs(*)
   integer, intent(in) :: nb

   select case (modelname)

      case('RT')
         call eval_obs_dctl_rt(nb,t,w,obs)

      case('FRT')
         call eval_obs_dctl_FRT(nb,t,w,obs)

      case('VOLT_VAR')
         call eval_obs_dctl_volt_var(nb,t,w,obs)

      case('SIM_MINMAXVOLT')
         call eval_obs_dctl_sim_minmaxvolt(nb,t,w,obs)

      case('SIM_MINMAXSPEED')
         call eval_obs_dctl_sim_minmaxspeed(nb,t,w,obs)

      case default

   end select
end subroutine eval_obs_dctl_model
