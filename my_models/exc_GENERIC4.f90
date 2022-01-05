!  MODEL NAME : exc_GENERIC3            
!  MODEL DESCRIPTION FILE : exc_GENERIC3.txt
!  Data :
!       prm(  1)=  G
!       prm(  2)=  Ta
!       prm(  3)=  Tb
!       prm(  4)=  Te
!       prm(  5)=  vfmin
!       prm(  6)=  vfmax
!       prm(  7)=  Kpss
!       prm(  8)=  Tw
!       prm(  9)=  T1
!       prm( 10)=  T2
!       prm( 11)=  T3
!       prm( 12)=  T4
!       prm( 13)=  C
!       prm( 14)=  if1lim
!       prm( 15)=  if2lim
!       prm( 16)=  Toel
!       prm( 17)=  Koel
!       prm( 18)=  L1
!       prm( 19)=  L2
!       prm( 20)=  L3
!  Parameters :
!       prm( 21)=  Vo   AVR voltage set-point
!  Output states :
!       x(  1)=  vf           field voltage
!  Internal states defined by user :
!       x(  2)=  deltaV                 output of main summing junction of AVR
!       x(  3)=  V1                     output of transient gain reduction = input of exciter
!       x(  4)=  x1pss                  PSS : output of washout filter
!       x(  5)=  x2pss                  PSS : output of first lead-lag filter
!       x(  6)=  x3pss                  PSS : output of second lead-lag filter
!       x(  7)=  dvpss                  PSS : output signal added to main summing junction
!       x(  8)=  xoel1                  OEL : if - if1lim
!       x(  9)=  deltaif                OEL : input of integrator
!       x( 10)=  xoel2                  OEL : output of integrator
!       x( 11)=  xoel3                  OEL : output of multiplier
!       x( 12)=  dvoel                  OEL : output signal added to main summing junction

!.........................................................................................................

subroutine exc_GENERIC4(nb,name,mode,nbxvar,nbzvar,nbdata,nbaddpar,prm,parname,nbobs, &
   obsname,advf,eqtyp,tc,t,v,p,q,omega,if,vf,x,z,f,obs)

#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"exc_GENERIC4" :: exc_GENERIC4
#endif

   use MODELING
   use SETTINGS, only : blocktol1
   use FUNCTIONS_IN_MODELS

   implicit none
   double precision, intent(in):: t,v,p,q,omega,if,vf
   double precision, intent(out):: f(*)
   double precision :: obs(*)
   double precision, intent(inout):: x(*),prm(*),tc(*)
   integer, intent(in):: nb,mode
   integer, intent(inout):: nbxvar,nbzvar,nbdata,nbaddpar,advf,nbobs,eqtyp(*),z(*)
   character(len=20), intent(in):: name
   character(len=10) :: parname(*),obsname(*)

   select case (mode)
   case (define_var_and_par)
      nbdata= 20
      nbaddpar=  1
      parname(  1)='G'
      parname(  2)='Ta'
      parname(  3)='Tb'
      parname(  4)='Te'
      parname(  5)='vfmin'
      parname(  6)='vfmax'
      parname(  7)='Kpss'
      parname(  8)='Tw'
      parname(  9)='T1'
      parname( 10)='T2'
      parname( 11)='T3'
      parname( 12)='T4'
      parname( 13)='C'
      parname( 14)='if1lim'
      parname( 15)='if2lim'
      parname( 16)='Toel'
      parname( 17)='Koel'
      parname( 18)='L1'
      parname( 19)='L2'
      parname( 20)='L3'
      parname( 21)='Vo'
      advf=  1
      nbxvar= 16
      nbzvar=  5

!........................................................................................
   case (define_obs)
      nbobs=  4
      obsname(  1)='dvpss'
      obsname(  2)='deltaif'
      obsname(  3)='dvoel'
      obsname(  4)='vf'

!........................................................................................
   case (evaluate_obs)
      obs(  1)=x(  7)              
      obs(  2)=x(  9)              
      obs(  3)=x( 12)              
      obs(  4)=x(  1)              

!........................................................................................
   case (initialize)

!Vo = [v]+[vf]/{G}
      prm( 21)= v+vf/prm(  1)

!deltaV =  [vf]/{G}
      x(  2)= x(  1)/prm(  1)

!V1 =  [vf]
      x(  3)= x(  1)

!x1pss =  0.d0
      x(  4)= 0.d0

!x2pss =  0.d0
      x(  5)= 0.d0

!x3pss =  0.d0
      x(  6)= 0.d0

!dvpss =  0.d0
      x(  7)= 0.d0

!xoel1 =  [if] - {if1lim}
      x(  8)= if - prm( 14)

!deltaif =  [if] - {if1lim}
      x(  9)= if - prm( 14)

!xoel2 =  {L1}
      x( 10)= prm( 18)

!xoel3 =  {L1}*{Koel}
      x( 11)= prm( 18)*prm( 17)

!dvoel =  0.d0
      x( 12)= 0.d0

!& algeq                     main summing junction of AVR
      eqtyp(  1)=0

!& tf1p1z                    transient gain reduction
      x( 13)=x(  2)
      eqtyp(  2)= 13
      tc(  2)=prm(  3)
      eqtyp(  3)=0

!& tf1plim                   exciter
      if(x(  1)>prm(  6))then
         z(  1)=1
         eqtyp(  4)=0
      elseif(x(  1)<prm(  5))then
         z(  1)=-1
         eqtyp(  4)=0
      else
         z(  1)=0
         eqtyp(  4)=  1
      endif
      tc(  4)=prm(  4)

!& tfder1p                   PSS : washout filter
      x( 14)=omega 
      eqtyp(  5)= 14
      tc(  5)=prm(  8)
      eqtyp(  6)=0

!& tf1p1z                    PSS : first lead-lag filter
      x( 15)=x(  4)
      eqtyp(  7)= 15
      tc(  7)=prm( 10)
      eqtyp(  8)=0

!& tf1p1z                    PSS : second lead-lag filter
      x( 16)=x(  5)
      eqtyp(  9)= 16
      tc(  9)=prm( 12)
      eqtyp( 10)=0

!& lim                       PSS : output limiter
      eqtyp( 11)=0
      if(x(  6)>prm( 13))then
         z(  2)=1
      elseif(x(  6)<(-prm( 13)))then
         z(  2)=-1
      else
         z(  2)=0
      endif

!& algeq                    OEL: difference between if and its permanent limit
      eqtyp( 12)=0

!& min1v1c                  upper limited by if2lim
      eqtyp( 13)=0
      if(x(  8)<prm( 15))then
         z(  3)=1
      else
         z(  3)=2
      endif

!& inlim                    OEL integrator (timer)
      if (prm( 16)>= 0.005)then
         tc( 14)=prm( 16)
      endif
      if (x( 10)>prm( 19))then
         z(  4)=1
         eqtyp( 14)=0
      elseif (x( 10)<prm( 18)) then
         z(  4)=-1
         eqtyp( 14)=0
      else
         z(  4)=0
         if (prm( 16)>= 0.005)then
            eqtyp( 14)= 10
         else
            eqtyp( 14)=0
         endif
      endif

!& algeq                    multiplication by KOEL
      eqtyp( 15)=0

!& lim                      output limiter of OEL
      eqtyp( 16)=0
      if(x( 11)>prm( 20))then
         z(  5)=1
      elseif(x( 11)<0.d0)then
         z(  5)=-1
      else
         z(  5)=0
      endif

!........................................................................................
   case (evaluate_eqs)

!& algeq                     main summing junction of AVR
      f(  1)=prm( 21)-v+x(  7)-x( 12)-x(  2)

!& tf1p1z                    transient gain reduction
      f(  2)=-x( 13)+x(  2)
      if (prm(  3)< 0.005)then
         f(  3)=prm(  1)*x(  2)-x(  3)
      else
         f(  3)=prm(  1)*(prm(  2)*x(  2)+(prm(  3)-prm(  2))*x( 13))-prm(  3)*x(  3)
      endif

!& tf1plim                   exciter
      select case (z(  1))
         case(0)
            f(  4)=-x(  1)+1.d0*x(  3)
         case(1)
            f(  4)=x(  1)-prm(  6)
         case(-1)
            f(  4)=x(  1)-prm(  5)
      end select

!& tfder1p                   PSS : washout filter
      f(  5)=-x( 14)+omega 
      if (prm(  8)< 0.005)then
         f(  6)=(prm(  7)/prm(  8))*omega -x(  4)
      else
         f(  6)=(prm(  7)/prm(  8))*(omega -x( 14))-x(  4)
      endif

!& tf1p1z                    PSS : first lead-lag filter
      f(  7)=-x( 15)+x(  4)
      if (prm( 10)< 0.005)then
         f(  8)=1.d0*x(  4)-x(  5)
      else
         f(  8)=1.d0*(prm(  9)*x(  4)+(prm( 10)-prm(  9))*x( 15))-prm( 10)*x(  5)
      endif

!& tf1p1z                    PSS : second lead-lag filter
      f(  9)=-x( 16)+x(  5)
      if (prm( 12)< 0.005)then
         f( 10)=1.d0*x(  5)-x(  6)
      else
         f( 10)=1.d0*(prm( 11)*x(  5)+(prm( 12)-prm( 11))*x( 16))-prm( 12)*x(  6)
      endif

!& lim                       PSS : output limiter
      select case (z(  2))
         case(0)
            f( 11)=x(  7)-x(  6)
         case(-1)
            f( 11)=x(  7)-(-prm( 13))
         case(1)
            f( 11)=x(  7)-prm( 13)
      end select

!& algeq                    OEL: difference between if and its permanent limit
      f( 12)=if-prm( 14)-x(  8)

!& min1v1c                  upper limited by if2lim
      select case (z(  3))
         case(1)
            f( 13)=x(  8)-x(  9)
         case(2)
            f( 13)=prm( 15)-x(  9)
      end select

!& inlim                    OEL integrator (timer)
      if (prm( 16)>= 0.005)then
         select case (z(  4))
            case(0)
               f( 14)=x(  9)
            case(1)
               f( 14)=x( 10)-prm( 19)
            case(-1)
               f( 14)=x( 10)-prm( 18)
         end select
      else
         select case (z(  4))
            case(0)
               f( 14)=x(  9)-x( 10)
            case(1)
               f( 14)=x( 10)-prm( 19)
            case(-1)
               f( 14)=x( 10)-prm( 18)
         end select
      endif

!& algeq                    multiplication by KOEL
      f( 15)=prm( 17)*x( 10)-x( 11)

!& lim                      output limiter of OEL
      select case (z(  5))
         case(0)
            f( 16)=x( 12)-x( 11)
         case(-1)
            f( 16)=x( 12)-0.d0
         case(1)
            f( 16)=x( 12)-prm( 20)
      end select

!........................................................................................
   case (update_disc)

!& algeq                     main summing junction of AVR

!& tf1p1z                    transient gain reduction

!& tf1plim                   exciter
      select case (z(  1))
         case(0)
            if(x(  1)<prm(  5))then
                  z(  1)=-1
                  eqtyp(  4)=0
            elseif(x(  1)>prm(  6))then
                  z(  1)= 1
                  eqtyp(  4)=0
            endif
         case(1)
            if(-x(  1)+1.d0*x(  3)<0.)then
                  z(  1)=0
                  eqtyp(  4)=  1
            endif
         case(-1)
            if(-x(  1)+1.d0*x(  3)>0.)then
                  z(  1)= 0
                  eqtyp(  4)=  1
            endif
      end select

!& tfder1p                   PSS : washout filter

!& tf1p1z                    PSS : first lead-lag filter

!& tf1p1z                    PSS : second lead-lag filter

!& lim                       PSS : output limiter
      select case (z(  2))
         case(0)
            if(x(  6)>prm( 13))then
               z(  2)=1
            elseif(x(  6)<(-prm( 13)))then
               z(  2)=-1
            endif
         case(-1)
            if(x(  6)>(-prm( 13)))then
               z(  2)=0
            endif
         case(1)
            if(x(  6)<prm( 13))then
               z(  2)=0
            endif
      end select

!& algeq                    OEL: difference between if and its permanent limit

!& min1v1c                  upper limited by if2lim
      select case (z(  3))
         case(1)
            if(x(  8)>prm( 15))then
               z(  3)=2
            endif
         case(2)
            if(prm( 15)>x(  8))then
               z(  3)=1
            endif
      end select

!& inlim                    OEL integrator (timer)
      if (prm( 16)>= 0.005)then
         select case (z(  4))
            case(0)
               if(x( 10)<prm( 18))then
                  z(  4)=-1
                  eqtyp( 14)=0
               elseif(x( 10)>prm( 19))then
                  z(  4)= 1
                  eqtyp( 14)=0
               endif
            case(1)
               if(x(  9)<0.)then
                  z(  4)=0
                  eqtyp( 14)= 10
               endif
            case(-1)
               if(x(  9)>0.)then
                  z(  4)=0
                  eqtyp( 14)= 10
               endif
         end select
      else
         select case (z(  4))
            case(0)
               if(x( 10)<prm( 18))then
                  z(  4)=-1
               elseif(x( 10)>prm( 19))then
                  z(  4)= 1
               endif
            case(1)
               if(x(  9)<prm( 19))then
                  z(  4)=0
               endif
            case(-1)
               if(x(  9)>prm( 18))then
                  z(  4)=0
               endif
         end select
      endif

!& algeq                    multiplication by KOEL

!& lim                      output limiter of OEL
      select case (z(  5))
         case(0)
            if(x( 11)>prm( 20))then
               z(  5)=1
            elseif(x( 11)<0.d0)then
               z(  5)=-1
            endif
         case(-1)
            if(x( 11)>0.d0)then
               z(  5)=0
            endif
         case(1)
            if(x( 11)<prm( 20))then
               z(  5)=0
            endif
      end select
   end select

end subroutine exc_GENERIC4
