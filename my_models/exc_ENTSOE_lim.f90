!  MODEL NAME : exc_ENTSOE_lim          
!  Data :
!       prm(  1)=  TW1
!       prm(  2)=  TW2
!       prm(  3)=  KS1
!       prm(  4)=  T1
!       prm(  5)=  T2
!       prm(  6)=  T3
!       prm(  7)=  T4
!       prm(  8)=  VSTMIN
!       prm(  9)=  VSTMAX
!       prm( 10)=  TA
!       prm( 11)=  TB
!       prm( 12)=  KE
!       prm( 13)=  TE
!       prm( 14)=  EMIN
!       prm( 15)=  EMAX
!       prm( 16)=  IFDN
!       prm( 17)=  TOEL
!       prm( 18)=  LOEL
!       prm( 19)=  UOEL
!       prm( 20)=  KOEL
!       prm( 21)=  OELLI
!  Parameters :
!       prm( 22)=  Vo  
!       prm( 23)=  ad_deltavoel  
!  Output states :
!       x(  1)=  vf           field voltage
!  Internal states defined by user :
!       x(  2)=  domega                
!       x(  3)=  pss1                  
!       x(  4)=  pss2                  
!       x(  5)=  pss3                  
!       x(  6)=  pss4                  
!       x(  7)=  dvpss                 
!       x(  8)=  avr1                  
!       x(  9)=  avr2                  
!       x( 10)=  deltaif1              
!       x( 11)=  deltaif2              
!       x( 12)=  intoel1               
!       x( 13)=  intoel2               
!       x( 14)=  deltavoel             

!.........................................................................................................

subroutine exc_ENTSOE_lim(nb,name,mode,nbxvar,nbzvar,nbdata,nbaddpar,prm,parname,nbobs, &
   obsname,advf,eqtyp,tc,t,v,p,q,omega,if,vf,x,z,f,obs)

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
      nbdata= 21
      nbaddpar=  2
      parname( 22)='Vo'
      parname( 23)='ad_deltavoel'
      advf=  1
      nbxvar= 19
      nbzvar=  5

!........................................................................................
   case (define_obs)
      nbobs=  3
      obsname(  1)='domega'
      obsname(  2)='dvpss'
      obsname(  3)='vf'

!........................................................................................
   case (evaluate_obs)
      obs(  1)=x(  2)              
      obs(  2)=x(  7)              
      obs(  3)=x(  1)              

!........................................................................................
   case (initialize)

!Vo  = v+(vf/{KE})
      prm( 22)= v+(vf/prm( 12))

!ad_deltavoel = 14
      prm( 23)= 14

!domega =  0.
      x(  2)= 0.

!pss1 =  0.
      x(  3)= 0.

!pss2 =  0.
      x(  4)= 0.

!pss3 =  0.
      x(  5)= 0.

!pss4 =  0.
      x(  6)= 0.

!dvpss =  0.
      x(  7)= 0.

!avr1 =  vf/{KE}
      x(  8)= vf/prm( 12)

!avr2 =  vf/{KE}
      x(  9)= vf/prm( 12)

!deltaif1 =  [if]-1.05*{IFDN}
      x( 10)= if-1.05*prm( 16)

!deltaif2 =  [if]-1.05*{IFDN}
      x( 11)= if-1.05*prm( 16)

!intoel1 =  {LOEL}
      x( 12)= prm( 18)

!intoel2 =  {KOEL}*{LOEL}
      x( 13)= prm( 20)*prm( 18)

!deltavoel =  0.d0
      x( 14)= 0.d0

!& algeq           PSS
      eqtyp(  1)=0

!& tfder1p
      x( 15)=x(  2)
      eqtyp(  2)= 15
      tc(  2)=prm(  1)
      eqtyp(  3)=0

!& tfder1p
      x( 16)=x(  3)
      eqtyp(  4)= 16
      tc(  4)=prm(  2)
      eqtyp(  5)=0

!& tf1p1z
      x( 17)=x(  4)
      eqtyp(  6)= 17
      tc(  6)=prm(  5)
      eqtyp(  7)=0

!& tf1p1z
      x( 18)=x(  5)
      eqtyp(  8)= 18
      tc(  8)=prm(  7)
      eqtyp(  9)=0

!& lim
      eqtyp( 10)=0
      if(x(  6)>prm(  9))then
         z(  1)=1
      elseif(x(  6)<prm(  8))then
         z(  1)=-1
      else
         z(  1)=0
      endif

!& algeq             main AVR summing junction
      eqtyp( 11)=0

!& tf1p1z            AVR : lead-lag filter
      x( 19)=x(  8)
      eqtyp( 12)= 19
      tc( 12)=prm( 11)
      eqtyp( 13)=0

!& tf1plim           exciter time constant with limits
      if(x(  1)>prm( 15))then
         z(  2)=1
         eqtyp( 14)=0
      elseif(x(  1)<prm( 14))then
         z(  2)=-1
         eqtyp( 14)=0
      else
         z(  2)=0
         eqtyp( 14)=  1
      endif
      tc( 14)=prm( 13)

!& algeq             OEL: difference between if and its permanent limit
      eqtyp( 15)=0

!& min1v1c                upper limited by 0.35*IFN
      eqtyp( 16)=0
      if(x( 10)<(0.35*prm( 16)))then
         z(  3)=1
      else
         z(  3)=2
      endif

!& inlim                  OEL integrator (timer)
      if (prm( 17)>= 0.005)then
         tc( 17)=prm( 17)
      endif
      if (x( 12)>prm( 19))then
         z(  4)=1
         eqtyp( 17)=0
      elseif (x( 12)<prm( 18)) then
         z(  4)=-1
         eqtyp( 17)=0
      else
         z(  4)=0
         if (prm( 17)>= 0.005)then
            eqtyp( 17)= 12
         else
            eqtyp( 17)=0
         endif
      endif

!& algeq                  multiplication by KOEL
      eqtyp( 18)=0

!& lim                    output limiter of OEL
      eqtyp( 19)=0
      if(x( 13)>0.)then
         z(  5)=1
      elseif(x( 13)<prm( 21))then
         z(  5)=-1
      else
         z(  5)=0
      endif

!........................................................................................
   case (evaluate_eqs)

!& algeq           PSS
      f(  1)=omega-1.-x(  2)

!& tfder1p
      f(  2)=-x( 15)+x(  2)
      if (prm(  1)< 0.005)then
         f(  3)=prm(  1)*x(  2)-x(  3)
      else
         f(  3)=prm(  1)*(x(  2)-x( 15))-x(  3)
      endif

!& tfder1p
      f(  4)=-x( 16)+x(  3)
      if (prm(  2)< 0.005)then
         f(  5)=prm(  2)*x(  3)-x(  4)
      else
         f(  5)=prm(  2)*(x(  3)-x( 16))-x(  4)
      endif

!& tf1p1z
      f(  6)=-x( 17)+x(  4)
      if (prm(  5)< 0.005)then
         f(  7)=prm(  3)*x(  4)-x(  5)
      else
         f(  7)=prm(  3)*(prm(  4)*x(  4)+(prm(  5)-prm(  4))*x( 17))-prm(  5)*x(  5)
      endif

!& tf1p1z
      f(  8)=-x( 18)+x(  5)
      if (prm(  7)< 0.005)then
         f(  9)=1.*x(  5)-x(  6)
      else
         f(  9)=1.*(prm(  6)*x(  5)+(prm(  7)-prm(  6))*x( 18))-prm(  7)*x(  6)
      endif

!& lim
      select case (z(  1))
         case(0)
            f( 10)=x(  7)-x(  6)
         case(-1)
            f( 10)=x(  7)-prm(  8)
         case(1)
            f( 10)=x(  7)-prm(  9)
      end select

!& algeq             main AVR summing junction
      f( 11)=x(  8)-x(  7)-x( 14)+v-prm( 22)

!& tf1p1z            AVR : lead-lag filter
      f( 12)=-x( 19)+x(  8)
      if (prm( 11)< 0.005)then
         f( 13)=1.*x(  8)-x(  9)
      else
         f( 13)=1.*(prm( 10)*x(  8)+(prm( 11)-prm( 10))*x( 19))-prm( 11)*x(  9)
      endif

!& tf1plim           exciter time constant with limits
      select case (z(  2))
         case(0)
            f( 14)=-x(  1)+prm( 12)*x(  9)
         case(1)
            f( 14)=x(  1)-prm( 15)
         case(-1)
            f( 14)=x(  1)-prm( 14)
      end select

!& algeq             OEL: difference between if and its permanent limit
      f( 15)=x( 10)-if+1.05*prm( 16)

!& min1v1c                upper limited by 0.35*IFN
      select case (z(  3))
         case(1)
            f( 16)=x( 10)-x( 11)
         case(2)
            f( 16)=(0.35*prm( 16))-x( 11)
      end select

!& inlim                  OEL integrator (timer)
      if (prm( 17)>= 0.005)then
         select case (z(  4))
            case(0)
               f( 17)=x( 11)
            case(1)
               f( 17)=x( 12)-prm( 19)
            case(-1)
               f( 17)=x( 12)-prm( 18)
         end select
      else
         select case (z(  4))
            case(0)
               f( 17)=x( 11)-x( 12)
            case(1)
               f( 17)=x( 12)-prm( 19)
            case(-1)
               f( 17)=x( 12)-prm( 18)
         end select
      endif

!& algeq                  multiplication by KOEL
      f( 18)=x( 13)-prm( 20)*x( 12)

!& lim                    output limiter of OEL
      select case (z(  5))
         case(0)
            f( 19)=x( 14)-x( 13)
         case(-1)
            f( 19)=x( 14)-prm( 21)
         case(1)
            f( 19)=x( 14)-0.
      end select

!........................................................................................
   case (update_disc)
      select case (z(  1))
         case(0)
            if(x(  6)>prm(  9))then
               z(  1)=1
            elseif(x(  6)<prm(  8))then
               z(  1)=-1
            endif
         case(-1)
            if(x(  6)>prm(  8))then
               z(  1)=0
            endif
         case(1)
            if(x(  6)<prm(  9))then
               z(  1)=0
            endif
      end select
      select case (z(  2))
         case(0)
            if(x(  1)<prm( 14))then
                  z(  2)=-1
                  eqtyp( 14)=0
            elseif(x(  1)>prm( 15))then
                  z(  2)= 1
                  eqtyp( 14)=0
            endif
         case(1)
            if(-x(  1)+prm( 12)*x(  9)<0.)then
                  z(  2)=0
                  eqtyp( 14)=  1
            endif
         case(-1)
            if(-x(  1)+prm( 12)*x(  9)>0.)then
                  z(  2)= 0
                  eqtyp( 14)=  1
            endif
      end select
      select case (z(  3))
         case(1)
            if(x( 10)>(0.35*prm( 16)))then
               z(  3)=2
            endif
         case(2)
            if((0.35*prm( 16))>x( 10))then
               z(  3)=1
            endif
      end select
      if (prm( 17)>= 0.005)then
         select case (z(  4))
            case(0)
               if(x( 12)<prm( 18))then
                  z(  4)=-1
                  eqtyp( 17)=0
               elseif(x( 12)>prm( 19))then
                  z(  4)= 1
                  eqtyp( 17)=0
               endif
            case(1)
               if(x( 11)<0.)then
                  z(  4)=0
                  eqtyp( 17)= 12
               endif
            case(-1)
               if(x( 11)>0.)then
                  z(  4)=0
                  eqtyp( 17)= 12
               endif
         end select
      else
         select case (z(  4))
            case(0)
               if(x( 12)<prm( 18))then
                  z(  4)=-1
               elseif(x( 12)>prm( 19))then
                  z(  4)= 1
               endif
            case(1)
               if(x( 11)<prm( 19))then
                  z(  4)=0
               endif
            case(-1)
               if(x( 11)>prm( 18))then
                  z(  4)=0
               endif
         end select
      endif
      select case (z(  5))
         case(0)
            if(x( 13)>0.)then
               z(  5)=1
            elseif(x( 13)<prm( 21))then
               z(  5)=-1
            endif
         case(-1)
            if(x( 13)>prm( 21))then
               z(  5)=0
            endif
         case(1)
            if(x( 13)<0.)then
               z(  5)=0
            endif
      end select
   end select

end subroutine exc_ENTSOE_lim
