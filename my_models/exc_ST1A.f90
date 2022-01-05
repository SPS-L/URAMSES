!  MODEL NAME : exc_ST1A                
!  MODEL FILE : exc_ST1A.txt
!  Compiled: 2017/05/11   8:20
!
!  
!  Data :
!       prm(  1)=  Kv
!       prm(  2)=  Rc
!       prm(  3)=  Xc
!       prm(  4)=  TR
!       prm(  5)=  UEL
!       prm(  6)=  VIMIN
!       prm(  7)=  VIMAX
!       prm(  8)=  VUEL    ERROR - Vuel is the UEL signal not a parameter!
!       prm(  9)=  TC
!       prm( 10)=  TB
!       prm( 11)=  TC1
!       prm( 12)=  TB1
!       prm( 13)=  KA
!       prm( 14)=  TA
!       prm( 15)=  VAMIN
!       prm( 16)=  VAMAX
!       prm( 17)=  VRMIN
!       prm( 18)=  VRMAX
!       prm( 19)=  KC
!       prm( 20)=  KF
!       prm( 21)=  TF
!       prm( 22)=  KLR
!       prm( 23)=  ILR
!  Parameters :
!       prm( 24)=  VREF  
!  Output states :
!       x(  1)=  vf           field voltage
!  Internal states defined by user :
!       x(  2)=  Vc1                   
!       x(  3)=  Vc                    
!       x(  4)=  deltaV                
!       x(  5)=  V1                    
!       x(  6)=  V2                    
!       x(  7)=  V3                    
!       x(  8)=  V4                    
!       x(  9)=  VA                    
!       x( 10)=  VLR                   
!       x( 11)=  VLRlim                
!       x( 12)=  V5                    
!       x( 13)=  V6                    
!       x( 14)=  V7                    
!       x( 15)=  uplim                 
!       x( 16)=  lolim                 
!       x( 17)=  VF                    
!       x( 18)=  VOEL                  

!.........................................................................................................

subroutine exc_ST1A(nb,name,mode,nbxvar,nbzvar,nbdata,nbaddpar,prm,parname,nbobs, &
   obsname,advf,eqtyp,tc,t,v,p,q,omega,if,vf,x,z,f,obs)

#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"exc_ST1A" :: exc_ST1A
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
      nbdata= 23
      nbaddpar=  1
      parname(  1)='Kv'
      parname(  2)='Rc'
      parname(  3)='Xc'
      parname(  4)='TR'
      parname(  5)='UEL'
      parname(  6)='VIMIN'
      parname(  7)='VIMAX'
      parname(  8)='VUEL'
      parname(  9)='TC'
      parname( 10)='TB'
      parname( 11)='TC1'
      parname( 12)='TB1'
      parname( 13)='KA'
      parname( 14)='TA'
      parname( 15)='VAMIN'
      parname( 16)='VAMAX'
      parname( 17)='VRMIN'
      parname( 18)='VRMAX'
      parname( 19)='KC'
      parname( 20)='KF'
      parname( 21)='TF'
      parname( 22)='KLR'
      parname( 23)='ILR'
      parname( 24)='VREF'
      advf=  1
      nbxvar= 21
      nbzvar=  7

!........................................................................................
   case (define_obs)
      nbobs= 15
      obsname(  1)='VA'
      obsname(  2)='Vc1'
      obsname(  3)='Vc'
      obsname(  4)='V1'
      obsname(  5)='V2'
      obsname(  6)='V3'
      obsname(  7)='V4'
      obsname(  8)='vf'
      obsname(  9)='p'
      obsname( 10)='q'
      obsname( 11)='v'
      obsname( 12)='VLRlim'
      obsname( 13)='deltaV'
      obsname( 14)='VREF'
      obsname( 15)='VOEL'

!........................................................................................
   case (evaluate_obs)
      obs(  1)=x(  9)              
      obs(  2)=x(  2)              
      obs(  3)=x(  3)              
      obs(  4)=x(  5)              
      obs(  5)=x(  6)              
      obs(  6)=x(  7)              
      obs(  7)=x(  8)              
      obs(  8)=x(  1)              
      obs(  9)=p                   
      obs( 10)=q                   
      obs( 11)=v                   
      obs( 12)=x( 11)              
      obs( 13)=x(  4)              
      obs( 14)=prm( 24)            
      obs( 15)=x( 18)              

!........................................................................................
   case (initialize)

!VREF     = vcomp([v],[p],[q],{Kv},{Rc},{Xc})+([vf]/{KA})
      prm( 24)= vcomp(v,p,q,prm(  1),prm(  2),prm(  3))+(vf/prm( 13))

!Vc1 =  vcomp([v],[p],[q],{Kv},{Rc},{Xc})
      x(  2)= vcomp(v,p,q,prm(  1),prm(  2),prm(  3))

!Vc =  vcomp([v],[p],[q],{Kv},{Rc},{Xc})
      x(  3)= vcomp(v,p,q,prm(  1),prm(  2),prm(  3))

!deltaV =  ([vf]/{KA})
      x(  4)= (x(  1)/prm( 13))

!V1 =  [vf]/{KA}
      x(  5)= x(  1)/prm( 13)

!V2 =  [vf]/{KA}
      x(  6)= x(  1)/prm( 13)

!V3 =  [vf]/{KA}
      x(  7)= x(  1)/prm( 13)

!V4 =  [vf]/{KA}
      x(  8)= x(  1)/prm( 13)

!VA =  [vf]
      x(  9)= x(  1)

!VLR =  {KLR}*([if]-{ILR})
      x( 10)= prm( 22)*(if-prm( 23))

!VLRlim =  0.
      x( 11)= 0.

!V5 =  [vf]
      x( 12)= x(  1)

!V6 =  [vf]
      x( 13)= x(  1)

!V7 =  [vf]
      x( 14)= x(  1)

!uplim =  [v]*{VRMAX}-{KC}*[if]
      x( 15)= v*prm( 18)-prm( 19)*if

!lolim =  [v]*{VRMIN}
      x( 16)= v*prm( 17)

!VF =  0.
      x( 17)= 0.

!VOEL =  99999.
      x( 18)= 99999.

!& algeq                    fictitious over excitation limit
      eqtyp(  1)=0

!& algeq                    line drop compensation
      eqtyp(  2)=0

!& tf1p                     voltage measurement time constant
      eqtyp(  3)=  3
      tc(  3)=prm(  4)

!& algeq                    summation point of AVR
      eqtyp(  4)=0

!& lim                      limiter on deltaV
      eqtyp(  5)=0
      if(x(  4)-prm(  7)>1.d-6)then
         z(  1)=1
      elseif(x(  4)-prm(  6)<1.d-6)then
         z(  1)=-1
      else
         z(  1)=0
      endif

!& max1v1c                  HV gate for VUEL if UEL=2
      eqtyp(  6)=0
      if(x(  5)-(prm(  8)*equal(prm(  5),2.d0) - 9999.*(1-equal(prm(  5),2.d0)))<1.d-6)then
         z(  2)=1
      else
         z(  2)=2
      endif

!& tf1p1z                   first lead-lag
      x( 19)=x(  6)
      eqtyp(  7)= 19
      tc(  7)=prm( 10)
      eqtyp(  8)=0

!& tf1p1z                   second lead-lag
      x( 20)=x(  7)
      eqtyp(  9)= 20
      tc(  9)=prm( 12)
      eqtyp( 10)=0

!& tf1plim                  amplifier
      if(x(  9)>prm( 16))then
         z(  3)=1
         eqtyp( 11)=0
      elseif(x(  9)<prm( 15))then
         z(  3)=-1
         eqtyp( 11)=0
      else
         z(  3)=0
         eqtyp( 11)=  9
      endif
      tc( 11)=prm( 14)

!& algeq                    immediate field current limiting signal
      eqtyp( 12)=0

!& lim                      ... lower limited to zero
      eqtyp( 13)=0
      if(x( 10)-99999.>1.d-6)then
         z(  4)=1
      elseif(x( 10)-0.<1.d-6)then
         z(  4)=-1
      else
         z(  4)=0
      endif

!& algeq                    ... and subtracted from main signal
      eqtyp( 14)=0

!& max1v1c                  HV gate for VUEL if UEL=3
      eqtyp( 15)=0
      if(x( 12)-(prm(  8)*equal(prm(  5),3.d0) - 9999.*(1-equal(prm(  5),3.d0)))<1.d-6)then
         z(  5)=1
      else
         z(  5)=2
      endif

!& min2v                    LV gate for VOEL
      eqtyp( 16)=0
      if(x( 13)<x( 18))then
         z(  6)=1
      else
         z(  6)=2
      endif

!& algeq                    variable lower limit of final limiter
      eqtyp( 17)=0

!& algeq                    variable upper limit of final limiter
      eqtyp( 18)=0

!& limvb                    final limiter
      eqtyp( 19)=0
      if(x( 14)-x( 15)>1.d-6)then
         z(  7)=1
      elseif(x( 14)-x( 16)<1.d-6)then
         z(  7)=-1
      else
         z(  7)=0
      endif

!& tfder1p                  derivative feedback
      x( 21)=x( 14)
      eqtyp( 20)= 21
      tc( 20)=prm( 21)
      eqtyp( 21)=0

!........................................................................................
   case (evaluate_eqs)

!& algeq                    fictitious over excitation limit
      f(  1)=x( 18)-99999.

!& algeq                    line drop compensation
      f(  2)=x(  2)-vcomp(v,p,q,prm(  1),prm(  2),prm(  3))

!& tf1p                     voltage measurement time constant
      f(  3)=(-x(  3)+1.d0*x(  2))

!& algeq                    summation point of AVR
      f(  4)=x(  4)-prm( 24)+x(  3)-equal(prm(  5),1.d0)*prm(  8)+x( 17)

!& lim                      limiter on deltaV
      select case (z(  1))
         case(0)
            f(  5)=x(  5)-x(  4)
         case(-1)
            f(  5)=x(  5)-prm(  6)
         case(1)
            f(  5)=x(  5)-prm(  7)
      end select

!& max1v1c                  HV gate for VUEL if UEL=2
      select case (z(  2))
         case(1)
            f(  6)=(prm(  8)*equal(prm(  5),2.d0) - 9999.*(1-equal(prm(  5),2.d0)))-x(  6)
         case(2)
            f(  6)=x(  5)-x(  6)
      end select

!& tf1p1z                   first lead-lag
      f(  7)=-x( 19)+x(  6)
      if (prm( 10)< 0.005)then
         f(  8)=1.*x(  6)-x(  7)
      else
         f(  8)=1.*(prm(  9)*x(  6)+(prm( 10)-prm(  9))*x( 19))-prm( 10)*x(  7)
      endif

!& tf1p1z                   second lead-lag
      f(  9)=-x( 20)+x(  7)
      if (prm( 12)< 0.005)then
         f( 10)=1.*x(  7)-x(  8)
      else
         f( 10)=1.*(prm( 11)*x(  7)+(prm( 12)-prm( 11))*x( 20))-prm( 12)*x(  8)
      endif

!& tf1plim                  amplifier
      select case (z(  3))
         case(0)
            f( 11)=-x(  9)+prm( 13)*x(  8)
         case(1)
            f( 11)=x(  9)-prm( 16)
         case(-1)
            f( 11)=x(  9)-prm( 15)
      end select

!& algeq                    immediate field current limiting signal
      f( 12)=x( 10)-prm( 22)*(if-prm( 23))

!& lim                      ... lower limited to zero
      select case (z(  4))
         case(0)
            f( 13)=x( 11)-x( 10)
         case(-1)
            f( 13)=x( 11)-0.
         case(1)
            f( 13)=x( 11)-99999.
      end select

!& algeq                    ... and subtracted from main signal
      f( 14)=x( 12)-x(  9)+x( 11)

!& max1v1c                  HV gate for VUEL if UEL=3
      select case (z(  5))
         case(1)
            f( 15)=(prm(  8)*equal(prm(  5),3.d0) - 9999.*(1-equal(prm(  5),3.d0)))-x( 13)
         case(2)
            f( 15)=x( 12)-x( 13)
      end select

!& min2v                    LV gate for VOEL
      select case (z(  6))
         case(1)
            f( 16)=x( 13)-x( 14)
         case(2)
            f( 16)=x( 18)-x( 14)
      end select

!& algeq                    variable lower limit of final limiter
      f( 17)=v*prm( 17)-x( 16)

!& algeq                    variable upper limit of final limiter
      f( 18)=v*prm( 18)-prm( 19)*if-x( 15)

!& limvb                    final limiter
      select case (z(  7))
         case(0)
            f( 19)=x(  1)-x( 14)
         case(-1)
            f( 19)=x(  1)-x( 16)
         case(1)
            f( 19)=x(  1)-x( 15)
      end select

!& tfder1p                  derivative feedback
      f( 20)=-x( 21)+x( 14)
      if (prm( 21)< 0.005)then
         f( 21)=prm( 20)/prm( 21)*x( 14)-x( 17)
      else
         f( 21)=prm( 20)/prm( 21)*(x( 14)-x( 21))-x( 17)
      endif

!........................................................................................
   case (update_disc)

!& algeq                    fictitious over excitation limit

!& algeq                    line drop compensation

!& tf1p                     voltage measurement time constant

!& algeq                    summation point of AVR

!& lim                      limiter on deltaV
      select case (z(  1))
         case(0)
            if(x(  4)-prm(  7)>1.d-6)then
               z(  1)=1
            elseif(x(  4)-prm(  6)<1.d-6)then
               z(  1)=-1
            endif
         case(-1)
            if(x(  4)-prm(  6)>1.d-6)then
               z(  1)=0
            endif
         case(1)
            if(x(  4)-prm(  7)<1.d-6)then
               z(  1)=0
            endif
      end select

!& max1v1c                  HV gate for VUEL if UEL=2
      select case (z(  2))
         case(1)
            if(x(  5)-(prm(  8)*equal(prm(  5),2.d0) - 9999.*(1-equal(prm(  5),2.d0)))>1.d-6)then
               z(  2)=2
            endif
         case(2)
            if((prm(  8)*equal(prm(  5),2.d0) - 9999.*(1-equal(prm(  5),2.d0)))-x(  5)>1.d-6)then
               z(  2)=1
            endif
      end select

!& tf1p1z                   first lead-lag

!& tf1p1z                   second lead-lag

!& tf1plim                  amplifier
      select case (z(  3))
         case(0)
            if(x(  9)<prm( 15))then
                  z(  3)=-1
                  eqtyp( 11)=0
            elseif(x(  9)>prm( 16))then
                  z(  3)= 1
                  eqtyp( 11)=0
            endif
         case(1)
            if(-x(  9)+prm( 13)*x(  8)<0.)then
                  z(  3)=0
                  eqtyp( 11)=  9
            endif
         case(-1)
            if(-x(  9)+prm( 13)*x(  8)>0.)then
                  z(  3)= 0
                  eqtyp( 11)=  9
            endif
      end select

!& algeq                    immediate field current limiting signal

!& lim                      ... lower limited to zero
      select case (z(  4))
         case(0)
            if(x( 10)-99999.>1.d-6)then
               z(  4)=1
            elseif(x( 10)-0.<1.d-6)then
               z(  4)=-1
            endif
         case(-1)
            if(x( 10)-0.>1.d-6)then
               z(  4)=0
            endif
         case(1)
            if(x( 10)-99999.<1.d-6)then
               z(  4)=0
            endif
      end select

!& algeq                    ... and subtracted from main signal

!& max1v1c                  HV gate for VUEL if UEL=3
      select case (z(  5))
         case(1)
            if(x( 12)-(prm(  8)*equal(prm(  5),3.d0) - 9999.*(1-equal(prm(  5),3.d0)))>1.d-6)then
               z(  5)=2
            endif
         case(2)
            if((prm(  8)*equal(prm(  5),3.d0) - 9999.*(1-equal(prm(  5),3.d0)))-x( 12)>1.d-6)then
               z(  5)=1
            endif
      end select

!& min2v                    LV gate for VOEL
      select case (z(  6))
         case(1)
            if(x( 13)-x( 18)>1.d-6)then
               z(  6)=2
            endif
         case(2)
            if(x( 18)-x( 13)>1.d-6)then
               z(  6)=1
            endif
      end select

!& algeq                    variable lower limit of final limiter

!& algeq                    variable upper limit of final limiter

!& limvb                    final limiter
      select case (z(  7))
         case(0)
            if(x( 14)-x( 15)>1.d-6)then
               z(  7)=1
            elseif(x( 14)-x( 16)<1.d-6)then
               z(  7)=-1
            endif
         case(-1)
            if(x( 14)-x( 16)>1.d-6)then
               z(  7)=0
            endif
         case(1)
            if(x( 14)-x( 15)<1.d-6)then
               z(  7)=0
            endif
      end select

!& tfder1p                  derivative feedback
   end select

end subroutine exc_ST1A
