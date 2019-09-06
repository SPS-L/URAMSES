!  MODEL NAME : twop_HQSVC              
!  Data :
!       prm(  1)=  G1            ! PSS : first band : input gain
!       prm(  2)=  T1            !                  : time constant
!       prm(  3)=  a             !                  : time constant ratio
!       prm(  4)=  K1            !                  : output gain
!       prm(  5)=  L1            !                  : output limit
!       prm(  6)=  G2            ! PSS : second band : input gain
!       prm(  7)=  T2            !                   : time constant
!       prm(  8)=  b             !                   : time constant ratio
!       prm(  9)=  K2            !                   : output gain
!       prm( 10)=  L2            !                   : output limit
!       prm( 11)=  Ltot          !     : total output limit
!       prm( 12)=  Kp            ! gain of proportional term of PI controller
!       prm( 13)=  Ki            ! gain of integral term of PI controller
!       prm( 14)=  Bp            ! droop
!       prm( 15)=  Bmax          ! maximum shunt susceptance (pu on base Qnom)
!       prm( 16)=  Bmin          ! minimum shunt susceptance (pu on base Qnom)
!       prm( 17)=  Qnom          ! nominal apparent power of compensator (Mvar)
!       prm( 18)=  Intmin        ! lower bound on integral term of PI controller
!       prm( 19)=  Intmax        ! upper bound on itegral term of PI controller
!       prm( 20)=  Tf            ! equivalent time constant of frequency measurement
!  Parameters :
!       prm( 21)=  Bo   initial susceptance
!       prm( 22)=  Vo   voltage set-point
!  Output states :
!       x(  1)=  ix1           real component of current at first bus
!       x(  2)=  iy1           imaginary component of current at first bus
!       x(  3)=  ix2           real component of current at second bus
!       x(  4)=  iy2           imaginary component of current at second bus
!  Internal states defined by user :
!       x(  5)=  B                      shunt suceptance
!       x(  6)=  Bunlim                 output of PI controller
!       x(  7)=  dV                     input of PI controller
!       x(  8)=  deltaVpss              output of PSS
!       x(  9)=  deltaVpssunlim         output of PSS before limiter
!       x( 10)=  dV1                    output of first band
!       x( 11)=  dV1unlim               output of first band before limiter
!       x( 12)=  dV11                   output of first lead-lag of band 1
!       x( 13)=  dV12                   output of second lead-lag of band 1
!       x( 14)=  dV2                    output of second band
!       x( 15)=  dV2unlim               output of first band before limiter
!       x( 16)=  dV21                   output of first lead-lag of band 2
!       x( 17)=  dV22                   output of second lead-lag of band 2
!       x( 18)=  f1                     frequency at bus 1
!       x( 19)=  df1                    frequency deviation at bus 1

!.........................................................................................................

subroutine twop_HQSVC(nb,name,mode,nbxvar,nbzvar,nbdata,nbaddpar,prm,parname,nbobs, &
   obsname,adiy1,adix1,adiy2,adix2,eqtyp,tc,t,omega1,omega2,sbase1,sbase2, &
   bus1,bus2,vx1,vy1,vx2,vy2,ix1,iy1,ix2,iy2,x,z,f,obs)

   use MODELING
   use FREQUENCY
   use ISLAND, only : isl
   use SETTINGS, only : blocktol1,omega_ref,pi
   use FUNCTIONS_IN_MODELS

   implicit none
   double precision, intent(in):: t,vx1,vy1,vx2,vy2,ix1,iy1,ix2,iy2,omega1,omega2,sbase1,sbase2
   double precision, intent(out):: f(*)
   double precision :: obs(*)
   double precision, intent(inout):: x(*),prm(*),tc(*)
   integer, intent(in):: nb,mode,bus1,bus2
   integer, intent(inout):: nbxvar,nbzvar,nbdata,nbaddpar,nbobs,eqtyp(*),z(*),adiy1,adix1,adiy2,adix2
   character(len=20), intent(in):: name
   character(len=10) :: parname(*),obsname(*)

   select case (mode)
   case (define_var_and_par)
      nbdata= 20
      nbaddpar=  2
      parname( 21)='Bo'
      parname( 22)='Vo'
      adix1=  1
      adiy1=  2
      adix2=  3
      adiy2=  4
      nbxvar= 26
      nbzvar=  5

!........................................................................................
   case (define_obs)
      nbobs=  3
      obsname(  1)='B'
      obsname(  2)='deltaVpss'
      obsname(  3)='f1'

!........................................................................................
   case (evaluate_obs)
      obs(  1)=x(  5)              
      obs(  2)=x(  8)              
      obs(  3)=x( 18)              

!........................................................................................
   case (initialize)

!Bo = ([vy2]*[ix2]-[vx2]*[iy2])*sbase2/(([vx2]**2+[vy2]**2)*{Qnom})
      prm( 21)= (vy2*ix2-vx2*iy2)*sbase2/((vx2**2+vy2**2)*prm( 17))

!Vo = dsqrt([vx1]**2+[vy1]**2)+{Bp}*{Bo}
      prm( 22)= dsqrt(vx1**2+vy1**2)+prm( 14)*prm( 21)

!B =  {Bo}
      x(  5)= prm( 21)

!Bunlim =  {Bo}
      x(  6)= prm( 21)

!dV =  0.d0
      x(  7)= 0.d0

!deltaVpss =  0.d0
      x(  8)= 0.d0

!deltaVpssunlim =  0.d0
      x(  9)= 0.d0

!dV1 =  0.d0
      x( 10)= 0.d0

!dV1unlim =  0.d0
      x( 11)= 0.d0

!dV11 =  0.d0
      x( 12)= 0.d0

!dV12 =  0.d0
      x( 13)= 0.d0

!dV2 =  0.d0
      x( 14)= 0.d0

!dV2unlim =  0.d0
      x( 15)= 0.d0

!dV21 =  0.d0
      x( 16)= 0.d0

!dV22 =  0.d0
      x( 17)= 0.d0

!f1 =  1.d0
      x( 18)= 1.d0

!df1 =  0.d0
      x( 19)= 0.d0

!& f_twop_bus1      ! measurement of frequency at bus 1
      x( 20)=vx1
      x( 21)=vy1
      eqtyp(  1)= 20
      eqtyp(  2)= 21
      eqtyp(  3)=0.

!& algeq            ! frequency deviation
      eqtyp(  4)=0

!& tf1p1z           ! PSS : band 1 : first lead-lag
      x( 22)=x( 19)
      eqtyp(  5)= 22
      tc(  5)=prm(  2)
      eqtyp(  6)=0

!& tf1p1z           !                second lead-lag
      x( 23)=x( 19)
      eqtyp(  7)= 23
      tc(  7)=(prm(  2)*prm(  3))
      eqtyp(  8)=0

!& algeq            !                dV1unlim function of dV11 and dV12
      eqtyp(  9)=0

!& lim              !                limiter
      eqtyp( 10)=0
      if(x( 11)>prm(  5))then
         z(  1)=1
      elseif(x( 11)<-prm(  5))then
         z(  1)=-1
      else
         z(  1)=0
      endif

!& tf1p1z           ! PSS : band 2 : first lead-lag
      x( 24)=x( 19)
      eqtyp( 11)= 24
      tc( 11)=prm(  7)
      eqtyp( 12)=0

!& tf1p1z           !                second lead-lag
      x( 25)=x( 19)
      eqtyp( 13)= 25
      tc( 13)=(prm(  7)*prm(  8))
      eqtyp( 14)=0

!& algeq            !                dV2unlim function of dV21 and dV22
      eqtyp( 15)=0

!& lim              !                limiter
      eqtyp( 16)=0
      if(x( 15)>prm( 10))then
         z(  2)=1
      elseif(x( 15)<-prm( 10))then
         z(  2)=-1
      else
         z(  2)=0
      endif

!& algeq            !       sum of two bands
      eqtyp( 17)=0

!& lim              !       final limiter
      eqtyp( 18)=0
      if(x(  9)>prm( 11))then
         z(  3)=1
      elseif(x(  9)<-prm( 11))then
         z(  3)=-1
      else
         z(  3)=0
      endif

!& algeq            ! main summation point
      eqtyp( 19)=0

!& pictllim         ! PI controller
      if(prm( 13)*x(  7)> 0.)then
         z(  4)=1
         eqtyp( 20)=0
         x( 26)=prm( 19)
      elseif(prm( 13)*x(  7)< 0.)then
         z(  4)=-1
         eqtyp( 20)=0
         x( 26)=prm( 18)
      else
         z(  4)=0
         eqtyp( 20)= 26
         x( 26)=x(  6)
      endif
      eqtyp( 21)=0

!& lim              ! final output limiter
      eqtyp( 22)=0
      if(x(  6)>prm( 15))then
         z(  5)=1
      elseif(x(  6)<prm( 16))then
         z(  5)=-1
      else
         z(  5)=0
      endif

!& algeq            ! current injected by susceptance at bus 2 - x comp
      eqtyp( 23)=0

!& algeq            ! current injected by susceptance at bus 2 - y comp
      eqtyp( 24)=0

!& algeq            ! zero current injected at bus 1 - x comp
      eqtyp( 25)=0

!& algeq            ! zero current injected at bus 1 - y comp
      eqtyp( 26)=0

!........................................................................................
   case (evaluate_eqs)

!& f_twop_bus1      ! measurement of frequency at bus 1
      f(  1)=(-x( 20)+vx1)/max(0.05,prm( 20))
      f(  2)=(-x( 21)+vy1)/max(0.05,prm( 20))
      f(  3)=((vy1-x( 21))*x( 20) - (vx1-x( 20))*x( 21))/(2.*pi*fnom*max(0.05,prm( 20))*(x( 20)**2+x( 21)**2))-x( 18)
      if(omega_ref=='COI')then
         f(  3)=f(  3)+omegacoi(isl(bus1),0)
      else
         f(  3)=f(  3)+1.d0
      endif

!& algeq            ! frequency deviation
      f(  4)=x( 18)-1.d0-x( 19)

!& tf1p1z           ! PSS : band 1 : first lead-lag
      f(  5)=-x( 22)+x( 19)
      if (prm(  2)< 0.005)then
         f(  6)=prm(  1)*x( 19)-x( 12)
      else
         f(  6)=prm(  1)*((prm(  2)/prm(  3))*x( 19)+(prm(  2)-(prm(  2)/prm(  3)))*x( 22))-prm(  2)*x( 12)
      endif

!& tf1p1z           !                second lead-lag
      f(  7)=-x( 23)+x( 19)
      if ((prm(  2)*prm(  3))< 0.005)then
         f(  8)=prm(  1)*x( 19)-x( 13)
      else
         f(  8)=prm(  1)*(prm(  2)*x( 19)+((prm(  2)*prm(  3))-prm(  2))*x( 23))-(prm(  2)*prm(  3))*x( 13)
      endif

!& algeq            !                dV1unlim function of dV11 and dV12
      f(  9)=prm(  4)*(x( 12)-x( 13))-x( 11)

!& lim              !                limiter
      select case (z(  1))
         case(0)
            f( 10)=x( 10)-x( 11)
         case(-1)
            f( 10)=x( 10)--prm(  5)
         case(1)
            f( 10)=x( 10)-prm(  5)
      end select

!& tf1p1z           ! PSS : band 2 : first lead-lag
      f( 11)=-x( 24)+x( 19)
      if (prm(  7)< 0.005)then
         f( 12)=prm(  6)*x( 19)-x( 16)
      else
         f( 12)=prm(  6)*((prm(  7)/prm(  8))*x( 19)+(prm(  7)-(prm(  7)/prm(  8)))*x( 24))-prm(  7)*x( 16)
      endif

!& tf1p1z           !                second lead-lag
      f( 13)=-x( 25)+x( 19)
      if ((prm(  7)*prm(  8))< 0.005)then
         f( 14)=prm(  6)*x( 19)-x( 17)
      else
         f( 14)=prm(  6)*(prm(  7)*x( 19)+((prm(  7)*prm(  8))-prm(  7))*x( 25))-(prm(  7)*prm(  8))*x( 17)
      endif

!& algeq            !                dV2unlim function of dV21 and dV22
      f( 15)=prm(  9)*(x( 16)-x( 17))-x( 15)

!& lim              !                limiter
      select case (z(  2))
         case(0)
            f( 16)=x( 14)-x( 15)
         case(-1)
            f( 16)=x( 14)--prm( 10)
         case(1)
            f( 16)=x( 14)-prm( 10)
      end select

!& algeq            !       sum of two bands
      f( 17)=x( 10)+x( 14)-x(  9)

!& lim              !       final limiter
      select case (z(  3))
         case(0)
            f( 18)=x(  8)-x(  9)
         case(-1)
            f( 18)=x(  8)--prm( 11)
         case(1)
            f( 18)=x(  8)-prm( 11)
      end select

!& algeq            ! main summation point
      f( 19)=prm( 22)-dsqrt(vx1**2+vy1**2)+x(  8)-prm( 14)*x(  6)-x(  7)

!& pictllim         ! PI controller
      select case (z(  4))
        case(0)
           f( 20)=prm( 13)*x(  7)
        case(-1)
         f( 20)=x( 26)-prm( 18)
        case(1)
         f( 20)=x( 26)-prm( 19)
      end select
      f( 21)=prm( 12)*x(  7)+x( 26)-x(  6)

!& lim              ! final output limiter
      select case (z(  5))
         case(0)
            f( 22)=x(  5)-x(  6)
         case(-1)
            f( 22)=x(  5)-prm( 16)
         case(1)
            f( 22)=x(  5)-prm( 15)
      end select

!& algeq            ! current injected by susceptance at bus 2 - x comp
      f( 23)=x(  3)-x(  5)*vy2*prm( 17)/sbase2

!& algeq            ! current injected by susceptance at bus 2 - y comp
      f( 24)=x(  4)+x(  5)*vx2*prm( 17)/sbase2

!& algeq            ! zero current injected at bus 1 - x comp
      f( 25)=x(  1)

!& algeq            ! zero current injected at bus 1 - y comp
      f( 26)=x(  2)

!........................................................................................
   case (update_disc)
      select case (z(  1))
         case(0)
            if(x( 11)>prm(  5))then
               z(  1)=1
            elseif(x( 11)<-prm(  5))then
               z(  1)=-1
            endif
         case(-1)
            if(x( 11)>-prm(  5))then
               z(  1)=0
            endif
         case(1)
            if(x( 11)<prm(  5))then
               z(  1)=0
            endif
      end select
      select case (z(  2))
         case(0)
            if(x( 15)>prm( 10))then
               z(  2)=1
            elseif(x( 15)<-prm( 10))then
               z(  2)=-1
            endif
         case(-1)
            if(x( 15)>-prm( 10))then
               z(  2)=0
            endif
         case(1)
            if(x( 15)<prm( 10))then
               z(  2)=0
            endif
      end select
      select case (z(  3))
         case(0)
            if(x(  9)>prm( 11))then
               z(  3)=1
            elseif(x(  9)<-prm( 11))then
               z(  3)=-1
            endif
         case(-1)
            if(x(  9)>-prm( 11))then
               z(  3)=0
            endif
         case(1)
            if(x(  9)<prm( 11))then
               z(  3)=0
            endif
      end select
      select case (z(  4))
         case(0)
            if(x( 26)>prm( 19))then
                  z(  4)=1
                  eqtyp( 20)=0
            elseif(x( 26)<prm( 18))then
                  z(  4)=-1
                  eqtyp( 20)=0
            endif
         case(1)
            if(prm( 13)*x(  7)<0.)then
                  z(  4)=0
                  eqtyp( 20)= 26
            endif
         case(-1)
            if(prm( 13)*x(  7)>0.)then
                  z(  4)=0
                  eqtyp( 20)= 26
            endif
      end select
      select case (z(  5))
         case(0)
            if(x(  6)>prm( 15))then
               z(  5)=1
            elseif(x(  6)<prm( 16))then
               z(  5)=-1
            endif
         case(-1)
            if(x(  6)>prm( 16))then
               z(  5)=0
            endif
         case(1)
            if(x(  6)<prm( 15))then
               z(  5)=0
            endif
      end select
   end select

end subroutine twop_HQSVC
