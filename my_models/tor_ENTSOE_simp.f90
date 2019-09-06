!  MODEL NAME : tor_ENTSOE_simp
!  Data :
!       prm(  1)=  R
!       prm(  2)=  T1
!       prm(  3)=  VMIN
!       prm(  4)=  VMAX
!       prm(  5)=  T2
!       prm(  6)=  T3
!  Parameters :
!       prm(  7)=  C
!  Output states :
!       x(  1)=  tm           mechanical torque
!  Internal states defined by user :
!       x(  2)=  dp1
!       x(  3)=  dp2
!       x(  4)=  Pm

!.........................................................................................................

subroutine tor_ENTSOE_simp(nb,name,mode,nbxvar,nbzvar,nbdata,nbaddpar,prm,parname,nbobs, &
   obsname,adtm,eqtyp,tc,t,p,tm,omega,x,z,f,obs)

#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"tor_ENTSOE_simp" :: tor_ENTSOE_simp
#endif

   use MODELING
   use SETTINGS, only : blocktol1
   use FUNCTIONS_IN_MODELS

   implicit none
   double precision, intent(in):: t,p,tm,omega
   double precision, intent(out):: f(*)
   double precision :: obs(*)
   double precision, intent(inout):: x(*),prm(*),tc(*)
   integer, intent(in):: nb,mode
   integer, intent(inout):: nbxvar,nbzvar,nbdata,nbaddpar,adtm,nbobs,eqtyp(*),z(*)
   character(len=20), intent(in):: name
   character(len=10) :: parname(*),obsname(*)

   select case (mode)
   case (define_var_and_par)
      nbdata=  6
      nbaddpar=  1
      parname(  7)='Tm0'
      adtm=  1
      nbxvar=  5
      nbzvar=  1

!........................................................................................
   case (define_obs)
      nbobs=  1
      obsname(  1)='Pm'

!........................................................................................
   case (evaluate_obs)
      obs(  1)=x(  4)

!........................................................................................
   case (initialize)

!C   =  [tm]*{R}
      prm(  7)=  tm*prm(  1)

!dp1 =  [tm]
      x(  2)= x(  1)

!dp2 =  [tm]
      x(  3)= x(  1)

!Pm =  [tm]
      x(  4)= x(  1)

!& algeq
      eqtyp(  1)=0

!& tf1plim
      if(x(  3)>prm(  4))then
         z(  1)=1
         eqtyp(  2)=0
      elseif(x(  3)<prm(  3))then
         z(  1)=-1
         eqtyp(  2)=0
      else
         z(  1)=0
         eqtyp(  2)=  3
      endif
      tc(  2)=prm(  2)

!& tf1p1z
      x(  5)=x(  3)
      eqtyp(  3)=  5
      tc(  3)=prm(  6)
      eqtyp(  4)=0

!& algeq
      eqtyp(  5)=0

!........................................................................................
   case (evaluate_eqs)

!& algeq
      f(  1)=prm(  1)*x(  2)-prm(  7)+omega-1.d0

!& tf1plim
      select case (z(  1))
         case(0)
            f(  2)=-x(  3)+1.*x(  2)
         case(1)
            f(  2)=x(  3)-prm(  4)
         case(-1)
            f(  2)=x(  3)-prm(  3)
      end select

!& tf1p1z
      f(  3)=-x(  5)+x(  3)
      if (prm(  6)< 0.005)then
         f(  4)=1.*x(  3)-x(  4)
      else
         f(  4)=1.*(prm(  5)*x(  3)+(prm(  6)-prm(  5))*x(  5))-prm(  6)*x(  4)
      endif

!& algeq
      f(  5)=x(  1)*omega-x(  4)

!........................................................................................
   case (update_disc)
      select case (z(  1))
         case(0)
            if(x(  3)<prm(  3))then
                  z(  1)=-1
                  eqtyp(  2)=0
            elseif(x(  3)>prm(  4))then
                  z(  1)= 1
                  eqtyp(  2)=0
            endif
         case(1)
            if(-x(  3)+1.*x(  2)<0.)then
                  z(  1)=0
                  eqtyp(  2)=  3
            endif
         case(-1)
            if(-x(  3)+1.*x(  2)>0.)then
                  z(  1)= 0
                  eqtyp(  2)=  3
            endif
      end select
   end select

end subroutine tor_ENTSOE_simp
