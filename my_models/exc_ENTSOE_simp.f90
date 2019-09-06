subroutine exc_ENTSOE_simp(nb,name,mode,nbxvar,nbzvar,nbdata,nbaddpar,prm,parname,nbobs, &
   obsname,advf,eqtyp,tc,t,v,p,q,omega,if,vf,x,z,f,obs)

#ifdef DLL
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"exc_ENTSOE_simp" :: exc_ENTSOE_simp
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
      nbdata= 15
      nbaddpar=  1
      parname( 16)='Vo'
      advf=  9
      nbxvar= 14
      nbzvar=  2

   case (define_obs)
      nbobs=  3
      obsname(  1)='domega'
      obsname(  2)='dvpss'
      obsname(  3)='vf'

   case (evaluate_obs)
      obs(  1)=x(  1)
      obs(  2)=x(  6)
      obs(  3)=x(  9)

   case (initialize)
      prm( 16)= v+(vf/prm( 12))
      x(  1)= 0.
      x(  2)= 0.
      x(  3)= 0.
      x(  4)= 0.
      x(  5)= 0.
      x(  6)= 0.
      x(  7)= vf/prm( 12)
      x(  8)= vf/prm( 12)
      x(  9)=vf
      eqtyp(  1)=0
      x( 10)=x(  1)
      eqtyp(  2)= 10
      tc(  2)=prm(  1)
      eqtyp(  3)=0
      x( 11)=x(  2)
      eqtyp(  4)= 11
      tc(  4)=prm(  2)
      eqtyp(  5)=0
      x( 12)=x(  3)
      eqtyp(  6)= 12
      tc(  6)=prm(  5)
      eqtyp(  7)=0
      x( 13)=x(  4)
      eqtyp(  8)= 13
      tc(  8)=prm(  7)
      eqtyp(  9)=0
      eqtyp( 10)=0
      if(x(  5)>prm(  9))then
         z(  1)=1
      elseif(x(  5)<prm(  8))then
         z(  1)=-1
      else
         z(  1)=0
      endif
      eqtyp( 11)=0
      x( 14)=x(  7)
      eqtyp( 12)= 14
      tc( 12)=prm( 11)
      eqtyp( 13)=0
      if(x(  9)>prm( 15))then
         z(  2)=1
         eqtyp( 14)=0
      elseif(x(  9)<prm( 14))then
         z(  2)=-1
         eqtyp( 14)=0
      else
         z(  2)=0
         eqtyp( 14)=  9
      endif
      tc( 14)=prm( 13)

   case (evaluate_eqs)
      f(  1)=omega-1.-x(  1)
      f(  2)=-x( 10)+x(  1)
      if (prm(  1)< 0.005)then
         f(  3)=prm(  1)*x(  1)-x(  2)
      else
         f(  3)=prm(  1)*(x(  1)-x( 10))-prm(  1)*x(  2)
      endif
      f(  4)=-x( 11)+x(  2)
      if (prm(  2)< 0.005)then
         f(  5)=prm(  2)*x(  2)-x(  3)
      else
         f(  5)=prm(  2)*(x(  2)-x( 11))-prm(  2)*x(  3)
      endif
      f(  6)=-x( 12)+x(  3)
      if (prm(  5)< 0.005)then
         f(  7)=prm(  3)*x(  3)-x(  4)
      else
         f(  7)=prm(  3)*(prm(  4)*x(  3)+(prm(  5)-prm(  4))*x( 12))-prm(  5)*x(  4)
      endif
      f(  8)=-x( 13)+x(  4)
      if (prm(  7)< 0.005)then
         f(  9)=1.*x(  4)-x(  5)
      else
         f(  9)=1.*(prm(  6)*x(  4)+(prm(  7)-prm(  6))*x( 13))-prm(  7)*x(  5)
      endif
      select case (z(  1))
         case(0)
            f( 10)=x(  6)-x(  5)
         case(-1)
            f( 10)=x(  6)-prm(  8)
         case(1)
            f( 10)=x(  6)-prm(  9)
      end select
      f( 11)=x(  7)-x(  6)+v-prm( 16)
      f( 12)=-x( 14)+x(  7)
      if (prm( 11)< 0.005)then
         f( 13)=1.*x(  7)-x(  8)
      else
         f( 13)=1.*(prm( 10)*x(  7)+(prm( 11)-prm( 10))*x( 14))-prm( 11)*x(  8)
      endif
      select case (z(  2))
         case(0)
            f( 14)=-x(  9)+prm( 12)*x(  8)
         case(1)
            f( 14)=x(  9)-prm( 15)
         case(-1)
            f( 14)=x(  9)-prm( 14)
      end select

   case (update_disc)
      select case (z(  1))
         case(0)
            if(x(  5)>prm(  9))then
               z(  1)=1
            elseif(x(  5)<prm(  8))then
               z(  1)=-1
            endif
         case(-1)
            if(x(  5)>prm(  8))then
               z(  1)=0
            endif
         case(1)
            if(x(  5)<prm(  9))then
               z(  1)=0
            endif
      end select
      select case (z(  2))
         case(0)
            if(x(  9)<prm( 14))then
                  z(  2)=-1
                  eqtyp( 14)=0
            elseif(x(  9)>prm( 15))then
                  z(  2)= 1
                  eqtyp( 14)=0
            endif
         case(1)
            if(-x(  9)+prm( 12)*x(  8)<0.)then
                  z(  2)=0
                  eqtyp( 14)=  9
            endif
         case(-1)
            if(-x(  9)+prm( 12)*x(  8)>0.)then
                  z(  2)= 0
                  eqtyp( 14)=  9
            endif
      end select
   end select

end subroutine exc_ENTSOE_simp
