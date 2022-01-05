module functions_in_models
   
contains
    
   double precision function ppower(vx,vy,ix,iy)

   !   active power

       double precision, intent(in):: vx,vy,ix,iy
    
       ppower=vx*ix+vy*iy
   end function ppower
 
    
   double precision function qpower(vx,vy,ix,iy)

   !   reactive power

       double precision, intent(in):: vx,vy,ix,iy
   
       qpower=vy*ix-vx*iy
   end function qpower
 
    
   double precision function vrectif(if,vin,kc)
    
   !   IEEE standard excitation models
   !   function modelling the voltage drop in rectifiers
    
      double precision, intent(in):: if,vin,kc
      double precision:: in
    
      in=kc*if/max(vin,1d-03)
      if (in <= 0.d0) then
         vrectif=vin
      elseif (in <= 0.433) then
         vrectif=vin-0.577*kc*if
      elseif (in < 0.75) then
         vrectif=dsqrt(0.75*vin**2-(kc*if)**2)
      elseif (in < 1.d0) then
         vrectif=1.732*(vin-kc*if)
      else
         vrectif=0.d0
      endif
   end function vrectif
   
   
double precision function vinrectif(if,vrectif,kc)
    
   !   IEEE standard excitation models
   !   inverse function of vrectif
    
      double precision, intent(in):: if,vrectif,kc
      double precision:: in,vinest
      integer:: nbtries
    
      nbtries=1
      vinest=vrectif+0.577*kc*if
      do 
         in=kc*if/max(vinest,1d-03)
         if (in <= 0.d0) then
            vinrectif=vrectif
         elseif (in <= 0.433) then
            vinrectif=vrectif+0.577*kc*if
         elseif (in < 0.75) then
            vinrectif=dsqrt((vrectif**2+(kc*if)**2)/0.75)  
         else
            vinrectif=(vrectif/1.732)+kc*if
         endif
         if(vinrectif==vinest.or.nbtries>5)exit
         nbtries=nbtries+1
         vinest=vinrectif
      end do
end function vinrectif
 
    
   double precision function vcomp(v,p,q,Kv,Rc,Xc)

   !  compensated voltage magnitude
    
      double precision, intent(in):: v,p,q,Kv,Rc,Xc
    
      vcomp=dsqrt((Kv*v**2+Rc*p+Xc*q)**2+(Xc*p-Rc*q)**2)/v
   end function vcomp
 
    
   double precision function satur(ve,ve1,se1,ve2,se2)

   !   IEEE standard excitation models
   !   saturation function
    
       double precision, intent(in):: ve,ve1,se1,ve2,se2
       double precision:: n
    
       if(ve1==ve2 .or. ve1<=0.d0 .or. ve2<=0.d0 .or. se1<=0.d0 .or. se2<=0.d0)then
          satur=0.d0
       else
          if(ve<=0.d0)then
             satur=0.d0
          else
             n=log10(se1/se2)/log10(ve1/ve2)
             satur=se1*(ve/ve1)**n
          endif
       endif
   end function satur
   
   
   double precision function equal(var1,var2)
   
   !  returns 1.d0 if var1 approaches var2 by less than 1d-6, 0.d0 otherwise
   
      double precision:: var1,var2
      
      if(dabs(var1-var2)<1.d-6)then
          equal=1.d0
      else
          equal=0.d0
      endif
      
   end function equal

   double precision function equalstr(str1,str2)
   
   !  returns 1.d0 if the non-blank part of str1 and str2 are the same, 0.d0 otherwise
   
      character:: str1,str2
      
      if(trim(str1)==trim(str2))then
          equalstr=1.d0
      else
          equalstr=0.d0
      endif
      
   end function equalstr
   
!   double precision function INI_indmach1(name,sbase,SNOM,RS,Lls,LSR,RR,Llr,A,B,LF,vx,vy,ix,iy,typprm)
!
!   use UNITS, only : log
!   use INI_indmach1_mod
!
!   implicit none
!
!   double precision, intent(in):: sbase,RS,Lls,LSR,Llr,RR,A,B,LF,vx,vy,ix,iy
!   double precision, intent(inout):: SNOM
!   character(len=*), intent(in):: typprm
!   character(len=20), intent(in):: name
!
!   double precision f(5),JAC(5,5),fmax,tol,LSS,LRR
!   integer IPVT(5),info,nbit
!
!   select case (typprm)
!
!   case ('')                             ! compute initial values of motor variables
!                                         ! y = ( Bsh, Tmo, psidr0, psiqr0, omega0 )
!      if (SNOM == 0.d0) then
!         if (LF == 0.d0) then
!            write(log,"(' Induction machine ',a20,': has both SNOM and LF equal to zero')")name
!            stop
!         else
!            SNOM=dabs(vy*iy+vx*ix)/LF
!         endif
!      else
!         SNOM=SNOM/sbase
!      endif
!
!      LSS=LSR+Lls
!      LRR=LSR+Llr
!
!      y(1)=0.d0
!      y(2)=(-vy*iy-vx*ix)/SNOM
!      y(3)=-(LSR/LSS)*vx
!      y(4)=(LSR/LSS)*vy
!      if ((vx*ix+vy*iy) < 0)then        ! motor case
!         y(5)=0.99
!      else                              ! generator case
!         y(5)=1.01
!      endif
!
!      fmax=999.
!      nbit=1
!      tol=0.0001d0
!
!      do while (fmax > tol .and. nbit < 11)
!
!         JAC(1,1)=-RS*vx+(LSS-LSR**2/LRR)*vy
!         JAC(1,2)=0.d0
!         JAC(1,3)=0.d0
!         JAC(1,4)=LSR/LRR
!         JAC(1,5)=0.d0
!         JAC(2,1)=RS*vy+(LSS-LSR**2/LRR)*vx
!         JAC(2,2)=0.d0
!         JAC(2,3)=-LSR/LRR
!         JAC(2,4)=0.d0
!         JAC(2,5)=0.d0
!         JAC(3,1)=-vx*LSR*RR/LRR
!         JAC(3,2)=0.d0
!         JAC(3,3)=-RR/LRR
!         JAC(3,4)=-1+y(5)
!         JAC(3,5)=y(4)
!         JAC(4,1)=vy*LSR*RR/LRR
!         JAC(4,2)=0.d0
!         JAC(4,3)=1-y(5)
!         JAC(4,4)=-RR/LRR
!         JAC(4,5)=-y(3)
!         JAC(5,1)=(-y(4)*vx-y(3)*vy)*LSR/LRR
!         JAC(5,2)=-A*y(5)**2-B*y(5)-(1-A-B)
!         JAC(5,3)=(ix-y(1)*vy)*LSR/LRR
!         JAC(5,4)=(-iy-y(1)*vx)*LSR/LRR
!         JAC(5,5)=-y(2)*(2*A*y(5)+B)
!
!         CALL DGETRF (5,5,JAC,5,IPVT,info)
!
!         f(1)=RS*(-iy/SNOM-y(1)*vx) + (LSS-LSR**2/LRR)*(-ix/SNOM+y(1)*vy) + y(4)*LSR/LRR - vy
!         f(2)=RS*(-ix/SNOM+y(1)*vy) - (LSS-LSR**2/LRR)*(-iy/SNOM-y(1)*vx) - y(3)*LSR/LRR - vx
!         f(3)=-y(3)*RR/LRR + (-iy/SNOM-y(1)*vx)*LSR*RR/LRR - (1-y(5))*y(4)
!         f(4)=-y(4)*RR/LRR + (-ix/SNOM+y(1)*vy)*LSR*RR/LRR + (1-y(5))*y(3)
!         f(5)=( y(4)*(-iy/SNOM-y(1)*vx)-y(3)*(-ix/SNOM+y(1)*vy) )*LSR/LRR - y(2)*(A*y(5)**2+B*y(5)+1-A-B)
!         f=-f
!
!         CALL DGETRS('N',5,1,JAC,5,IPVT,f,5,info)
!
!         y=y+f
!
!         fmax=0.d0
!         fmax=maxval(dabs(f(1:5)))
!
!         nbit=nbit+1
!      enddo
!
!      if(fmax > tol)then
!         write(log,"('INDMACH1 ',a20,': iterative initialization of states failed.')")name
!         stop
!      endif
!
!   case('bsh')
!       INI_indmach1=y(1)
!   case('tm0')
!       INI_indmach1=y(2)
!   case('psidr0')
!       INI_indmach1=y(3)
!   case('psiqr0')
!       INI_indmach1=y(4)
!   case('omega0')
!       INI_indmach1=y(5)
!   case default
!       write(log,"('INI_INDMACH1 called with wrong value of typprm')")
!       stop
!   end select
!   
!end function INI_indmach1

double precision function INI_AIR_COND1(name,sbase,SNOM,RS,Lls,LSR,RR,Llr,A,B,LF,alfa,TMNOM,vx,vy,omegaref,ix,iy,typprm)

    use UNITS, only : log

    implicit none

    double precision, intent(in):: sbase,RS,Lls,LSR,Llr,RR,A,B,LF,TMNOM,vx,vy,ix,iy,omegaref
    double precision, intent(inout):: SNOM, alfa
    character(len=*), intent(in):: typprm
    character(len=20), intent(in):: name

    double precision f(5),JAC(5,5),fmax,tol,LSS,LRR
    double precision, save :: y(5)=0.d0
    integer IPVT(5),info,nbit

    select case (typprm)

    case ('')                             ! compute initial values of motor variables
                                         ! y = ( Bsh, Tmo, psidr0, psiqr0, omega0 )
      if (SNOM == 0.d0) then
         if (LF == 0.d0) then
            write(log,"(' Induction machine ',a20,': has both SNOM and LF equal to zero')")name
            stop
         else
            SNOM=dabs(vy*iy+vx*ix)/LF
         endif
      else
         SNOM=SNOM/sbase
      endif

     LSS=LSR+Lls
     LRR=LSR+Llr
     
    if ((vx*ix+vy*iy)<0.d0) then
        alfa=1.d0
    else
        alfa=0.d0
    endif

    if (alfa<0.5d0)then
      y(1)=0.d0 ! bsh0=0
      y(2)=TMNOM ! nominal torque in pu (motor base) used as initial value
      y(3)=0.d0 !psidr0
      y(4)= vy*LRR/(LSR*omegaref) !psiqr0
      y(5)=0.d0 ! omega0
    else


      y(1)=0.d0
      y(2)=(-vy*iy-vx*ix)/SNOM
      y(3)=-(LSR/LSS)*vx
      y(4)=(LSR/LSS)*vy
      if ((vx*ix+vy*iy) < 0)then        ! motor case
         y(5)=0.99
      else                              ! generator case
         y(5)=1.01
      endif

      fmax=999.
      nbit=1
      tol=0.0001d0

      do while (fmax > tol .and. nbit < 11)

         JAC(1,1)=-RS*vx+(LSS-LSR**2/LRR)*vy
         JAC(1,2)=0.d0
         JAC(1,3)=0.d0
         JAC(1,4)=LSR/LRR
         JAC(1,5)=0.d0
         JAC(2,1)=RS*vy+(LSS-LSR**2/LRR)*vx
         JAC(2,2)=0.d0
         JAC(2,3)=-LSR/LRR
         JAC(2,4)=0.d0
         JAC(2,5)=0.d0
         JAC(3,1)=-vx*LSR*RR/LRR
         JAC(3,2)=0.d0
         JAC(3,3)=-RR/LRR
         JAC(3,4)=-1+y(5)
         JAC(3,5)=y(4)
         JAC(4,1)=vy*LSR*RR/LRR
         JAC(4,2)=0.d0
         JAC(4,3)=1-y(5)
         JAC(4,4)=-RR/LRR
         JAC(4,5)=-y(3)
         JAC(5,1)=(-y(4)*vx-y(3)*vy)*LSR/LRR
         JAC(5,2)=-A*y(5)**2-B*y(5)-(1-A-B)
         JAC(5,3)=(ix-y(1)*vy)*LSR/LRR
         JAC(5,4)=(-iy-y(1)*vx)*LSR/LRR
         JAC(5,5)=-y(2)*(2*A*y(5)+B)

         CALL DGETRF (5,5,JAC,5,IPVT,info)

         f(1)=RS*(-iy/SNOM-y(1)*vx) + (LSS-LSR**2/LRR)*(-ix/SNOM+y(1)*vy) + y(4)*LSR/LRR - vy
         f(2)=RS*(-ix/SNOM+y(1)*vy) - (LSS-LSR**2/LRR)*(-iy/SNOM-y(1)*vx) - y(3)*LSR/LRR - vx
         f(3)=-y(3)*RR/LRR + (-iy/SNOM-y(1)*vx)*LSR*RR/LRR - (1-y(5))*y(4)
         f(4)=-y(4)*RR/LRR + (-ix/SNOM+y(1)*vy)*LSR*RR/LRR + (1-y(5))*y(3)
         f(5)=( y(4)*(-iy/SNOM-y(1)*vx)-y(3)*(-ix/SNOM+y(1)*vy) )*LSR/LRR - y(2)*(A*y(5)**2+B*y(5)+1-A-B)
         f=-f

         CALL DGETRS('N',5,1,JAC,5,IPVT,f,5,info)

         y=y+f

         fmax=0.d0
         fmax=maxval(dabs(f(1:5)))

         nbit=nbit+1
      enddo

      if(fmax > tol)then
         write(log,"('INDMACH1 ',a20,': iterative initialization of states failed.')")name
         stop
      endif
    endif
    case('bsh')
       INI_AIR_COND1=y(1)
    case('tm0')
       INI_AIR_COND1=y(2)
    case('psidr0')
       INI_AIR_COND1=y(3)
    case('psiqr0')
       INI_AIR_COND1=y(4)
    case('omega0')
       INI_AIR_COND1=y(5)
    case default
       write(log,"('INI_AIR_COND1 called with wrong value of typprm')")
       stop
    end select

end function INI_AIR_COND1

end module functions_in_models
