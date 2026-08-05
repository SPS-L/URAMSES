module inj_AIR_COND1_mod
    
    double precision :: y(5)=0.d0

end module inj_AIR_COND1_mod

    double precision function INI_AIR_COND1(name,sbase,SNOM,RS,Lls,LSR,RR,Llr,A,B,LF,alfa,TMNOM,vx,vy,omegaref,ix,iy,typprm)

   use SETTINGS, only: write_msg_warning,write_msg_and_stop
   use inj_AIR_COND1_mod

    implicit none

    double precision, intent(in):: sbase,RS,Lls,LSR,Llr,RR,A,B,LF,TMNOM,vx,vy,ix,iy,omegaref
    double precision, intent(inout):: SNOM, alfa
    character(len=*), intent(in):: typprm
    character(len=20), intent(in):: name

    double precision f(5),JAC(5,5),fmax,tol,LSS,LRR
    integer IPVT(5),info,nbit

    select case (typprm)

    case ('')                             ! compute initial values of motor variables
                                         ! y = ( Bsh, Tmo, psidr0, psiqr0, omega0 )
      if (SNOM == 0.d0) then
         if (LF == 0.d0) then
             call write_msg_and_stop('INI_AIR_COND1', ' Induction machine '//trim(name)//': has both SNOM and LF equal to zero')
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
         call write_msg_and_stop('INI_AIR_COND1', ' Induction machine '//trim(name)//': iterative initialization of states failed.')
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
        call write_msg_and_stop('INI_AIR_COND1', ' INI_AIR_COND1 called with wrong value of typprm')
    end select

    end function INI_AIR_COND1

