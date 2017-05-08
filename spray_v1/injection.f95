!
!=====INJECTION INLET PARAMETERS
!
      subroutine injection
      include "common.inc"
!
      i=iinj
      j=jinj
      k=kinj
!
!-----CONE HALF-ANGLE
!
!       beta=atan(0.7*(dng(i,j,k)/denl(i,j,k))**0.5)
!
!-----SIZE OF INJECTOR CELLS
!
      rs=rinj+xinj*tan(beta*(pi/180.))
      dinj=rs/nradic
!
!-----VELOCITIES AT INLET
!
      uinj=cd*(2*(pinj-p(i,j,k))&
           /dnl(i,j,k))**0.5
      uqin=uinj&
           *((-(1-dng(i,j,k)/dnl(i,j,k))&
           +((1-dng(i,j,k)/dnl(i,j,k))**2&
           +4*dng(i,j,k)*rs**2&
           /(dnl(i,j,k)*rinj**2))**0.5)&
           /(2*dng(i,j,k)*rs**2&
           /(dnl(i,j,k)*rinj**2)))
!
      vqin=uqin*tan(beta*(pi/180.))
      wqin=uqin*tan(beta*(pi/180.))
      ruqin=vqin
!
      ugin=uqin*0.95
      vgin=vqin*0.95
      wgin=wqin*0.95
      rugin=vgin
!
!-----LIQUID VOLUME FRACTION AT INLET
!
      void=1-(uinj/uqin)*(rinj**2/rs**2)
      q3(i-1,j,k)=(1-void)/cnt
!
!-----MOMENTS AT INLET BASED ON GAMMA DISTRIBUTION
!
      q2(i-1,j,k)=q3(i-1,j,k)/r32in
      q1(i-1,j,k)=q2(i-1,j,k)*(pkin+2)/(r32in*(pkin+1))
      q0(i-1,j,k)=q1(i-1,j,k)*(pkin+2)/(r32in*pkin)
      q4(i-1,j,k)=((pkin+3)/(pkin+2))*r32in*q3(i-1,j,k)
      q5(i-1,j,k)=((pkin+4)/(pkin+2))*r32in*q4(i-1,j,k)
      q6(i-1,j,k)=((pkin+5)/(pkin+2))*r32in*q5(i-1,j,k)
      q7(i-1,j,k)=((pkin+6)/(pkin+2))*r32in*q6(i-1,j,k)
!
!-----INLET MOMENTS BASED ON LAGUERRE DISTRIBUTION
!     USING GAMMA DISTRIBUTION MOMENTS
!
      q0in=0;q1in=0;q2in=0;q3in=0
      q4in=0;q5in=0;q6in=0;q7in=0
!
      do iqq=0,iqmax-1
        power=iqq
        rlb=0.;rub=0.
        call moment(qa,power,rlb,rub,i-1,j,k)
        if(iqq.eq.0)q0in=qa
        if(iqq.eq.1)q1in=qa
        if(iqq.eq.2)q2in=qa
        if(iqq.eq.3)q3in=qa
        if(iqq.eq.4)q4in=qa
        if(iqq.eq.5)q5in=qa
        if(iqq.eq.6)q6in=qa
        if(iqq.eq.7)q7in=qa
      enddo
      qcheck=q3in/q3(i-1,j,k)
!
!-----NORMALISATION INLET FLUXES FOR LIQUID AND GAS PHASE
!
      fluxn1=pi*rs**2*ugin
      fluxn2=pi*rs**2*uqin
!
!-----NUMBER OF TIME STEPS
!
      ntsmx=tmax/delt
!
!-----CREATE GRID
!
      call grid
!
!-----ININTIALISATION DATA
!
      print*,'=='
      print*,'INITIALISATION PARAMETERS'
      print*,'=='
      print*,'ia, ja, ka'
      print*,ia,ja,ka
      print*,'q0in, q1in, q2in'
      print*,q0in,q1in,q2in
      print*,'q3in, q4in, q5in'
      print*,q3in,q4in,q5in
      print*,'q6in, q7in, r32in'
      print*,q6in,q7in,(q3in/q2in)
      print*,'rlb, rub, Q check [-> 1]'
      print*,rlb,rub,qcheck
      print*,'beta, uqin, delP'
      print*,beta,uqin,(pinj-p(i,j,k))
      print*,'delt, tmax, ntsmx'
      print*,delt,tmax,ntsmx
      print*,'=='
!
      return
      end