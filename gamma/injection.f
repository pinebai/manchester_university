!
!=====INJECTION INLET PARAMETERS
!
      subroutine injection
      include "common.inc"
!
!-----SIZE OF INJECTOR CELLS
!
      rs=rinj+xinj*tan(beta*(pi/180.))
      dinj=rs/nradic
!
!-----VELOCITIES AT INLET
!
      i=iinj
      j=jinj
      k=kinj
      uinj=cd*(2*(pinj-p(i,j,k))
     &     /dnl(i,j,k))**0.5
      ulin=uinj
     &     *((-(1-dng(i,j,k)/dnl(i,j,k))
     &     +((1-dng(i,j,k)/dnl(i,j,k))**2
     &     +4*dng(i,j,k)*rs**2
     &     /(dnl(i,j,k)*rinj**2))**0.5)
     &     /(2*dng(i,j,k)*rs**2
     &     /(dnl(i,j,k)*rinj**2)))
!
      vlin=ulin*tan(beta*(pi/180.))
      wlin=ulin*tan(beta*(pi/180.))
!
      ugin=ulin*0.9
      vgin=vlin*0.9
      wgin=wlin*0.9
!
!-----MOMENTS AT INLET:
!
      void=1-(uinj/ulin)*(rinj**2/rs**2)
      q3in=(1-void)/cnt
!
!-----GAMMA DISTRIBUTION
!
      q2in=q3in/r32in
      q1in=q2in*(pkin+2)/(r32in*(pkin+1))
!
!-----TIME STEP
!
      delt=0.99*min(rs/ulin,dinj/vlin,dinj/wlin)
!
!-----NUMBER OF TIME STEPS
!
      ntsmx=tmax/delt
!
!-----DOMAIN SIZE
!
      xl=tmax*ulin*1.4
      yl=tmax*vlin*2.2
      zl=tmax*wlin*2.2
!
!-----CREATE GRID
!
      call grid
!
!-----ININTIALISATION DATA
!
      print*,'=================================================='
      print*,'ia, ja, ka'
      print*,ia,ja,ka
      print*,'nradic'
      print*,nradic
      print*,'beta, q3in, ulin'
      print*,beta,q3in,ulin
      print*,'delt, tmax'
      print*,delt,tmax
      print*,'=================================================='
      
      return
      end