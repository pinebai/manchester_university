!
!=====INJECTION INLET PARAMETERS
!
      subroutine inj
      include "common.inc"
!
!-----DURATION
!
      tmax=0.5e-3
!
!-----CONE HALF-ANGLE
!
      beta=13.*(pi/180.)
!
!-----INJECTOR RADIUS
!
      rinj=1.e-4
      xinj=0.5e-3
      rs=rinj+xinj*tan(beta)
!
!-----NUMBER OF INJECTION CELLS
!
      ninjc=3
!
!-----SIZE OF INJECTOR CELLS
!
      dxc=rs/ninjc!xinj
      dyc=rs/ninjc
!
!-----REFERENCE LOCATION OF INJECTOR
!
      iinj=imn+2
      jinj=jmn+2
!
!-----CREATE GRID
!
      call grd
!
!-----INJECTOR PRESSURE
!
      pinj=10.e+6
!
!-----INJECTOR DRAG COEFFICIENT
!
      cd=0.7
!
!-----VELOCITIES AT INLET
!
      i=iinj
      j=jinj
      uinj=cd*(2*(pinj-p(imx,jmx))
     &     /dnl(i,j))**0.5
      ulin=uinj
     &     *((-(1-dng(i,j)/dnl(i,j))
     &     +((1-dng(i,j)/dnl(i,j))**2
     &     +4*dng(i,j)*rs**2
     &     /(dnl(i,j)*rinj**2))**0.5)
     &     /(2*dng(i,j)*rs**2
     &     /(dnl(i,j)*rinj**2)))
      vlin=ulin*tan(beta)
!
!-----SMR AT INLET
!
      r32in=7.e-6
!
!-----PARAMETER K AT INLET
!
      pkin=3.
!
!-----MOMENTS AT INLET
!
      void=1.-(uinj/ulin)*(rinj**2/rs**2)
      q3in=(1.-void)/cnt
      q2in=q3in/r32in
      q1in=q2in/(r32in*((pkin+1)/(pkin+2)))
!
!-----TIME STEP
!
      delt=0.5*dyc/ulin
!
!-----NUMBER OF TIME STEPS
!
      ntsmx=tmax/delt
!
!-----ON SCREEN INITIALISATION READOUT
!
      print*,'-------------------------------'
      print*,'   ia',ia
      print*,'   ja',ja
      print*,'ninjc',ninjc
      print*,'ntsmx',ntsmx
      print*,'itcmx',itcmx
      print*,'  dyc',dyc
      print*,'   xl',xl
      print*,'   yl',yl
      print*,' ulin',ulin
      print*,' vlin',vlin
      print*,' q3in',q3in
      print*,'r32in',r32in
      print*,' pkin',pkin
      print*,' delt',delt
      print*,' pgas',p(imn,jmn)
      print*,' pinj',pinj
      print*,' beta',beta
      print*,'-------------------------------'

      return
      end