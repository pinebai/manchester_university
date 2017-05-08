!
!=====GAS MOMENTUM
!
      subroutine calcug
      include "common.inc"
!
      ivel=1;iph=1;iq=3
!
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=dng0(i,j,k)*(1.-(4/3.)*pi*q30(i,j,k))
            ff(i,j,k)=dng(i,j,k)*(1.-(4/3.)*pi*q3(i,j,k))
            dd(i,j,k)=dvsg(i,j,k)*(1.-(4/3.)*pi*q3(i,j,k))
          enddo
        enddo
      enddo
!
      call zerocoef
      call convection(ug,vg,wg)
      call diffusion
      call hybrid
      do ib=2,6
        call nbc(ib,ug)
      enddo
      call dbc(1,ug,0.)
      call drag(ug,uq3)
      call pressure
      call vstress
      call temporal(ug0)
      call solve(ug,ugin)
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            apu(i,j,k)=ap(i,j,k)
          enddo
        enddo
      enddo
!
      return
      end
!
!=====Q0 MOMENTUM
!
      subroutine calcuq0
      include "common.inc"
!
      ivel=1;iph=2;iq=0
!
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=q00(i,j,k)
            ff(i,j,k)=q0(i,j,k)
            dd(i,j,k)=0.7*dvsl(i,j,k)*q0(i,j,k)/dnl(i,j,k)
          enddo
        enddo
      enddo
!
      call zerocoef
      call convection(uq0,vq0,wq0)
      call diffusion
      call hybrid
      do ib=2,6
        call nbc(ib,uq0)
      enddo
      call dbc(1,uq0,0.)
      call injbc(uq0)
      call spbc
      call drag(ug,uq0)
      call vstress
      call source(uq0,bq0)
      call temporal(uq00)
      call solve(uq0,uqin)
!
      return
      end
!
!=====Q1 MOMENTUM
!
      subroutine calcuq1
      include "common.inc"
!
      ivel=1;iph=2;iq=1
!
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=q10(i,j,k)
            ff(i,j,k)=q1(i,j,k)
            dd(i,j,k)=0.7*dvsl(i,j,k)*q1(i,j,k)/dnl(i,j,k)
          enddo
        enddo
      enddo
!
      call zerocoef
      call convection(uq1,vq1,wq1)
      call diffusion
      call hybrid
      do ib=2,6
        call nbc(ib,uq1)
      enddo
      call dbc(1,uq1,0.)
      call injbc(uq1)
      call spbc
      call drag(ug,uq1)
      call vstress
      call source(uq1,bq1)
      call temporal(uq10)
      call solve(uq1,uqin)
!
      return
      end
!
!=====Q2 MOMENTUM
!
      subroutine calcuq2
      include "common.inc"
!
      ivel=1;iph=2;iq=2
!
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=q20(i,j,k)
            ff(i,j,k)=q2(i,j,k)
            dd(i,j,k)=0.7*dvsl(i,j,k)*q2(i,j,k)/dnl(i,j,k)
          enddo
        enddo
      enddo
!
      call zerocoef
      call convection(uq2,vq2,wq2)
      call diffusion
      call hybrid
      do ib=2,6
        call nbc(ib,uq2)
      enddo
      call dbc(1,uq2,0.)
      call injbc(uq2)
      call spbc
      call drag(ug,uq2)
      call vstress
      call source(uq2,bq2)
      call temporal(uq20)
      call solve(uq2,uqin)
!
      return
      end
!
!=====Q3 MOMENTUM
!
      subroutine calcuq3
      include "common.inc"
!
      ivel=1;iph=2;iq=3
!
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=(4/3.)*pi*q30(i,j,k)*dnl0(i,j,k)
            ff(i,j,k)=(4/3.)*pi*q3(i,j,k)*dnl(i,j,k)
            dd(i,j,k)=0.7*dvsl(i,j,k)*(4/3.)*pi*q3(i,j,k)
          enddo
        enddo
      enddo
!
      call zerocoef
      call convection(uq3,vq3,wq3)
      call diffusion
      call hybrid
      do ib=2,6
        call nbc(ib,uq3)
      enddo
      call dbc(1,uq3,0.)
      call injbc(uq3)
      call spbc
      call drag(ug,uq3)
      call vstress
      call temporal(uq30)
      call solve(uq3,uqin)
!
      return
      end
!
!=====Q4 MOMENTUM
!
      subroutine calcuq4
      include "common.inc"
!
      ivel=1;iph=2;iq=4
!
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=q40(i,j,k)
            ff(i,j,k)=q4(i,j,k)
            dd(i,j,k)=0.7*dvsl(i,j,k)*q4(i,j,k)/dnl(i,j,k)
          enddo
        enddo
      enddo
!
      call zerocoef
      call convection(uq4,vq4,wq4)
      call diffusion
      call hybrid
      do ib=2,6
        call nbc(ib,uq4)
      enddo
      call dbc(1,uq4,0.)
      call injbc(uq4)
      call spbc
      call drag(ug,uq4)
      call vstress
      call source(uq4,bq4)
      call temporal(uq40)
      call solve(uq4,uqin)
!
      return
      end
!
!=====Q5 MOMENTUM
!
      subroutine calcuq5
      include "common.inc"
!
      ivel=1;iph=2;iq=5
!
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=q50(i,j,k)
            ff(i,j,k)=q5(i,j,k)
            dd(i,j,k)=0.7*dvsl(i,j,k)*q5(i,j,k)/dnl(i,j,k)
          enddo
        enddo
      enddo
!
      call zerocoef
      call convection(uq5,vq5,wq5)
      call diffusion
      call hybrid
      do ib=2,6
        call nbc(ib,uq5)
      enddo
      call dbc(1,uq5,0.)
      call injbc(uq5)
      call spbc
      call drag(ug,uq5)
      call vstress
      call source(uq5,bq5)
      call temporal(uq50)
      call solve(uq5,uqin)
!
      return
      end
!
!=====Q6 MOMENTUM
!
      subroutine calcuq6
      include "common.inc"
!
      ivel=1;iph=2;iq=6
!
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=q60(i,j,k)
            ff(i,j,k)=q6(i,j,k)
            dd(i,j,k)=0.7*dvsl(i,j,k)*q6(i,j,k)/dnl(i,j,k)
          enddo
        enddo
      enddo
!
      call zerocoef
      call convection(uq6,vq6,wq6)
      call diffusion
      call hybrid
      do ib=2,6
        call nbc(ib,uq6)
      enddo
      call dbc(1,uq6,0.)
      call injbc(uq6)
      call spbc
      call drag(ug,uq6)
      call vstress
      call source(uq6,bq6)
      call temporal(uq60)
      call solve(uq6,uqin)
!
      return
      end
!
!=====Q7 MOMENTUM
!
      subroutine calcuq7
      include "common.inc"
!
      ivel=1;iph=2;iq=7
!
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=q70(i,j,k)
            ff(i,j,k)=q7(i,j,k)
            dd(i,j,k)=0.7*dvsl(i,j,k)*q7(i,j,k)/dnl(i,j,k)
          enddo
        enddo
      enddo
!
      call zerocoef
      call convection(uq7,vq7,wq7)
      call diffusion
      call hybrid
      do ib=2,6
        call nbc(ib,uq7)
      enddo
      call dbc(1,uq7,0.)
      call injbc(uq7)
      call spbc
      call drag(ug,uq7)
      call vstress
      call source(uq7,bq7)
      call temporal(uq70)
      call solve(uq7,uqin)
!
      return
      end