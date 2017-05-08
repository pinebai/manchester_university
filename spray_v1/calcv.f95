!
!=====GAS MOMENTUM
!
      subroutine calcvg
      include "common.inc"
!
      ivel=2;iph=1;iq=3
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
      do ib=1,6
        call nbc(ib,vg)
      enddo
      call drag(vg,vq3)
      call pressure
      call vstress
      call temporal(vg0)
      call solve(vg,vgin)
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            apv(i,j,k)=ap(i,j,k)
          enddo
        enddo
      enddo
!
      return
      end
!
!=====Q0 MOMENTUM
!
      subroutine calcvq0
      include "common.inc"
!
      ivel=2;iph=2;iq=0
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
        call nbc(ib,vq0)
      enddo
      call dbc(1,vq0,0.)
      call injbc(vq0)
      call spbc
      call drag(vg,vq0)
      call vstress
      call source(vq3,bq0)
      call source(vq0,cq0)
      call temporal(vq00)
      call solve(vq0,vqin)
!
      return
      end
!
!=====Q1 MOMENTUM
!
      subroutine calcvq1
      include "common.inc"
!
      ivel=2;iph=2;iq=1
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
        call nbc(ib,vq1)
      enddo
      call dbc(1,vq1,0.)
      call injbc(vq1)
      call spbc
      call drag(vg,vq1)
      call vstress
      call source(vq3,bq1)
      call source(vq1,cq1)
      call temporal(vq10)
      call solve(vq1,vqin)
!
      return
      end
!
!=====Q2 MOMENTUM
!
      subroutine calcvq2
      include "common.inc"
!
      ivel=2;iph=2;iq=2
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
        call nbc(ib,vq2)
      enddo
      call dbc(1,vq2,0.)
      call injbc(vq2)
      call spbc
      call drag(vg,vq2)
      call vstress
      call source(vq3,bq2)
      call source(vq2,cq2)
      call temporal(vq20)
      call solve(vq2,vqin)
!
      return
      end
!
!=====Q3 MOMENTUM
!
      subroutine calcvq3
      include "common.inc"
!
      ivel=2;iph=2;iq=3
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
        call nbc(ib,vq3)
      enddo
      call dbc(1,vq3,0.)
      call injbc(vq3)
      call spbc
      call drag(vg,vq3)
      call vstress
      call temporal(vq30)
      call solve(vq3,vqin)
!
      return
      end
!
!=====Q4 MOMENTUM
!
      subroutine calcvq4
      include "common.inc"
!
      ivel=2;iph=2;iq=4
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
        call nbc(ib,vq4)
      enddo
      call dbc(1,vq4,0.)
      call injbc(vq4)
      call spbc
      call drag(vg,vq4)
      call vstress
      call source(vq3,bq4)
      call source(vq4,cq4)
      call temporal(vq40)
      call solve(vq4,vqin)
!
      return
      end
!
!=====Q5 MOMENTUM
!
      subroutine calcvq5
      include "common.inc"
!
      ivel=2;iph=2;iq=5
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
        call nbc(ib,vq5)
      enddo
      call dbc(1,vq5,0.)
      call injbc(vq5)
      call spbc
      call drag(vg,vq5)
      call vstress
      call source(vq3,bq5)
      call source(vq5,cq5)
      call temporal(vq50)
      call solve(vq5,vqin)
!
      return
      end
!
!=====Q6 MOMENTUM
!
      subroutine calcvq6
      include "common.inc"
!
      ivel=2;iph=2;iq=6
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
        call nbc(ib,vq6)
      enddo
      call dbc(1,vq6,0.)
      call injbc(vq6)
      call spbc
      call drag(vg,vq6)
      call vstress
      call source(vq3,bq6)
      call source(vq6,cq6)
      call temporal(vq60)
      call solve(vq6,vqin)
!
      return
      end
!
!=====Q7 MOMENTUM
!
      subroutine calcvq7
      include "common.inc"
!
      ivel=2;iph=2;iq=7
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
        call nbc(ib,vq7)
      enddo
      call dbc(1,vq7,0.)
      call injbc(vq7)
      call spbc
      call drag(vg,vq7)
      call vstress
      call source(vq3,bq7)
      call source(vq7,cq7)
      call temporal(vq70)
      call solve(vq7,vqin)
!
      return
      end