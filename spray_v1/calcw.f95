!
!=====GAS MOMENTUM
!
      subroutine calcwg
      include "common.inc"
!
      ivel=3;iph=1;iq=3
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
        call nbc(ib,wg)
      enddo
      call drag(wg,wq3)
      call pressure
      call vstress
      call temporal(wg0)
      call solve(wg,wgin)
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            apw(i,j,k)=ap(i,j,k)
          enddo
        enddo
      enddo
!
      return
      end
!
!=====Q0 MOMENTUM
!
      subroutine calcwq0
      include "common.inc"
!
      ivel=3;iph=2;iq=0
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
        call nbc(ib,wq0)
      enddo
      call dbc(1,wq0,0.)
      call injbc(wq0)
      call spbc
      call drag(wg,wq0)
      call vstress
      call source(wq3,bq0)
      call source(wq0,cq0)
      call temporal(wq00)
      call solve(wq0,wqin)
!
      return
      end
!
!=====Q1 MOMENTUM
!
      subroutine calcwq1
      include "common.inc"
!
      ivel=3;iph=2;iq=1
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
        call nbc(ib,wq1)
      enddo
      call dbc(1,wq1,0.)
      call injbc(wq1)
      call spbc
      call drag(wg,wq1)
      call vstress
      call source(wq3,bq1)
      call source(wq1,cq1)
      call temporal(wq10)
      call solve(wq1,wqin)
!
      return
      end
!
!=====Q2 MOMENTUM
!
      subroutine calcwq2
      include "common.inc"
!
      ivel=3;iph=2;iq=2
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
        call nbc(ib,wq2)
      enddo
      call dbc(1,wq2,0.)
      call injbc(wq2)
      call spbc
      call drag(wg,wq2)
      call vstress
      call source(wq3,bq2)
      call source(wq2,cq2)
      call temporal(wq20)
      call solve(wq2,wqin)
!
      return
      end
!
!=====Q3 MOMENTUM
!
      subroutine calcwq3
      include "common.inc"
!
      ivel=3;iph=2;iq=3
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
        call nbc(ib,wq3)
      enddo
      call dbc(1,wq3,0.)
      call injbc(wq3)
      call spbc
      call drag(wg,wq3)
      call vstress
      call temporal(wq30)
      call solve(wq3,wqin)
!
      return
      end
!
!=====Q4 MOMENTUM
!
      subroutine calcwq4
      include "common.inc"
!
      ivel=3;iph=2;iq=4
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
        call nbc(ib,wq4)
      enddo
      call dbc(1,wq4,0.)
      call injbc(wq4)
      call spbc
      call drag(wg,wq4)
      call vstress
      call source(wq3,bq4)
      call source(wq4,cq4)
      call temporal(wq40)
      call solve(wq4,wqin)
!
      return
      end
!
!=====Q5 MOMENTUM
!
      subroutine calcwq5
      include "common.inc"
!
      ivel=3;iph=2;iq=5
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
        call nbc(ib,wq5)
      enddo
      call dbc(1,wq5,0.)
      call injbc(wq5)
      call spbc
      call drag(wg,wq5)
      call vstress
      call source(wq3,bq5)
      call source(wq5,cq5)
      call temporal(wq50)
      call solve(wq5,wqin)
!
      return
      end
!
!=====Q6 MOMENTUM
!
      subroutine calcwq6
      include "common.inc"
!
      ivel=3;iph=2;iq=6
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
        call nbc(ib,wq6)
      enddo
      call dbc(1,wq6,0.)
      call injbc(wq6)
      call spbc
      call drag(wg,wq6)
      call vstress
      call source(wq3,bq6)
      call source(wq6,cq6)
      call temporal(wq60)
      call solve(wq6,wqin)
!
      return
      end
!
!=====Q7 MOMENTUM
!
      subroutine calcwq7
      include "common.inc"
!
      ivel=3;iph=2;iq=7
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
        call nbc(ib,wq7)
      enddo
      call dbc(1,wq7,0.)
      call injbc(wq7)
      call spbc
      call drag(wg,wq7)
      call vstress
      call source(wq3,bq7)
      call source(wq7,cq7)
      call temporal(wq70)
      call solve(wq7,wqin)
!
      return
      end