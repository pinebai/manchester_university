!
!=====GAS MOMENTUM
!
      subroutine calcvg
      include "common.inc"
      
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=dng0(i,j,k)*(1.-cnt*q30(i,j,k))
            ff(i,j,k)=dng(i,j,k)*(1.-cnt*q3(i,j,k))
            dd(i,j,k)=dvsg(i,j,k)*(1.-cnt*q3(i,j,k))
          enddo
        enddo
      enddo
      ivel=2;iph=1;iq=3;iequ=6
      call zerocoef
      call convection(ug,vg,wg)
      call diffusion
      call hybrid
C       call quick(vg)
C       do ib=1,6
C         call nbc(ib,vg)
C       enddo
      do ib=2,6
        call nbc(ib,vg)
      enddo
      call dbc(1,vg,0.)
      call injbc2(vg,nradic)
      call spbc2(nradic)
      call drag(vg,vl3)
      call pressure
      call temporal(vg0)
      call solve(vg)
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            apv(i,j,k)=ap(i,j,k)
          enddo
        enddo
      enddo

      return
      end
!
!=====Q1 MOMENTUM
!
      subroutine calcvl1
      include "common.inc"
      
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=dnl0(i,j,k)*q10(i,j,k)
            ff(i,j,k)=dnl(i,j,k)*q1(i,j,k)
            dd(i,j,k)=dvsl(i,j,k)*q1(i,j,k)
          enddo
        enddo
      enddo
      ivel=2;iph=2;iq=1;iequ=8
      call zerocoef
      call convection(ul1,vl1,wl1)
      call diffusion
      call hybrid
C       call quick(vl1)
      do ib=2,6
        call nbc(ib,vl1)
      enddo
      call dbc(1,vl1,0.)
      call injbc2(vl1,nradic)
      call spbc2(nradic)
      call drag(vg,vl1)
      call source(vl1,bq1)
      call temporal(vl10)
      call solve(vl1)

      return
      end
!
!=====Q2 MOMENTUM
!
      subroutine calcvl2
      include "common.inc"
      
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=dnl0(i,j,k)*q20(i,j,k)
            ff(i,j,k)=dnl(i,j,k)*q2(i,j,k)
            dd(i,j,k)=dvsl(i,j,k)*q2(i,j,k)
          enddo
        enddo
      enddo
      ivel=2;iph=2;iq=2;iequ=9
      call zerocoef
      call convection(ul2,vl2,wl2)
      call diffusion
      call hybrid
C       call quick(vl2)
      do ib=2,6
        call nbc(ib,vl2)
      enddo
      call dbc(1,vl2,0.)
      call injbc2(vl2,nradic)
      call spbc2(nradic)
      call drag(vg,vl2)
      call source(vl2,bq2)
      call temporal(vl20)
      call solve(vl2)

      return
      end
!
!=====Q3 MOMENTUM
!
      subroutine calcvl3
      include "common.inc"
      
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=dnl0(i,j,k)*cnt*q30(i,j,k)
            ff(i,j,k)=dnl(i,j,k)*cnt*q3(i,j,k)
            dd(i,j,k)=dvsl(i,j,k)*cnt*q3(i,j,k)
          enddo
        enddo
      enddo
      ivel=2;iph=2;iq=3;iequ=10
      call zerocoef
      call convection(ul3,vl3,wl3)
      call diffusion
      call hybrid
C       call quick(vl3)
      do ib=2,6
        call nbc(ib,vl3)
      enddo
      call dbc(1,vl3,0.)
      call injbc2(vl3,nradic)
      call spbc2(nradic)
      call drag(vg,vl3)
      call temporal(vl30)
      call solve(vl3)

      return
      end