!
!=====GAS MOMENTUM
!
      subroutine calcwg
      include "common.inc"
      
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff00(i,j,k)=dng00(i,j,k)*(1-cnt*q300(i,j,k))
            ff0(i,j,k)=dng0(i,j,k)*(1-cnt*q30(i,j,k))
            ff(i,j,k)=dng(i,j,k)*(1-cnt*q3(i,j,k))
            dd(i,j,k)=dvsg(i,j,k)*(1-cnt*q3(i,j,k))
          enddo
        enddo
      enddo
      ivel=3;iph=1;iq=3;iequ=11
      call zerocoef
      call convection(ug,vg,wg)
      call diffusion
C       call hybrid
      call quick(wg)
      do ib=1,6
        call nbc(ib,wg)
      enddo
      call drag(wg,wl3)
      call psource
      call temporal(wg00,wg0)
      call solve(wg)
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            apw(i,j,k)=ap(i,j,k)
          enddo
        enddo
      enddo

      return
      end
!
!=====Q0 MOMENTUM
!
      subroutine calcwl0
      include "common.inc"
      
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff00(i,j,k)=dnl00(i,j,k)*q000(i,j,k)
            ff0(i,j,k)=dnl0(i,j,k)*q00(i,j,k)
            ff(i,j,k)=dnl(i,j,k)*q0(i,j,k)
            dd(i,j,k)=dvsl(i,j,k)*q0(i,j,k)
          enddo
        enddo
      enddo
      ivel=3;iph=2;iq=0;iequ=12
      call zerocoef
      call convection(ul0,vl0,wl0)
      call diffusion
      call hybrid
C       call quick(wl0)
      do ib=2,6
        call nbc(ib,wl0)
      enddo
      call dbc(1,wl0,0.)
      call injbc(wl0)
      call spbc
      call drag(wg,wl0)
      call temporal(wl000,wl00)
      call solve(wl0)

      return
      end
!
!=====Q1 MOMENTUM
!
      subroutine calcwl1
      include "common.inc"
      
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff00(i,j,k)=dnl00(i,j,k)*q100(i,j,k)
            ff0(i,j,k)=dnl0(i,j,k)*q10(i,j,k)
            ff(i,j,k)=dnl(i,j,k)*q1(i,j,k)
            dd(i,j,k)=dvsl(i,j,k)*q1(i,j,k)
          enddo
        enddo
      enddo
      ivel=3;iph=2;iq=1;iequ=13
      call zerocoef
      call convection(ul1,vl1,wl1)
      call diffusion
      call hybrid
C       call quick(wl1)
      do ib=2,6
        call nbc(ib,wl1)
      enddo
      call dbc(1,wl1,0.)
      call injbc(wl1)
      call spbc
      call drag(wg,wl1)
      call temporal(wl100,wl10)
      call solve(wl1)

      return
      end
!
!=====Q2 MOMENTUM
!
      subroutine calcwl2
      include "common.inc"
      
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff00(i,j,k)=dnl00(i,j,k)*q200(i,j,k)
            ff0(i,j,k)=dnl0(i,j,k)*q20(i,j,k)
            ff(i,j,k)=dnl(i,j,k)*q2(i,j,k)
            dd(i,j,k)=dvsl(i,j,k)*q2(i,j,k)
          enddo
        enddo
      enddo
      ivel=3;iph=2;iq=2;iequ=14
      call zerocoef
      call convection(ul2,vl2,wl2)
      call diffusion
      call hybrid
C       call quick(wl2)
      do ib=2,6
        call nbc(ib,wl2)
      enddo
      call dbc(1,wl2,0.)
      call injbc(wl2)
      call spbc
      call drag(wg,wl2)
      call temporal(wl200,wl20)
      call solve(wl2)

      return
      end
!
!=====Q3 MOMENTUM
!
      subroutine calcwl3
      include "common.inc"
      
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff00(i,j,k)=dnl00(i,j,k)*cnt*q300(i,j,k)
            ff0(i,j,k)=dnl0(i,j,k)*cnt*q30(i,j,k)
            ff(i,j,k)=dnl(i,j,k)*cnt*q3(i,j,k)
            dd(i,j,k)=dvsl(i,j,k)*cnt*q3(i,j,k)
          enddo
        enddo
      enddo
      ivel=3;iph=2;iq=3;iequ=15
      call zerocoef
      call convection(ul3,vl3,wl3)
      call diffusion
      call hybrid
C       call quick(wl3)
      do ib=2,6
        call nbc(ib,wl3)
      enddo
      call dbc(1,wl3,0.)
      call injbc(wl3)
      call spbc
      call drag(wg,wl3)
      call temporal(wl300,wl30)
      call solve(wl3)

      return
      end