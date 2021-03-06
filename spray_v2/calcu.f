!
!=====GAS MOMENTUM
!
      subroutine calcug
      include "common.inc"

      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff00(i,j,k)=dng00(i,j,k)*(1.-cnt*q300(i,j,k))
            ff0(i,j,k)=dng0(i,j,k)*(1.-cnt*q30(i,j,k))
            ff(i,j,k)=dng(i,j,k)*(1.-cnt*q3(i,j,k))
            dd(i,j,k)=dvsg(i,j,k)*(1.-cnt*q3(i,j,k))
          enddo
        enddo
      enddo
      ivel=1;iph=1;iq=3;iequ=1
      call zerocoef
      call convection(ug,vg,wg)
      call diffusion
C       call hybrid
      call quick(ug)
      do ib=1,6
        call nbc(ib,ug)
      enddo
      call drag(ug,ul3)
      call psource
      call temporal(ug00,ug0)
      call solve(ug)
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            apu(i,j,k)=ap(i,j,k)
          enddo
        enddo
      enddo

      return
      end
!
!=====Q0 MOMENTUM
!
      subroutine calcul0
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
      ivel=1;iph=2;iq=0;iequ=2
      call zerocoef
      call convection(ul0,vl0,wl0)
      call diffusion
      call hybrid
C       call quick(ul0)
      do ib=2,6
        call nbc(ib,ul0)
      enddo
      call dbc(1,ul0,0.)
      call injbc(ul0)
      call spbc
      call drag(ug,ul0)
      call temporal(ul000,ul00)
      call solve(ul0)

      return
      end
!
!=====Q1 MOMENTUM
!
      subroutine calcul1
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
      ivel=1;iph=2;iq=1;iequ=3
      call zerocoef
      call convection(ul1,vl1,wl1)
      call diffusion
      call hybrid
C       call quick(ul1)
      do ib=2,6
        call nbc(ib,ul1)
      enddo
      call dbc(1,ul1,0.)
      call injbc(ul1)
      call spbc
      call drag(ug,ul1)
      call temporal(ul100,ul10)
      call solve(ul1)

      return
      end
!
!=====Q2 MOMENTUM
!
      subroutine calcul2
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
      ivel=1;iph=2;iq=2;iequ=4
      call zerocoef
      call convection(ul2,vl2,wl2)
      call diffusion
      call hybrid
C       call quick(ul2)
      do ib=2,6
        call nbc(ib,ul2)
      enddo
      call dbc(1,ul2,0.)
      call injbc(ul2)
      call spbc
      call drag(ug,ul2)
      call temporal(ul200,ul20)
      call solve(ul2)

      return
      end
!
!=====Q3 MOMENTUM
!
      subroutine calcul3
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
      ivel=1;iph=2;iq=3;iequ=5
      call zerocoef
      call convection(ul3,vl3,wl3)
      call diffusion
      call hybrid
C       call quick(ul3)
      do ib=2,6
        call nbc(ib,ul3)
      enddo
      call dbc(1,ul3,0.)
      call injbc(ul3)
      call spbc
      call drag(ug,ul3)
      call temporal(ul300,ul30)
      call solve(ul3)
      
      return
      end