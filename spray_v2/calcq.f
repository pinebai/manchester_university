!
!=====Q0 MOMENT
!
      subroutine calcq0
      include "common.inc"
      
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff00(i,j,k)=dnl00(i,j,k)
            ff0(i,j,k)=dnl0(i,j,k)
            ff(i,j,k)=dnl(i,j,k)
          enddo
        enddo
      enddo
      ivel=0;iph=2;iq=0;iequ=17
      call zerocoef
      call convection(ul0,vl0,wl0)
C       call hybrid
      call quick(q0)
      do ib=2,6
        call nbc(ib,q0)
      enddo
      call dbc(1,q0,0.)
      call injbc(q0)
      call spbc
      call temporal(q000,q00)
      call solve(q0)
      call cut(q0,q0in)

      return
      end
!
!=====Q1 MOMENT
!
      subroutine calcq1
      include "common.inc"
      
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff00(i,j,k)=dnl00(i,j,k)
            ff0(i,j,k)=dnl0(i,j,k)
            ff(i,j,k)=dnl(i,j,k)
          enddo
        enddo
      enddo
      ivel=0;iph=2;iq=1;iequ=18
      call zerocoef
      call convection(ul1,vl1,wl1)
C       call hybrid
      call quick(q1)
      do ib=2,6
        call nbc(ib,q1)
      enddo
      call dbc(1,q1,0.)
      call injbc(q1)
      call spbc
      call temporal(q100,q10)
      call solve(q1)
      call cut(q1,q1in)

      return
      end
!
!=====Q2 MOMENT
!
      subroutine calcq2
      include "common.inc"
      
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff00(i,j,k)=dnl00(i,j,k)
            ff0(i,j,k)=dnl0(i,j,k)
            ff(i,j,k)=dnl(i,j,k)
          enddo
        enddo
      enddo
      ivel=0;iph=2;iq=2;iequ=19
      call zerocoef
      call convection(ul2,vl2,wl2)
C       call hybrid
      call quick(q2)
      do ib=2,6
        call nbc(ib,q2)
      enddo
      call dbc(1,q2,0.)
      call injbc(q2)
      call spbc
      call temporal(q200,q20)
      call solve(q2)
      call cut(q2,q2in)

      return
      end
!
!=====Q3 MOMENT
!
      subroutine calcq3
      include "common.inc"
      
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff00(i,j,k)=dnl00(i,j,k)
            ff0(i,j,k)=dnl0(i,j,k)
            ff(i,j,k)=dnl(i,j,k)
          enddo
        enddo
      enddo
      ivel=0;iph=2;iq=3;iequ=20
      call zerocoef
      call convection(ul3,vl3,wl3)
C       call hybrid
      call quick(q3)
      do ib=2,6
        call nbc(ib,q3)
      enddo
      call dbc(1,q3,0.)
      call injbc(q3)
      call spbc
      call temporal(q300,q30)
      call solve(q3)
      call cut(q3,q3in)
      
      return
      end