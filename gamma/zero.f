!
!=====MAIN ARRAY ZERO
!
      subroutine zeromain
      include "common.inc"
      
      do i=1,ia
        do j=1,ja
          do k=1,ka
            bq1(i,j,k)=0
            bq2(i,j,k)=0
            pk(i,j,k)=0
            r21(i,j,k)=0
            r32(i,j,k)=0
            qm1(i,j,k)=0
            flux(i,j,k)=0
            p(i,j,k)=0
            dw(i,j,k)=0
            de(i,j,k)=0
            ds(i,j,k)=0
            dn(i,j,k)=0
            db(i,j,k)=0
            dt(i,j,k)=0
            dpx(i,j,k)=0
            dpy(i,j,k)=0
            dpz(i,j,k)=0
            apu(i,j,k)=great
            apv(i,j,k)=great
            apw(i,j,k)=great
!
!-----CURRENT TIME STEP
!
            ug(i,j,k)=0
            vg(i,j,k)=0
            wg(i,j,k)=0
            ul1(i,j,k)=0
            ul2(i,j,k)=0
            ul3(i,j,k)=0
            vl1(i,j,k)=0
            vl2(i,j,k)=0
            vl3(i,j,k)=0
            wl1(i,j,k)=0
            wl2(i,j,k)=0
            wl3(i,j,k)=0
            q1(i,j,k)=0
            q2(i,j,k)=0
            q3(i,j,k)=0
            dng(i,j,k)=0
            dnl(i,j,k)=0
!
!-----PREVIOUS TIME STEP
!
            ug0(i,j,k)=0
            vg0(i,j,k)=0
            wg0(i,j,k)=0
            ul10(i,j,k)=0
            ul20(i,j,k)=0
            ul30(i,j,k)=0
            vl10(i,j,k)=0
            vl20(i,j,k)=0
            vl30(i,j,k)=0
            wl10(i,j,k)=0
            wl20(i,j,k)=0
            wl30(i,j,k)=0
            q10(i,j,k)=0
            q20(i,j,k)=0
            q30(i,j,k)=0
            dng0(i,j,k)=0
            dnl0(i,j,k)=0
          enddo
        enddo
      enddo

      do iequ=1,20
        res(iequ)=0
      enddo
      
      gsmax0=0;gsmax=0;gdelta=0

      return
      end
!
!=====INTERMEDIATE ARRAY ZERO
!
      subroutine zerocoef
      include "common.inc"
      
      do i=1,ia
        do j=1,ja
          do k=1,ka
            aw(i,j,k)=0
            ae(i,j,k)=0
            as(i,j,k)=0
            an(i,j,k)=0
            ab(i,j,k)=0
            at(i,j,k)=0
            ap(i,j,k)=0
            su(i,j,k)=0
            sp(i,j,k)=0
          enddo
        enddo
      enddo

      return
      end
!
!=====ZERO PARTICULAR ARRAY
!
      subroutine zerovar(phi)
      include "common.inc"
      real phi(ia,ja,ka)
      
      do i=1,ia
        do j=1,ja
          do k=1,ka
            phi(i,j,k)=0
          enddo
        enddo
      enddo

      return
      end
!
!=====ZERO BOUNDARY COEFFICIENTS
!
      subroutine zerobcoef
      include "common.inc"
!
!-----WEST
!
        i=imn
        do j=jmn,jmx
          do k=kmn,kmx
            aw(i,j,k)=0
          enddo
        enddo
!
!-----EAST
!
        i=imx
        do j=jmn,jmx
          do k=kmn,kmx
            ae(i,j,k)=0
          enddo
        enddo
!
!-----SOUTH
!
        j=jmn
        do i=imn,imx
          do k=kmn,kmx
            as(i,j,k)=0
          enddo
        enddo
!
!-----NORTH
!
        j=jmx
        do i=imn,imx
          do k=kmn,kmx
           an(i,j,k)=0
          enddo
        enddo
!
!-----BOTTOM
!
        k=kmn
        do i=imn,imx
          do j=jmn,jmx
            ab(i,j,k)=0
          enddo
        enddo
!
!-----TOP
!
        k=kmx
        do i=imn,imx
          do j=jmn,jmx
            at(i,j,k)=0
          enddo
        enddo

      return
      end