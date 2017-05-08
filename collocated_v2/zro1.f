!
!=====MAIN ARRAY ZERO
!
      subroutine zro1
      include "common.inc"
      
      do i=imin,imax
        do j=jmin,jmax
          ug(i,j)=0.
          ugx(i,j)=0.
          vg(i,j)=0.
          vgx(i,j)=0.
          ul1(i,j)=0.
          vl1(i,j)=0.
          ul2(i,j)=0.
          vl2(i,j)=0.
          ul3(i,j)=0.
          vl3(i,j)=0.
          p(i,j)=0.
          pc(i,j)=0.
          pk(i,j)=0.
          qm1(i,j)=0.
          q0(i,j)=0.
          q1(i,j)=0.
          q10(i,j)=0.
          q2(i,j)=0.
          q20(i,j)=0.
          q3(i,j)=0.
          q30(i,j)=0.
          r21(i,j)=0.
          r32(i,j)=0.
          du(i,j)=0.
          dv(i,j)=0.
          dng(i,j)=0.
          dng0(i,j)=0.
          dnl(i,j)=0.
          dnl0(i,j)=0.
          drg(i,j)=0.
          drg1a(i,j)=0.
          drg1b(i,j)=0.
          drg2a(i,j)=0.
          drg2b(i,j)=0.
          drg3a(i,j)=0.
          drg3b(i,j)=0.
          drgu1(i,j)=0.
          drgv1(i,j)=0.
          drgu2(i,j)=0.
          drgv2(i,j)=0.
          drgu3(i,j)=0.
          drgv3(i,j)=0.
        enddo
      enddo

      return
      end