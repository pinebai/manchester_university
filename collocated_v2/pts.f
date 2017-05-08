!
!=====SAVE PREVIOUS TIME STEP VALUES
!
      subroutine pts

      include "common.inc"
      
      call idx
      do i=imn,imx
        do j=jmn,jmx
!
!-----U - VELOCITIES
!
          ug0(i,j)=ug(i,j)
          ul10(i,j)=ul1(i,j)
          ul20(i,j)=ul2(i,j)
          ul30(i,j)=ul3(i,j)
!
!-----V - VELOCITIES
!
          vg0(i,j)=vg(i,j)
          vl10(i,j)=vl1(i,j)
          vl20(i,j)=vl2(i,j)
          vl30(i,j)=vl3(i,j)
!
!-----DENSITIES
!
          dng0(i,j)=dng(i,j)
          dnl0(i,j)=dnl(i,j)
!
!-----MOMENTS
!
          q10(i,j)=q1(i,j)
          q20(i,j)=q2(i,j)
          q30(i,j)=q3(i,j)
        enddo
      enddo

      return
      end