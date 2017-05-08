!
!=====SAVE PREVIOUS TIME STEP VALUES
!
      subroutine store
      include "common.inc"
      
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
!
!-----U - VELOCITIES
!
            ug0(i,j,k)=ug(i,j,k)
            ul10(i,j,k)=ul1(i,j,k)
            ul20(i,j,k)=ul2(i,j,k)
            ul30(i,j,k)=ul3(i,j,k)
!
!-----V - VELOCITIES
!
            vg0(i,j,k)=vg(i,j,k)
            vl10(i,j,k)=vl1(i,j,k)
            vl20(i,j,k)=vl2(i,j,k)
            vl30(i,j,k)=vl3(i,j,k)
!
!-----W - VELOCITIES
!
            wg0(i,j,k)=wg(i,j,k)
            wl10(i,j,k)=wl1(i,j,k)
            wl20(i,j,k)=wl2(i,j,k)
            wl30(i,j,k)=wl3(i,j,k)
!
!-----DENSITIES
!
            dng0(i,j,k)=dng(i,j,k)
            dnl0(i,j,k)=dnl(i,j,k)
!
!-----MOMENTS
!
            q10(i,j,k)=q1(i,j,k)
            q20(i,j,k)=q2(i,j,k)
            q30(i,j,k)=q3(i,j,k)
          enddo
        enddo
      enddo

      return
      end