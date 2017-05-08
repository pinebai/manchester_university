!
!=====INTERMEDIATE ARRAY ZERO
!
      subroutine zro2
      include "common.inc"
      
      do i=imin,imax
        do j=jmin,jmax
          aw(i,j)=0.
          ae(i,j)=0.
          as(i,j)=0.
          an(i,j)=0.
          ap(i,j)=0.
          ap0(i,j)=0.
          su(i,j)=0.
          sp(i,j)=0.
          drg(i,j)=0.
        enddo
      enddo

      return
      end