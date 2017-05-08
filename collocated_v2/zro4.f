!
!=====ZERO PARTICULAR ARRAY
!
      subroutine zro4(phi)
      include "common.inc"

      real phi(ia,ja)
      
      do i=imin,imax
        do j=jmin,jmax
          phi(i,j)=0.
        enddo
      enddo

      return
      end