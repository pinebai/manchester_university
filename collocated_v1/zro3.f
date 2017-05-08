!
!=====ZERO DIFFUSION FOR MOMENT EQUATIONS
!
      subroutine zro3
      include "common.inc"

      do i=imn,imx
        do j=jmn,jmx
          dw(i,j)=0
          de(i,j)=0
          ds(i,j)=0
          dn(i,j)=0
        enddo
      enddo

      return
      end