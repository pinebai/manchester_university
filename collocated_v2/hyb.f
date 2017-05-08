!
!=====HYBRID SCHEME
!
      subroutine hyb
      include "common.inc"

      do i=imn,imx
        do j=jmn,jmx
          aw(i,j)=amax1(fw(i,j),(dw(i,j)+fw(i,j)/2),0.)
          ae(i,j)=amax1(-fe(i,j),(de(i,j)-fe(i,j)/2),0.)
          as(i,j)=amax1(fs(i,j),(ds(i,j)+fs(i,j)/2),0.)
          an(i,j)=amax1(-fn(i,j),(dn(i,j)-fn(i,j)/2),0.)
        enddo
      enddo

      return
      end