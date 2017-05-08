!
!=====CONVECTION COEFFICIENTS
!
      subroutine convection(uvel,vvel,wvel)
      include "common.inc"
      real uvel(ia,ja,ka),vvel(ia,ja,ka),wvel(ia,ja,ka)
      
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            fw(i,j,k)=area(1,i,j,k)*f2(1,i,j,k,ff,uvel)
            fe(i,j,k)=area(2,i,j,k)*f2(2,i,j,k,ff,uvel)
            fs(i,j,k)=area(3,i,j,k)*f2(3,i,j,k,ff,vvel)
            fn(i,j,k)=area(4,i,j,k)*f2(4,i,j,k,ff,vvel)
            fb(i,j,k)=area(5,i,j,k)*f2(5,i,j,k,ff,wvel)
            ft(i,j,k)=area(6,i,j,k)*f2(6,i,j,k,ff,wvel)
          enddo
        enddo
      enddo
      
      return
      end