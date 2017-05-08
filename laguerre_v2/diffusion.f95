!
!=====DIFFUSION COEFFICIENTS
!
      subroutine diffusion
      include "common.inc"
      
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            dw(i,j,k)=area(1,i,j,k)*f(1,i,j,k,dd)/(xc(i)-xc(i-1))
            de(i,j,k)=area(2,i,j,k)*f(2,i,j,k,dd)/(xc(i+1)-xc(i))
            ds(i,j,k)=area(3,i,j,k)*f(3,i,j,k,dd)/(yc(j)-yc(j-1))
            dn(i,j,k)=area(4,i,j,k)*f(4,i,j,k,dd)/(yc(j+1)-yc(j))
            db(i,j,k)=area(5,i,j,k)*f(5,i,j,k,dd)/(zc(k)-zc(k-1))
            dt(i,j,k)=area(6,i,j,k)*f(6,i,j,k,dd)/(zc(k+1)-zc(k))
          enddo
        enddo
      enddo
      
      return
      end