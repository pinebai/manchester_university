!
!-----CONITINUTY IMBALANCE TERM USING PSEUDO-VELOCITIES
!
      subroutine rcc(uvel,vvel)
      include "common.inc"

      real uvel(ia,ja),vvel(ia,ja)
      
      do i=imn,imx
        do j=jmn,jmx
          su(i,j)=-(0.5*uvel(i+1,j)*dng(i+1,j)
     &            *dase(i,j)*(1-cnt*q3(i+1,j))
     &            +0.25*du(i,j)*dng(i,j)*dase(i,j)*(1-cnt*q3(i,j))
     &            *(p(i+2,j)-3*p(i+1,j)+3*p(i,j)-p(i-1,j)))
     &            +(0.5*uvel(i-1,j)*dng(i-1,j)
     &            *dasw(i,j)*(1-cnt*q3(i-1,j))
     &            +0.25*du(i,j)*dng(i,j)*dasw(i,j)*(1-cnt*q3(i,j))
     &            *(p(i+1,j)-3*p(i,j)+3*p(i-1,j)-p(i-2,j)))
     &            -(0.5*vvel(i,j+1)*dng(i,j+1)
     &            *dasn(i,j)*(1-cnt*q3(i,j+1))
     &            +0.25*dv(i,j)*dng(i,j)*dasn(i,j)*(1-cnt*q3(i,j))
     &            *(p(i,j+2)-3*p(i,j+1)+3*p(i,j)-p(i,j-1)))
     &            +(0.5*vvel(i,j-1)*dng(i,j-1)
     &            *dass(i,j)*(1-cnt*q3(i,j-1))
     &            +0.25*dv(i,j)*dng(i,j)*dass(i,j)*(1-cnt*q3(i,j))
     &            *(p(i,j+1)-3*p(i,j)+3*p(i,j-1)-p(i,j-2)))
     &            +(dng0(i,j)*qq0(i,j)-dng(i,j)*qq(i,j))
     &            *dvs(i,j)/delt
        enddo
      enddo

      return
      end