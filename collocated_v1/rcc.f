!
!-----CONITINUTY IMBALANCE TERM USING PSEUDO-VELOCITIES
!
      subroutine rcc(uvel,vvel)
      include "common.inc"

      real uvel(ia,ja),vvel(ia,ja)
      
      do i=imn,imx
        do j=jmn,jmx
          su(i,j)=-(0.5*uvel(i+1,j)*dng(i+1,j)
     &            *areae(i,j)*(1-cnt*q3(i+1,j))
     &            +0.25*du(i,j)*dng(i,j)*areae(i,j)*(1-cnt*q3(i,j))
     &            *(p(i+2,j)-3*p(i+1,j)+3*p(i,j)-p(i-1,j)))
     &            +(0.5*uvel(i-1,j)*dng(i-1,j)
     &            *areaw(i,j)*(1-cnt*q3(i-1,j))
     &            +0.25*du(i,j)*dng(i,j)*areaw(i,j)*(1-cnt*q3(i,j))
     &            *(p(i+1,j)-3*p(i,j)+3*p(i-1,j)-p(i-2,j)))
     &            -(0.5*vvel(i,j+1)*dng(i,j+1)
     &            *arean(i,j)*(1-cnt*q3(i,j+1))
     &            +0.25*dv(i,j)*dng(i,j)*arean(i,j)*(1-cnt*q3(i,j))
     &            *(p(i,j+2)-3*p(i,j+1)+3*p(i,j)-p(i,j-1)))
     &            +(0.5*vvel(i,j-1)*dng(i,j-1)
     &            *areas(i,j)*(1-cnt*q3(i,j-1))
     &            +0.25*dv(i,j)*dng(i,j)*areas(i,j)*(1-cnt*q3(i,j))
     &            *(p(i,j+1)-3*p(i,j)+3*p(i,j-1)-p(i,j-2)))
     &            +(dng0(i,j)*qq0(i,j)-dng(i,j)*qq(i,j))
     &            *vol(i,j)/delt
        enddo
      enddo

      return
      end