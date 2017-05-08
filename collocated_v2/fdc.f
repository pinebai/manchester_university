!
!=====CONVECTION - DIFFUSION COEFFICIENTS
!
      subroutine fdc(uvel,vvel)
      include "common.inc"
      
      real uvel(ia,ja),vvel(ia,ja)
      
      do i=imn,imx
        do j=jmn,jmx
!
!-----CONVECTION TERMS  F = RHO * A * ALPHA * U
!
          fw(i,j)=dasw(i,j)
     &            *0.5*(den(i-1,j)*qq(i-1,j)*uvel(i-1,j)
     &            +den(i,j)*qq(i,j)*uvel(i,j))
          fe(i,j)=dase(i,j)
     &            *0.5*(den(i,j)*qq(i,j)*uvel(i,j)
     &            +den(i+1,j)*qq(i+1,j)*uvel(i+1,j))
          fs(i,j)=dass(i,j)
     &            *0.5*(den(i,j-1)*qq(i,j-1)*vvel(i,j-1)
     &            +den(i,j)*qq(i,j)*vvel(i,j))
          fn(i,j)=dasn(i,j)
     &            *0.5*(den(i,j)*qq(i,j)*vvel(i,j)
     &            +den(i,j+1)*qq(i,j+1)*vvel(i,j+1))
!
!-----DIFFUSION TERMS  D = MU * A * ALPHA * DPHI / DX
!
          dw(i,j)=dasw(i,j)
     &            *0.5*(dvis(i-1,j)*qq(i-1,j)
     &            +dvis(i,j)*qq(i,j))
     &            /dzs(i)
          de(i,j)=dase(i,j)
     &            *0.5*(dvis(i,j)*qq(i,j)
     &            +dvis(i+1,j)*qq(i+1,j))
     &            /dzs(i)
          ds(i,j)=dass(i,j)
     &            *0.5*(dvis(i,j-1)*qq(i,j-1)
     &            +dvis(i,j)*qq(i,j))
     &            /drs(j)
          dn(i,j)=dasn(i,j)
     &            *0.5*(dvis(i,j)*qq(i,j)
     &            +dvis(i,j+1)*qq(i,j+1))
     &            /drs(j)
        enddo
      enddo
      
      return
      end