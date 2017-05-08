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
          fw(i,j)=areaw(i,j)
     &            *0.5*(den(i-1,j)*qq(i-1,j)*uvel(i-1,j)
     &            +den(i,j)*qq(i,j)*uvel(i,j))
          fe(i,j)=areae(i,j)
     &            *0.5*(den(i,j)*qq(i,j)*uvel(i,j)
     &            +den(i+1,j)*qq(i+1,j)*uvel(i+1,j))
          fs(i,j)=areas(i,j)
     &            *0.5*(den(i,j-1)*qq(i,j-1)*vvel(i,j-1)
     &            +den(i,j)*qq(i,j)*vvel(i,j))
          fn(i,j)=arean(i,j)
     &            *0.5*(den(i,j)*qq(i,j)*vvel(i,j)
     &            +den(i,j+1)*qq(i,j+1)*vvel(i,j+1))
!
!-----DIFFUSION TERMS  D = MU * A * ALPHA * DPHI / DX
!
          dw(i,j)=areaw(i,j)
     &            *0.5*(dvs(i-1,j)*qq(i-1,j)
     &            +dvs(i,j)*qq(i,j))
     &            /dxs(i)
          de(i,j)=areae(i,j)
     &            *0.5*(dvs(i,j)*qq(i,j)
     &            +dvs(i+1,j)*qq(i+1,j))
     &            /dxs(i)
          ds(i,j)=areas(i,j)
     &            *0.5*(dvs(i,j-1)*qq(i,j-1)
     &            +dvs(i,j)*qq(i,j))
     &            /dys(j)
          dn(i,j)=arean(i,j)
     &            *0.5*(dvs(i,j)*qq(i,j)
     &            +dvs(i,j+1)*qq(i,j+1))
     &            /dys(j)
        enddo
      enddo
      
      return
      end