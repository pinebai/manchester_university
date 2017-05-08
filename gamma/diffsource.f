!
!=====VISCOUS STRESS SOURCE
!
      subroutine vstress(uvel,vvel,wvel)
      include "common.inc"
      real uvel(ia,ja,ka),vvel(ia,ja,ka),wvel(ia,ja,ka)
      
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
!
!-----U - COMPONENT
!
            if(ivel.eq.1)then
!
!-----TERM 1
!
              dudxw=(dd(i,j,k)*uvel(i,j,k)-dd(i-1,j,k)*uvel(i-1,j,k))
     &          /(dxc(i)-dxc(i-1))
              dudxe=(dd(i+1,j,k)*uvel(i+1,j,k)-dd(i,j,k)*uvel(i,j,k))
     &          /(dxc(i+1)-dxc(i))
              
              du2d2x=(dudxe-dudxw)/dx(i)
!
!-----TERM 2
!
              vsfw=f2(3,i-1,j,k,dd,vvel)
              vsfp=f2(3,i,j,k,dd,vvel)
              vsfe=f2(3,i+1,j,k,dd,vvel)

              vsw=fxw(i)*vsfw+(1-fxw(i))*vsfp
              vse=(1-fxe(i))*vsfp+fxe(i)*vsfe
              dvdxs=(vse-vsw)/dx(i)

              vnfw=f2(4,i-1,j,k,dd,vvel)
              vnfp=f2(4,i,j,k,dd,vvel)
              vnfe=f2(4,i+1,j,dd,k,vvel)

              vnw=fxw(i)*vnfw+(1-fxw(i))*vnfp
              vne=(1-fxe(i))*vnfp+fxe(i)*vnfe
              dvdxn=(vne-vnw)/dx(i)

              dv2dydx=(dvdxn-dvdxs)/dy(j)
!
!-----TERM 3
!
              wbfw=f2(5,i-1,j,k,dd,wvel)
              wbfp=f2(5,i,j,k,dd,wvel)
              wbfe=f2(5,i+1,j,dd,k,wvel)

              wbw=fxw(i)*wbfw+(1-fxw(i))*wbfp
              wbe=(1-fxe(i))*wbfp+fxe(i)*wbfe
              dwdxb=(wbe-wbw)/dx(i)

              wtfw=f2(6,i-1,j,k,dd,wvel)
              wtfp=f2(6,i,j,k,dd,wvel)
              wtfe=f2(6,i+1,j,k,dd,wvel)

              wtw=fxw(i)*wtfw+(1-fxw(i))*wtfp
              wte=(1-fxe(i))*wtfp+fxe(i)*wtfe
              dwdxt=(wte-wtw)/dx(i)

              dw2dzdx=(dwdxt-dwdxb)/dz(k)
!
!-----SOURCE
!
              sum=du2d2x+dv2dydx+dw2dzdx
!
!-----V - COMPONENT
!
            elseif(ivel.eq.2)then

          enddo
        enddo
      enddo
 
      return
      end