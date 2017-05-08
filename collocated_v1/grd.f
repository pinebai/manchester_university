!
!=====DESCRIBE GRID
!
      subroutine grd
      include "common.inc"
!
!-----POSITION OF CELL CENTERS
!
      alpha=1.
100   continue
      i=imin
      xc(i)=0.
      do i=imin+1,iinj+ninjc
        xc(i)=xc(i-1)+dxc
      enddo
      do i=imn+ninjc+1,imax
        xc(i)=xc(i-1)
     &        +dxc*(1+i/real(imax))**alpha
      enddo
      if(xc(imax).lt.xl)then
        alpha=alpha+0.1
        goto 100
      endif
      
      alpha=1.
200   continue
      j=jmin
      yc(j)=0.
      do j=jmin+1,jinj+ninjc
        yc(j)=yc(j-1)+dyc
      enddo
      do j=jmn+ninjc+1,jmax
        yc(j)=yc(j-1)
     &        +dyc*(1+j/real(jmax))**alpha
      enddo
      if(yc(jmax).lt.yl)then
        alpha=alpha+0.1
        goto 200
      endif
!
!-----CONTROL SURFACE LENGTHS
!
      do i=imin+1,imax-1
        dxs(i)=0.5*(xc(i+1)+xc(i))
     &         -0.5*(xc(i)+xc(i-1))
      enddo
      dxs(imin)=dxs(imin+1)
      dxs(imax)=dxs(imax-1)
      
      do j=jmin+1,jmax-1
        dys(j)=0.5*(yc(j+1)+yc(j))
     &         -0.5*(yc(j)+yc(j-1))
      enddo
      dys(jmin)=dys(jmin+1)
      dys(jmax)=dys(jmax-1)
!
!-----CONTROL SURFACE AREAS
!
      do i=imin+1,imax-1
        do j=jmin+1,jmax-1
          areaw(i,j)=0.5*((yc(j)+0.5*dys(j))**2-(yc(j)-0.5*dys(j))**2)
          areae(i,j)=0.5*((yc(j)+0.5*dys(j))**2-(yc(j)-0.5*dys(j))**2)
          areas(i,j)=(yc(j)-0.5*dys(j))
     &               *((xc(i)+dxs(i))**2-(xc(i)-dxs(i))**2)
          arean(i,j)=(yc(j)+0.5*dys(j))
     &               *((xc(i)+dxs(i))**2-(xc(i)-dxs(i))**2)
        enddo
      enddo
!
!-----CONTROL SURFACE VOLUMES
!
      do i=imin+1,imax-1
        do j=jmin+1,jmax-1
          vol(i,j)=0.5*((yc(j)+0.5*dys(j))**2-(yc(j)-0.5*dys(j))**2)
     &             *((xc(i)+dxs(i))-(xc(i)-dxs(i)))
        enddo
      enddo
      
      return
      end