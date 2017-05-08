!
!=====GRID
!
      subroutine grid
      include "common.inc"
!
!================================================
!     CV FACE LOCATIONS DEFINED BY USER
!================================================
!
!=====COORDINATES OF CV-FACES
!
!
!-----i driection
!
      alpha=0
100   continue
      x(imn-1)=0.
      x(imn)=xinj
      do i=imn+1,imx
        x(i)=x(i-1)+xinj
     &       *(1.+(i-imn+1)/real(imx-imn+1))**alpha
      enddo
      if(x(imx).lt.xl)then
        alpha=alpha+0.01
        goto 100
      endif
!
!-----j driection
!
      alpha=0
200   continue
      y(jinj)=0.
      do j=jinj+1,jinj+nradic
        y(j)=y(j-1)+dinj
      enddo
      do j=jinj+nradic+1,jmx+1
        y(j)=y(j-1)+dinj
     &       *(1.+(j-jinj-1)/real(jmx-jinj-1))**alpha
      enddo
      if(y(jmx).lt.yl)then
        alpha=alpha+0.01
        goto 200
      endif
      do j=jinj-1,jmn-1,-1
        y(j)=-y(2*jinj-j)
      enddo
!
!-----k driection
!
      z(kinj)=0.
      do k=kinj+1,kinj+nradic
        z(k)=z(k-1)+dinj
      enddo
      do k=kinj+nradic+1,kmx+1
        z(k)=z(k-1)+dinj
     &       *(1.+(k-kinj-1)/real(kmx-kinj-1))**alpha
      enddo
      do k=kinj-1,kmn-1,-1
        z(k)=-z(2*kinj-k)
      enddo
!
!================================================
!================================================
!
!-----COORDINATES OF CV-CENTERS
!
      xc(imn-1)=x(imn-1)
      do i=imn,imx
        xc(i)=(x(i)+x(i-1))/2
      end do
      xc(imx+1)=x(imx)
!
      yc(jmn-1)=y(jmn-1)
      do j=jmn,jmx
        yc(j)=(y(j)+y(j-1))/2
      end do
      yc(jmx+1)=y(jmx)
!
      zc(kmn-1)=z(kmn-1)
      do k=kmn,kmx
        zc(k)=(z(k)+z(k-1))/2
      end do
      zc(kmx+1)=z(kmx)
!
!-----INTERPOLATION FACTORS
!
      do i=imn,imx
        fxw(i)=1-(x(i-1)-xc(i-1))/(xc(i)-xc(i-1))
      end do
      do i=imn,imx
        fxe(i)=1-(xc(i+1)-x(i))/(xc(i+1)-xc(i))
      end do
!
      do j=jmn,jmx
        fys(j)=1-(y(j-1)-yc(j-1))/(yc(j)-yc(j-1))
      end do
      do j=jmn,jmx
        fyn(j)=1-(yc(j+1)-y(j))/(yc(j+1)-yc(j))
      end do
!
      do k=kmn,kmx
        fzb(k)=1-(z(k-1)-zc(k-1))/(zc(k)-zc(k-1))
      end do
      do k=kmn,kmx
        fzt(k)=1-(zc(k+1)-z(k))/(zc(k+1)-zc(k))
      end do
!
!-----LENGTHS
!
      do i=imn,imx
        dx(i)=x(i)-x(i-1)
      enddo
!
      do j=jmn,jmx
        dy(j)=y(j)-y(j-1)
      enddo
!
      do k=kmn,kmx
        dz(k)=z(k)-z(k-1)
      enddo
      
!
!-----AREAS
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            area(1,i,j,k)=dy(j)*dz(k)
            area(2,i,j,k)=dy(j)*dz(k)
            area(3,i,j,k)=dx(i)*dz(k)
            area(4,i,j,k)=dx(i)*dz(k)
            area(5,i,j,k)=dx(i)*dy(j)
            area(6,i,j,k)=dx(i)*dy(j)
          enddo
        enddo
      enddo
!
!-----VOLUMES
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            vol(i,j,k)=dx(i)*dy(j)*dz(k)
          enddo
        enddo
      enddo

      return
      end