!
!=====DESCRIBE GRID
!
      subroutine grd
      include "common.inc"
!
!-----POSITION OF CELL CENTERS
!
      alpha=1.
101   continue
      i=imin
      zc(i)=0.
      do i=imin+1,iinj+ninjc
        zc(i)=zc(i-1)+dzc
      enddo
      do i=imn+ninjc+1,imax
        zc(i)=zc(i-1)
     &        +dzc*(1+i/real(imax))**alpha
      enddo
      if(zc(imax).lt.zl)then
        alpha=alpha+0.1
        goto 101
      endif
      
      alpha=1.
202   continue
      j=jmin
      rc(j)=0.
      do j=jmin+1,jinj+ninjc
        rc(j)=rc(j-1)+drc
      enddo
      do j=jmn+ninjc+1,jmax
        rc(j)=rc(j-1)
     &        +drc*(1+j/real(jmax))**alpha
      enddo
      if(rc(jmax).lt.rl)then
        alpha=alpha+0.1
        goto 202
      endif
!
!-----CONTROL SURFACE LENGTHS
!
      do i=imin+1,imax-1
        dzs(i)=0.5*(zc(i+1)+zc(i))
     &         -0.5*(zc(i)+zc(i-1))
      enddo
      dzs(imin)=dzs(imin+1)
      dzs(imax)=dzs(imax-1)
      
      do j=jmin+1,jmax-1
        drs(j)=0.5*(rc(j+1)+rc(j))
     &         -0.5*(rc(j)+rc(j-1))
      enddo
      drs(jmin)=drs(jmin+1)
      drs(jmax)=drs(jmax-1)
!
!-----CONTROL SURFACE AREAS
!
      do i=imin+1,imax-1
        do j=jmin+1,jmax-1
          dasw(i,j)=0.5*((rc(j)+0.5*drs(j))**2-(rc(j)-0.5*drs(j))**2)
          dase(i,j)=0.5*((rc(j)+0.5*drs(j))**2-(rc(j)-0.5*drs(j))**2)
          dass(i,j)=(rc(j)-0.5*drs(j))
     &               *((zc(i)+dzs(i))**2-(zc(i)-dzs(i))**2)
          dasn(i,j)=(rc(j)+0.5*drs(j))
     &               *((zc(i)+dzs(i))**2-(zc(i)-dzs(i))**2)
        enddo
      enddo
!
!-----CONTROL SURFACE VOLUMES
!
      do i=imin+1,imax-1
        do j=jmin+1,jmax-1
          dvs(i,j)=0.5*((rc(j)+0.5*drs(j))**2-(rc(j)-0.5*drs(j))**2)
     &             *((zc(i)+dzs(i))-(zc(i)-dzs(i)))
        enddo
      enddo
      
      return
      end