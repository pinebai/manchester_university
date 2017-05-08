!
!=====NEUMANN CONDITIONS
!
      subroutine ncs(ib,phi)
      include "common.inc"
      
      real phi(ia,ja)
!
!-----WEST
!
      if(ib.eq.1)then
        i=imn
        do j=jmn,jmx
          phi(i-1,j)=phi(i,j)
          aw(i,j)=0.
        enddo
!
!-----EAST
!
      elseif(ib.eq.2)then
        i=imx
        do j=jmn,jmx
          phi(i+1,j)=phi(i,j)
          ae(i,j)=0.
        enddo
!
!-----SOUTH
!
      elseif(ib.eq.3)then
        j=jmn
        do i=imn,imx
          phi(i,j-1)=phi(i,j)
          as(i,j)=0.
        enddo
!
!-----NORTH
!
      elseif(ib.eq.4)then
        j=jmx
        do i=imn,imx
          phi(i,j+1)=phi(i,j)
          an(i,j)=0.
        enddo
      endif

      return
      end