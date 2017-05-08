!
!=====TRANSPORT AND THERMODYNAMIC PROPERTIES
!
      subroutine prp
      include "common.inc"
      
      do i=imin,imax
        do j=jmin,jmax
!
!-----GAS PROPERTIES
!
          dvsg(i,j)=0.00002
          dng(i,j)=1.177
          tg(i,j)=300.
          hcpg(i,j)=1004.9
          hcvg(i,j)=717.8
          p(i,j)=dng(i,j)*(hcpg(i,j)-hcvg(i,j))*tg(i,j)
!
!-----LIQUID PROPERTIES
!
          dvsl(i,j)=0.0027
          dnl(i,j)=814.
        enddo
      enddo

      return
      end