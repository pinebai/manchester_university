!
!=====TRANSPORT AND THERMODYNAMIC PROPERTIES
!
      subroutine properties
      include "common.inc"
      
      do i=1,ia
        do j=1,ja
          do k=1,ka
!
!-----GAS PROPERTIES
!
            dvsg(i,j,k)=0.00002
            dng(i,j,k)=1.177
            tg(i,j,k)=300.
            hcpg(i,j,k)=1004.9
            hcvg(i,j,k)=717.8
            p(i,j,k)=dng(i,j,k)*(hcpg(i,j,k)-hcvg(i,j,k))*tg(i,j,k)
!
!-----LIQUID PROPERTIES
!
            dvsl(i,j,k)=0.0027
            dnl(i,j,k)=814.
            st(i,j,k)=0.0283
          enddo
        enddo
      enddo

      return
      end