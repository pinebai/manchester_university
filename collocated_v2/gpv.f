!
!-----SOLVE GAS PSEUDO-VELOCITIES
!
      subroutine gpv
      include "common.inc"
      
      do i=imn,imx
        do j=jmn,jmx
          ugxx(i,j)=(awu(i,j)*ug(i-1,j)+aeu(i+1,j)*ug(i+1,j)
     &              +asu(i,j)*ug(i,j-1)+anu(i,j)*ug(i,j+1)+suu(i,j))
     &              /(apu(i,j)+tiny)
          vgxx(i,j)=(awv(i,j)*vg(i-1,j)+aev(i+1,j)*vg(i+1,j)
     &              +asv(i,j)*vg(i,j-1)+anv(i,j)*vg(i,j+1)+suv(i,j))
     &              /(apv(i,j)+tiny)
        enddo
      enddo

      return
      end