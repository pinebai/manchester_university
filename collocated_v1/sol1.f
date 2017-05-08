!
!=====TDMA SOLVER
!
      subroutine sol1(phi)
      include "common.inc"

      real phi(ia,ja),c(ia,ja),aa(ia,ja),cc(ia,ja)

      do m=1,3
        do j=jmn,jmx
          i=imn
          c(i,j)=as(i,j)*phi(i,j-1)+an(i,j)*phi(i,j+1)+su(i,j)
     &           +aw(i,j)*phi(i-1,j)
          aa(i,j)=ae(i,j)/(ap(i,j)+tiny)
          cc(i,j)=c(i,j)/(ap(i,j)+tiny)
          do i=imn+1,imx-1
            c(i,j)=as(i,j)*phi(i,j-1)+an(i,j)*phi(i,j+1)+su(i,j)
            aa(i,j)=ae(i,j)/(ap(i,j)-aw(i,j)*aa(i-1,j)+tiny)
            cc(i,j)=(aw(i,j)*cc(i-1,j)+c(i,j))
     &              /(ap(i,j)-aw(i,j)*aa(i-1,j)+tiny)
          enddo
          i=imx
          c(i,j)=as(i,j)*phi(i,j-1)+an(i,j)*phi(i,j+1)+su(i,j)
     &           +ae(i,j)*phi(i+1,j)
          aa(i,j)=0.
          cc(i,j)=(aw(i,j)*cc(i-1,j)+c(i,j))
     &            /(ap(i,j)-aw(i,j)*aa(i-1,j)+tiny)
        enddo
        do j=jmn,jmx
          do i=imx,imn,-1
            phi(i,j)=phi(i+1,j)*aa(i,j)+cc(i,j)
          enddo
        enddo
      enddo

      return
      end