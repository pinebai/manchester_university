!
!=====POWER-LAW SCHEME
!
      subroutine pls
      
      include "common.inc"

      do i=imn,imx
        do j=jmn,jmx
          aw(i,j)=dw(i,j)
     &            *amax1(0.,(1-0.1*abs(fw(i,j)/(dw(i,j)+tiny)))**5)
     &            +amax1(0.,fw(i,j))
          ae(i,j)=de(i,j)
     &            *amax1(0.,(1-0.1*abs(fe(i,j)/(de(i,j)+tiny)))**5)
     &            +amax1(0.,-fe(i,j))
          as(i,j)=ds(i,j)
     &            *amax1(0.,(1-0.1*abs(fs(i,j)/(ds(i,j)+tiny)))**5)
     &            +amax1(0.,fs(i,j))
          an(i,j)=dn(i,j)
     &            *amax1(0.,(1-0.1*abs(fn(i,j)/(dn(i,j)+tiny)))**5)
     &            +amax1(0.,-fn(i,j))
        enddo
      enddo

      return
      end