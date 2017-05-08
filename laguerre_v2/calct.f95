!
!=====DROPLET TEMPERATURE
!
      subroutine calctl
      include "common.inc"
!
      ivel=0;iph=2;iq=3
!
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=dnl0(i,j,k)*(1.-(4/3.)*pi*q3(i,j,k))
            ff(i,j,k)=dnl(i,j,k)*(1.-(4/3.)*pi*q3(i,j,k))
            dd(i,j,k)=dvsl(i,j,k)*(1.-(4/3.)*pi*q3(i,j,k))
          enddo
        enddo
      enddo
!
      call zerocoef
      call convection(uq3,vq3,wq3)
      call hybrid
      if(iquick.eq.1)call quick(tl)
      do ib=2,6
        call nbc(ib,tl)
      enddo
      call dbc(1,tl,0.)
      call injbc(tl)
      call spbc
!       call temporal(tl0)
      call solve(tl,tlin)
!
      return
      end