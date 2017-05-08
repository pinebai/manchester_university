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
            dng(i,j,k)=13.2
            tg(i,j,k)=300.
            hcpg(i,j,k)=1004.9
            hcvg(i,j,k)=717.8
            p(i,j,k)=dng(i,j,k)*(hcpg(i,j,k)-hcvg(i,j,k))*tg(i,j,k)
!
!-----LIQUID PROPERTIES
!
            dvsl(i,j,k)=0.0027
            dnl(i,j,k)=840.
            st(i,j,k)=0.0283
          enddo
        enddo
      enddo

      return
      end
!
!=====CHARACTERISTIC NUMBERS
!
      function we(vrel,rad,i,j,k)
      include "common.inc"
      
      we=dng(i,j,k)*vrel**2*rad/st(i,j,k)
      
      return
      end
!
      function re(vrel,rad,i,j,k)
      include "common.inc"
      
      re=2*dng(i,j,k)*vrel*rad/dvsg(i,j,k)
      
      return
      end
!
      function oh(rad,i,j,k)
      include "common.inc"
      
      oh=dvsl(i,j,k)/(dnl(i,j,k)*rad*st(i,j,k))**0.5
      
      return
      end
!
      function wecrit(vrel,i,j,k)
      include "common.inc"
      
      wecrit=12*(1.+1.077*(oh(vrel,i,j,k)/7.))
      
      return
      end