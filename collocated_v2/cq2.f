!
!=====SECOND MOMENT
!
      subroutine cq2
      include "common.inc"
      
      ivel=0
      iph=2
      iq=2
      iap=1
      call zro2
      call var
      call fdc(ul2,vl2)
      call zro3
      call hyb
      do ib=1,4
        call ncs(ib,q2)
      enddo
      call ics(q2)
      call mcf(q20)
      call sol1(q2)

      return
      end