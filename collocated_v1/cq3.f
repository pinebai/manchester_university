!
!=====THIRD MOMENT
!
      subroutine cq3
      include "common.inc"
      
      ivel=0
      iph=2
      iq=3
      iap=1
      call zro2
      call var
      call fdc(ul3,vl3)
      call zro3
      call hyb
      do ib=1,4
        call ncs(ib,q3)
      enddo
      call ics(q3)
      call mcf(q30)
      call sol1(q3)

      return
      end