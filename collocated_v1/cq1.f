!
!=====FIRST MOMENT
!
      subroutine cq1
      include "common.inc"
      
      ivel=0
      iph=2
      iq=1
      iap=1
      call zro2
      call var
      call fdc(ul1,vl1)
      call zro3
      call hyb
      do ib=1,4
        call ncs(ib,q1)
      enddo
      call ics(q1)
      call mcf(q10)
      call sol1(q1)

      return
      end