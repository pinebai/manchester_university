!
!=====FIRST MOMENT V - LIQUID VELOCITY
!
      subroutine cvl1

      include "common.inc"
      
      ivel=2
      iph=2
      iq=1
      iap=1
      call zro2
      call var
      call fdc(ul1,vl1)
      call zro3
      call hyb
      do ib=1,4
        call ncs(ib,vl1)
      enddo
      call ics(vl1)
      call idg(vg,vl1)
      call mcf(vl10)
      call sol1(vl1)

      return
      end