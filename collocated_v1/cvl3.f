!
!=====THIRD MOMENT V - LIQUID VELOCITY
!
      subroutine cvl3

      include "common.inc"
      
      ivel=2
      iph=2
      iq=3
      iap=1
      call zro2
      call var
      call fdc(ul3,vl3)
      call zro3
      call hyb
      do ib=1,4
        call ncs(ib,vl3)
      enddo
      call ics(vl3)
      call idg(vg,vl3)
      call mcf(vl30)
      call sol1(vl3)

      return
      end