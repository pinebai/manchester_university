!
!=====SECOND MOMENT V - LIQUID VELOCITY
!
      subroutine cvl2

      include "common.inc"
      
      ivel=2
      iph=2
      iq=2
      iap=1
      call zro2
      call var
      call fdc(ul2,vl2)
      call zro3
      call hyb
      do ib=1,4
        call ncs(ib,vl2)
      enddo
      call ics(vl2)
      call idg(vg,vl2)
      call mcf(vl20)
      call sol1(vl2)

      return
      end