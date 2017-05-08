!
!=====THIRD MOMENT U - LIQUID VELOCITY
!
      subroutine cul3

      include "common.inc"
      
      ivel=1
      iph=2
      iq=3
      iap=1
      call zro2
      call var
      call fdc(ul3,vl3)
      call zro3
      call hyb
      do ib=1,4
        call ncs(ib,ul3)
      enddo
      call ics(ul3)
      call idg(ug,ul3)
      call mcf(ul30)
      call sol1(ul3)

      return
      end