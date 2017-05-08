!
!=====SECOND MOMENT U - LIQUID VELOCITY
!
      subroutine cul2

      include "common.inc"
      
      ivel=1
      iph=2
      iq=2
      iap=1
      call zro2
      call var
      call fdc(ul2,vl2)
      call zro3
      call hyb
      do ib=1,4
        call ncs(ib,ul2)
      enddo
      call ics(ul2)
      call idg(ug,ul2)
      call mcf(ul20)
      call sol1(ul2)

      return
      end