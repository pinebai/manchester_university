!
!=====FIRST MOMENT U - LIQUID VELOCITY
!
      subroutine cul1

      include "common.inc"
      
      ivel=1
      iph=2
      iq=1
      iap=1
      call zro2
      call var
      call fdc(ul1,vl1)
      call zro3
      call hyb
      do ib=1,4
        call ncs(ib,ul1)
      enddo
      call ics(ul1)
      call idg(ug,ul1)
      call mcf(ul10)
      call sol1(ul1)

      return
      end