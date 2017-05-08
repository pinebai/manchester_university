!
!=====U - COMPONENT GAS VELOCITY
!
      subroutine cug
      include "common.inc"
      
      ivel=1
      iph=1
      iq=3
      iap=1
      call zro2
      call var
      call fdc(ug,vg)
      call hyb
      do ib=1,4
        call ncs(ib,ug)
      enddo
      call ics(ug)
      call idg(ug,ul3)
      call mcf(ug0)
      call sol1(ug)

      return
      end