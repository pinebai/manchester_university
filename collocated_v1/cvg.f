!
!=====V - COMPONENT GAS VELOCITY
!
      subroutine cvg
      include "common.inc"
      
      ivel=2
      iph=1
      iq=3
      iap=1
      call zro2
      call var
      call fdc(ug,vg)
      call hyb
      do ib=1,4
        call ncs(ib,vg)
      enddo
      call ics(vg)
      call idg(vg,vl3)
      call mcf(vg0)
      call sol1(vg)

      return
      end