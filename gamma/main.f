!
!=====DROPLET-SIZE MOMENTS MODEL
!     USING GAMMA DISTRIBUTION CLOSURE MODEL
!
!     DOMINIC JONES
!     19.09.2006
!
      program main
      include "common.inc"
!
!-----INITIALISE
!
      call intialise
      do nts=1,ntsmx
        time=nts*delt
!
!-----LIQUID PHASE
!
        call store
        do ito=1,itomx
            call calcq1
            call calcq2
            call calcq3
            call calcul1
            call calcul2
            call calcul3
            call calcvl1
            call calcvl2
            call calcvl3
            call calcwl1
            call calcwl2
            call calcwl3
            call converge(2)
!
!-----GAS PHASE
!
            call calcug
            call calcvg
            call calcwg
            call calcp
            call converge(1)
!
!-----GLOBAL CONVERGENCE
!
          call converge(0)
          if(gdelta.lt.sormx.and.ito.gt.10)goto 300
        enddo
300     continue
!
!-----OUTPUT
!
        call output1
      enddo
!
      continue
      stop
      end