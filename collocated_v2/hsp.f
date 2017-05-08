!
!=====DROPLET-SIZE MOMENTS MODEL
!
      program hsp
      include "common.inc"
!
!-----INITIALISE
!
      call int
!
!-----TIME STEP
!
      do nts=0,ntsmx
!
!-----STORE TIME STEP DATA
!
        call pts
        print*,'time',nts/real(ntsmx)
        npt=npt+1
!
!-----ITERATION CYCLE
!
C         do itc=0,itcmx
!
!-----GAS PRESSURE
!
C           call cpr(2)
!
!-----GAS VELOCITIES
!
          call cug
          call cvg
!
!-----PRESSURE CORRECTION
!
C           call cpr(1)
!
!-----VELOCITY CORRECTION
!
          call cuv
!
!-----CONVERGENCE
!
C         enddo
!
!-----LIQUID MOMENTS
!
        call cq1
        call cq2
        call cq3
        call cqk
!
!-----LIQUID VELOCITIES
!
        call cul1
        call cvl1
        call cul2
        call cvl2
        call cul3
        call cvl3
!
!-----CONVERGENCE
!
C         enddo
!
!-----OUTPUT
!
        if(npt.gt.9)then
          print*,'- output'
          call opt
          npt=0
        endif
!
!-----TIME STEP
!
      enddo
!
!-----END PROGRAM
!
      stop
      end