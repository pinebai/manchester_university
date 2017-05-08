!
!=====INITIALISE
!
      subroutine int
      include "common.inc"
!
!-----ARRAY INDEX
!
      imin=1
      jmin=1
      imax=ia
      jmax=ja
      call idx
!
!-----CONSTANTS
!
      npt=0
      tiny=1.e-30
      great=1.e+30
      pi=4.*atan(1.)
      cnt=(4/3.)*pi
!
!-----PRESSURE-VELOCITY RELAXATION
!
      ur=0.2
!
!-----NUMBER OF ITERATIONS
!
      itcmx=10
!
!-----PARAMETER K CUT-OFF VALUES
!
      pkmin=1.1
      pkmax=10.
!
!-----DOMAIN SIZE
!
      xl=0.05
      yl=0.01
!
!-----ZERO MAIN ARRAYS
!
      call zro1
!
!-----THERMODYAMIC PROPERTIES
!
      call prp
!
!-----INJECTION MODEL
!
      call inj
      
      return
      end