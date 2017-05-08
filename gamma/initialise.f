!
!=====INITIALISE
!
      subroutine intialise
      include "common.inc"
!
!-----SOLVED CONSTANTS
!
      pi=4.*atan(1.)
      cnt=(4/3.)*pi
!
!-----ARRAY INDEX
!
      imn=2
      jmn=2
      kmn=2
      imx=ia-1
      jmx=ja-1
      kmx=ka-1
!
!-----MONITOR INDEX
!
      imon=imn+5
      jmon=(jmx+jmn)/2
      kmon=(kmx+kmn)/2
!
!-----INJECTOR INDEX REFERENCE
!
      iinj=imn+1
      jinj=(jmx+jmn)/2
      kinj=(kmx+kmn)/2
!
!-----ZERO MAIN ARRAYS
!
      call zeromain
!
!-----THERMODYAMIC PROPERTIES
!
      call properties
!
!-----INJECTION MODEL
!
      call injection
      
      return
      end