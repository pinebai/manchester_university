!----------------------------------------------------------
!
!     DROPLET-SIZE MOMENTS MODEL
!     USING LAGUERRE DISTRIBUTION CLOSURE MODEL.
!
!     CONTAINING:
!     STANDARD AND MODIFIED DRAG MODELS.
!
!     EMPLOYING:
!     FULLY COLLOCATED VELOCITY SCHEME.
!
!     DOMINIC P. JONES  13.10.2006
!
!----------------------------------------------------------
!
      program main
      real resliq(100,100),resgas(100,100)
      include "common.inc"
!
!-----monitor residuals
!
      open(9,file='residuals.txt')
15    format(2(i5),2(e15.5))
!
!-----INITIALISE
!
      call intialise
      do nts=1,ntsmx
        time=nts*delt
!
!-----COMMENCE TIME STEP
!
        call store
        do ito=1,itomx
!
!-----LIQUID PHASE
!
          iph=2
          call calcq0
          call calcq1
          call calcq2
          if(iqmax.ge.4)call calcq3
          if(iqmax.ge.5)call calcq4
          if(iqmax.ge.6)call calcq5
          call calcqd
          call calcuq0
          call calcvq0
          call calcwq0
          call calcuq1
          call calcvq1
          call calcwq1
          call calcuq2
          call calcvq2
          call calcwq2
          if(iqmax.ge.4)then
            call calcuq3
            call calcvq3
            call calcwq3
          endif
          if(iqmax.ge.5)then
            call calcuq4
            call calcvq4
            call calcwq4
          endif
          if(iqmax.ge.6)then
            call calcuq5
            call calcvq5
            call calcwq5
          endif
          if(iqmax.ge.7)then
            call calcuq6
            call calcvq6
            call calcwq6
          endif
          if(iqmax.ge.8)then
            call calcuq7
            call calcvq7
            call calcwq7
          endif
!           call converge
!           resliq(nts,ito)=sliq
!
!-----GAS PHASE
!
          iph=1
          call calcug
          call calcvg
          call calcwg
          call calcp
!           call converge
!           resgas(nts,ito)=sgas
!
!-----GLOBAL CONVERGENCE
!
          print*,'TIME,NTS,ITO',time,nts,ito
!           iph=0
!           call converge
!           write(9,15)nts,ito,resliq(nts,ito),resgas(nts,ito)
!           if(gsmax.lt.sormx.and.ito.ge.1)goto 300
        enddo
! 300     continue
!
!-----OUTPUT
!
        call output1
!         call output2
        call penetration
      enddo
!
      continue
      close(9)
      stop
      end