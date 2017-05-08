!
!=====DEFINE VARIABLES
!
      subroutine var
      include "common.inc"
      
      do i=imin,imax
        do j=jmin,jmax
!
!-----GAS PHASE DENSITY AND VSICOSITY
!
          if(iph.eq.1)then
            den0(i,j)=dng0(i,j)
            den(i,j)=dng(i,j)
            dvis(i,j)=dvsg(i,j)
!
!-----LIQUID PHASE DENSITY AND VSICOSITY
!
          elseif(iph.eq.2)then
            den0(i,j)=dnl0(i,j)
            den(i,j)=dnl(i,j)
            dvis(i,j)=dvsl(i,j)
          endif
!
!-----Q1 MOMENT (FOR MOMENTUM EQU)
!
          if(iq.eq.1.and.ivel.ge.1)then
            qq0(i,j)=q10(i,j)
            qq(i,j)=q1(i,j)
!
!-----Q2 MOMENT (FOR MOMENTUM EQU)
!
          elseif(iq.eq.2.and.ivel.ge.1)then
            qq0(i,j)=q20(i,j)
            qq(i,j)=q2(i,j)
!
!-----GAS VOLUME FRACTION (FOR MOMENTUM EQU)
!
          elseif(iq.eq.3.and.ivel.ge.1.and.iph.eq.1)then
            qq0(i,j)=1.-cnt*q30(i,j)
            qq(i,j)=1.-cnt*q3(i,j)
!
!-----LIQUID VOLUME FRACTION (FOR MOMENTUM EQU)
!
          elseif(iq.eq.3.and.ivel.ge.1.and.iph.eq.2)then
            qq0(i,j)=cnt*q30(i,j)
            qq(i,j)=cnt*q3(i,j)
!
!-----Q1 MOMENT COEFFICIENT (FOR MOMENT EQU)
!
          elseif(iq.eq.1.and.ivel.eq.0)then
            qq0(i,j)=1.
            qq(i,j)=1.
!
!-----Q2 MOMENT COEFFICIENT (FOR MOMENT EQU)
!
          elseif(iq.eq.2.and.ivel.eq.0)then
            qq0(i,j)=1.
            qq(i,j)=1.
!
!-----GAS VOLUME FRACTION (FOR PRESSURE CORRECTION EQU)
!
          elseif(iq.eq.3.and.ivel.eq.0.and.iph.eq.1)then
            qq0(i,j)=1.-cnt*q30(i,j)
            qq(i,j)=1.-cnt*q3(i,j)
!
!-----LIQUID VOLUME FRACTION COEFFICIENT (FOR LIQ CONT EQU)
!
          elseif(iq.eq.3.and.ivel.eq.0.and.iph.eq.2)then
            qq0(i,j)=cnt
            qq(i,j)=cnt
          endif
        enddo
      enddo

      return
      end