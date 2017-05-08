!        
!=====COLLISIONS MODEL
!
      subroutine collisions
      include "common.inc"

      order=1.e-4
!
!-----NOTE: FCOLL IS DIFFERENT FROM TESTED MODEL BY BECK
!
      fcoll=0.15
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            term1=0
            call dmom2(i,j,k)
            if(pk.gt.1.and.i.gt.imn)then
!
!-----SET AVERAGE RELATIVE VELOCITY BETWEEN DROPLETS
!
              ulsd=ul1(i,j,k)
              vlsd=vl1(i,j,k)
              wlsd=wl1(i,j,k)
              ulpart=(ul3(i,j,k)-ulsd)**2
              vlpart=(vl3(i,j,k)-vlsd)**2
              wlpart=(wl3(i,j,k)-wlsd)**2
              vrel=sqrt(ulpart+vlpart+wlpart)
!
!-----NUMBER OF COLLISIONS PER UNIT VOLUME PERUNIT TIME (*DENSITY)
!
              collno=fcoll*0.5
     &               *2*pi*vrel*(q0d*q2(i,j,k)+q1(i,j,k)**2)
!
!-----CRITICAL RADII GIVEN BY CRITICAL WEBER NUMBER
!     ASSUMING THE CRITICAL WEBER NUMBERS ARE THOSE FOR DODECANE
!     NAMELY WEA=2.5, WEB=10, WEC=30
!
              if(vrel.gt.tiny.and.collno.gt.tiny)then
                term=st(i,j,k)/(2*dng(i,j,k)*vrel*vrel)
                rcrita=2.5*term
                rcritb=10.*term
                rcritc=30.*term
!
!-----PROBABILITY OF R<RCRITA
!     NOTE: GMFT IS THE CUMULATIVE GAMMA FUNCTION
!
                proba=gmft(parkt,(pk+2)*rcrita/r32)
!
!-----PROBABILITY OF R<RCRITB
!
                probb=gmft(parkt,(pk+2)*rcritb/r32)
!
!-----PROBABILITY OF R<RCRITC
!
                probc=gmft(parkt,(pk+2)*rcritc/r32)
!
!-----PROBABILITY OF COALESCENCE IN REGION A
!
                pcoala=0.6*proba*(2-proba)
!
!-----PROBABILITY OF COALESCENCE IN REGION B
!
                pcoalb=0.5*(probc-probb)*(2-probb-probc)
!
!-----PROBABILITY OF COALESCENCE IN REGION C
!
                pcoalc=0.2*(1-probc)**2
!
!-----PROBABILITY OF SEPARATION IN REGION C
!
                psep=0.8*(1-probc)**2
!
!-----AVERAGE RADII FOR EACH REGIME
!
                rasq=rcrita**2
                rb=0.5*(rcritb+rcritc)
                rbsq=rb*rb
                rcsq=rcritc**2
!
!-----COLLISION RESULTS
!
                dq2a=-0.41*rasq
                dq2b=-0.41*rbsq
                dq2c=-0.41*rcsq
                dq2sep=0.71*rcsq
                dq1a=-0.74*rcrita
                dq1b=-0.74*rb
                dq1c=-0.74*rcritc
                dq1sep=1.684*rcritc
!
!=====Q1 MOMENT
!
                if(iq.eq.1)then
                  term1=collno*(pcoala*dq1a+pcoalb*dq1b
     &                  +pcoalc*dq1c+psep*dq1sep)
!
!=====Q2 MOMENT
!
                elseif(iq.eq.2)then
                  term1=collno*(pcoala*dq2a+pcoalb*dq2b
     &                  +pcoalc*dq2c+psep*dq2sep)
                endif
!
!=====COLLISIONS SOURCE TERM
!
                su(i,j,k)=su(i,j,k)+term1*dnl(i,j,k)*vol(i,j,k)
              endif
            endif
          enddo
        enddo
      enddo
         
      return
      end