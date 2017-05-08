!
!=====PRESSURE SOURCE TERM
!
      subroutine pressure
      include "common.inc"
      
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            if(iph.eq.1)then
              if(ivel.eq.1)then
                su(i,j,k)=su(i,j,k)
     &                    -(1.-cnt*q3(i,j,k))*dpx(i,j,k)
              elseif(ivel.eq.2)then
                su(i,j,k)=su(i,j,k)
     &                    -(1.-cnt*q3(i,j,k))*dpy(i,j,k)
              elseif(ivel.eq.3)then
                su(i,j,k)=su(i,j,k)
     &                    -(1.-cnt*q3(i,j,k))*dpz(i,j,k)
              endif
            endif
          enddo
        enddo
      enddo
 
      return
      end
!
!=====INTER-PHASE DRAG SOURCE TERM
!
      subroutine drag(velg,vell)
      include "common.inc"
      real velg(ia,ja,ka),vell(ia,ja,ka)
      
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            term1=0;term2=0;term3=0
!
!-----RELATIVE VELOCITY
!
            velu=ul3(i,j,k)-ug(i,j,k)
            velv=vl3(i,j,k)-vg(i,j,k)
            velw=wl3(i,j,k)-wg(i,j,k)
            vrel=(velu**2+velv**2+velw**2)**0.5
            vel=vell(i,j,k)-velg(i,j,k)
!
!-----REYONLDS NUMBER
!
            re=2*dng(i,j,k)*vrel*r32(i,j,k)/dvsg(i,j,k)
!
!=====THREE MOMENT SCHEME
!
            if(pk(i,j,k).gt.pkmin.and.pk(i,j,k).lt.pkmax)then
!
              if(re.lt.1000)then
                if(iq.eq.1)then
                  term1=4.5*dvsg(i,j,k)*qm1(i,j,k)
                elseif(iq.eq.2)then
                  term1=4.5*dvsg(i,j,k)*q0(i,j,k)
                elseif(iq.eq.3)then
                 term1=4.5*dvsg(i,j,k)*q1(i,j,k)
                endif
                term2=1.35*(dvsg(i,j,k)/2.)**0.313
     &            *(dng(i,j,k)*vrel)**0.687
     &            *exp(gammln(pk(i,j,k)+iq-1.313)
     &            -gammln(pk(i,j,k)))
     &            *(r32(i,j,k)/(pk(i,j,k)+2.))**(iq-1.313)
!
              elseif(re.gt.1000)then
                if(iq.eq.1)then
                  term1=0.159*dng(i,j,k)*qm1(i,j,k)
                elseif(iq.eq.2)then
                  term1=0.159*dng(i,j,k)*q0(i,j,k)
                elseif(iq.eq.3)then
                 term1=0.159*dng(i,j,k)*q1(i,j,k)
                endif
                term2=0
              endif
            endif
!
!=====DRAG SOURCE TERM
!
            term3=(term1+term2)*vel
            if(iq.eq.3)term3=(4/3.)*pi*term3
            if(iph.eq.1.and.i.gt.imn)then
              su(i,j,k)=su(i,j,k)+term3
            elseif(iph.eq.2.and.i.gt.imn)then
              su(i,j,k)=su(i,j,k)-term3
            endif
          enddo
        enddo
      enddo
      
      return
      end
!
!=====BREAK-UP SOURCE TERM
!
      subroutine breakup
      include "common.inc"
!
!-----STRIPPING CONSTANT
!
      cs=1.0 ! higher cs will reduce source
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            term1=0;term2=0;term3=0;term4s=0;term4b=0;term4=0
            if(pk(i,j,k).gt.pkmin.and.pk(i,j,k).lt.pkmax
     &        .and.q3(i,j,k).lt.q3in*0.05)then
!
!=====RELATIVE VELOCITY
!
              if(iq.eq.1)then
                velu=ul1(i,j,k)-ug(i,j,k)
                velv=vl1(i,j,k)-vg(i,j,k)
                velw=wl1(i,j,k)-wg(i,j,k)
              elseif(iq.eq.2)then
                velu=ul2(i,j,k)-ug(i,j,k)
                velv=vl2(i,j,k)-vg(i,j,k)
                velw=wl2(i,j,k)-wg(i,j,k)
              endif
C               velu=ul3(i,j,k)-ug(i,j,k)
C               velv=vl3(i,j,k)-vg(i,j,k)
C               velw=wl3(i,j,k)-wg(i,j,k)
              vrel=(velu**2+velv**2+velw**2)**0.5
              vrelmx=(ulin**2+vlin**2+wlin**2)**0.5
              if(vrel.gt.vrelmx)vrel=vrelmx
!
!=====INTEGRAL LIMITS
!
              rbag=3*st(i,j,k)/(dng(i,j,k)*vrel**2)
              rstrip=st(i,j,k)**2/(2*dng(i,j,k)*vrel**3*dvsg(i,j,k))
!
!-----FOR BAG BREAK-UP
!
              rblb=rbag
              rbub=rstrip
!
              xblb=(pk(i,j,k)+2.)*(rblb/r32(i,j,k))
              xbub=(pk(i,j,k)+2)*(rbub/r32(i,j,k))
!
!-----FOR STRIPPING BREAK-UP
!
              rslb=rstrip
              rsub=100.e-6
!
              xslb=(pk(i,j,k)+2.)*(rslb/r32(i,j,k))
              xsub=(pk(i,j,k)+2.)*(rsub/r32(i,j,k))

!
!=====TRUNCATED MOMENTS
!
!
!-----FOR BAG BREAK-UP
!
              qm0p5b=qtrunc(-1,xblb,xbub,i,j,k)
              q0p5b=qtrunc(1,xblb,xbub,i,j,k)
!
!-----FOR STRIPPING BREAK-UP
!
              q0s=qtrunc(0,xslb,xsub,i,j,k)
              q1s=qtrunc(2,xslb,xsub,i,j,k)
              q1p5s=qtrunc(3,xslb,xsub,i,j,k)
!
!=====COMMON BREAK-UP TERMS
!
              cterm1=6.2
              cterm2=cs/vrel
              cterm3=dnl(i,j,k)/dng(i,j,k)
              cterm4=dvsl(i,j,k)/(2*dnl(i,j,k)*vrel)
              cterm5=pi*(dnl(i,j,k)/(2*st(i,j,k)))**0.5
!
!=====Q1 MOMENT
!
              if(iq.eq.1)then
!
!-----BAG BREAK-UP
!
                term1=3.*qm0p5b/cterm5
!
!-----STRIPPING BREAK-UP
!
                term2=q1s/(cterm1**2*cterm2*cterm3*cterm4)
                term3=q0s/(cterm2*cterm3**0.5)

!
!=====Q2 MOMENT
!
              elseif(iq.eq.2)then
!
!-----BAG BREAK-UP
!
                term1=q0p5b/cterm5
!
!-----STRIPPING BREAK-UP
!
                term2=q1p5s/(cterm1*cterm2*cterm3**0.75*cterm4**0.5)
                term3=q1s/(cterm2*cterm3**0.5)
              endif
!
!=====BREAK-UP SOURCE TERM
!
              term4s=dnl(i,j,k)*(term2-term3) ! stripping source
              term4b=dnl(i,j,k)*(term1)       ! bag source
              term4=term4s+term4b
              if(term4.gt.0.and.i.gt.imn)then
                su(i,j,k)=su(i,j,k)+term4
                if(iq.eq.1)bq1(i,j,k)=term4
                if(iq.eq.2)bq2(i,j,k)=term4
              endif
            endif
          enddo
        enddo
      enddo

      return
      end
!        
!=====COLLISION SOURCE TERM
!
      subroutine collisions
      include "common.inc"
!
!-----NOTE: FCOLL IS DIFFERENT FROM TESTED MODEL BY BECK
!
      fcoll=0.15
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            term1=0
            if(pk(i,j,k).gt.pkmin.and.pk(i,j,k).lt.pkmax)then
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
     &          *2*pi*vrel*(q0(i,j,k)*q2(i,j,k)+q1(i,j,k)**2)
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
                xba=(pk(i,j,k)+2.)*rcrita/r32(i,j,k)
                xbb=(pk(i,j,k)+2.)*rcritb/r32(i,j,k)
                xbc=(pk(i,j,k)+2.)*rcritc/r32(i,j,k)
!
!-----PROBABILITY OF R<RCRITA
!
                proba=gammp(pk(i,j,k),xba)
!
!-----PROBABILITY OF R<RCRITB
!
                probb=gammp(pk(i,j,k),xbb)
!
!-----PROBABILITY OF R<RCRITC
!
                probc=gammp(pk(i,j,k),xbc)
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
                if(term1.gt.0.and.i.gt.imn)then
                  su(i,j,k)=su(i,j,k)+term1*dnl(i,j,k)
                endif
              endif
            endif
          enddo
        enddo
      enddo
         
      return
      end
!
!=====ADDITION TO GENERAL SOURCE TERM
!
      subroutine source(phi1,phi2)
      include "common.inc"
      real phi1(ia,ja,ka),phi2(ia,ja,ka)
      
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            su(i,j,k)=su(i,j,k)+phi1(i,j,k)*phi2(i,j,k)
          enddo
        enddo
      enddo
 
      return
      end