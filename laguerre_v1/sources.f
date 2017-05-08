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
     &                    -(1.-cnt*q3(i,j,k))*dpx(i,j,k)*vol(i,j,k)
              elseif(ivel.eq.2)then
                su(i,j,k)=su(i,j,k)
     &                    -(1.-cnt*q3(i,j,k))*dpy(i,j,k)*vol(i,j,k)
              elseif(ivel.eq.3)then
                su(i,j,k)=su(i,j,k)
     &                    -(1.-cnt*q3(i,j,k))*dpz(i,j,k)*vol(i,j,k)
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
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            qa=0;qb=0;terma=0;termb=0;term1=0
            if((iqmax.eq.4
     &        .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order
     &        .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order)
     &        .or.(iqmax.eq.5
     &        .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order
     &        .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order
     &        .and.q4(i,j,k).gt.q4in*order)
     &        .or.(iqmax.eq.6
     &        .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order
     &        .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order
     &        .and.q4(i,j,k).gt.q4in*order.and.q5(i,j,k).gt.q5in*order)
     &        .or.(iqmax.eq.7
     &        .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order
     &        .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order
     &        .and.q4(i,j,k).gt.q4in*order.and.q5(i,j,k).gt.q5in*order
     &        .and.q6(i,j,k).gt.q6in*order)
     &        .or.(iqmax.eq.8
     &        .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order
     &        .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order
     &        .and.q4(i,j,k).gt.q4in*order.and.q5(i,j,k).gt.q5in*order
     &        .and.q6(i,j,k).gt.q6in*order
     &        .and.q7(i,j,k).gt.q7in*order))then
!
!-----NORMALISATION RADIUS
!
              if(rexpd.gt.tiny)rdrag=r32(i,j,k)
              if(rexpd.lt.tiny)rdrag=1.
!
!-----RELATIVE VELOCITY COEFFICIENTS
!
              urel=(vell(i,j,k)-velg(i,j,k))/rdrag
!
              termu=(uq3(i,j,k)-ug(i,j,k))/rdrag
              termv=(vq3(i,j,k)-vg(i,j,k))/rdrag
              termw=(wq3(i,j,k)-wg(i,j,k))/rdrag
              urelh=sqrt(termu**2+termv**2+termw**2)
!
!-----U /= U(r) [B. YUE THESIS; SECTION 4.5.1]
!
              if(rexpd.lt.tiny)then
                if(iq.eq.0)then
                  qa=qm2(i,j,k)
                  qb=qm1(i,j,k)
                elseif(iq.eq.1)then
                  qa=qm1(i,j,k)
                  qb=q0(i,j,k)
                elseif(iq.eq.2)then
                  qa=q0(i,j,k)
                  qb=q1(i,j,k)
                elseif(iq.eq.3)then
                  qa=q1(i,j,k)
                  qb=q2(i,j,k)
                elseif(iq.eq.4)then
                  qa=q2(i,j,k)
                  qb=q3(i,j,k)
                elseif(iq.eq.5)then
                  qa=q3(i,j,k)
                  qb=q4(i,j,k)
                elseif(iq.eq.6)then
                  qa=q4(i,j,k)
                  qb=q5(i,j,k)
                elseif(iq.eq.7)then
                  qa=q5(i,j,k)
                  qb=q6(i,j,k)
                endif
!
                terma=4.5*dvsg(i,j,k)*urel*qa
                termb=1.35*urel*(dng(i,j,k)*urelh*qb)**0.687
     &            *((dvsg(i,j,k)*qa)/2.)**0.313
!
!-----U = U(r)
!
              elseif(rexpd.gt.tiny)then
                if(iq.eq.0)then
                  qa=qa0(i,j,k)
                  qb=qb0(i,j,k)
                elseif(iq.eq.1)then
                  qa=qa1(i,j,k)
                  qb=qb1(i,j,k)
                elseif(iq.eq.2)then
                  qa=qa2(i,j,k)
                  qb=qb2(i,j,k)
                elseif(iq.eq.3)then
                  qa=qa3(i,j,k)
                  qb=qb3(i,j,k)
                elseif(iq.eq.4)then
                  qa=qa4(i,j,k)
                  qb=qb4(i,j,k)
                elseif(iq.eq.5)then
                  qa=qa5(i,j,k)
                  qb=qb5(i,j,k)
                elseif(iq.eq.6)then
                  qa=qa6(i,j,k)
                  qb=qb6(i,j,k)
                elseif(iq.eq.7)then
                  qa=qa7(i,j,k)
                  qb=qb7(i,j,k)
                endif
!
                terma=4.5*dvsg(i,j,k)*urel*qa
                termb=1.35*urel*(dng(i,j,k)*urelh)**0.687
     &            *qb*(dvsg(i,j,k)/2.)**0.313
              endif
!
!-----SOURCE TERM
!
              term1=(terma+termb)/dnl(i,j,k)
              if(iq.eq.3)term1=(4/3.)*pi*dnl(i,j,k)*term1
!
              if(iph.eq.1)su(i,j,k)=su(i,j,k)+term1*vol(i,j,k)
              if(iph.eq.2)su(i,j,k)=su(i,j,k)-term1*vol(i,j,k)
            endif
          enddo
        enddo
      enddo
      
      return
      end
!
!=====BREAK-UP SOURCE TERM [REITZ AND DIWAKAR] [1]
!
      subroutine breakup1
      include "common.inc"
!
!-----BAG CONSTANTS
!
      cb1=6.
      cb2=pi
!
!-----STRIPPING CONSTANTS
!
      cs1=0.5
      cs2=20.
!
      if(iq.eq.1)call zerovar(bq1)
      if(iq.eq.2)call zerovar(bq2)
      order=1.e-5
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            term1=0;term2=0;term3=0
            if(q1(i,j,k).gt.q1in*order
     &        .and.q2(i,j,k).gt.q2in*order
     &        .and.q3(i,j,k).gt.q3in*order
     &        .and.pk(i,j,k).gt.pkmin
     &        .and.pk(i,j,k).lt.pkmax)then
!
!-----RELATIVE VELOCITY
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
              vrel=(velu**2+velv**2+velw**2)**0.5
!
!=====BAG BREAK-UP
!
!
!-----INTEGRAL LIMITS
!
              rlb=min(cb1*st(i,j,k)
     &          /(dng(i,j,k)*vrel**2),50.e-6)
              rub=min(2.*cs1**2*st(i,j,k)**2
     &          /(dng(i,j,k)*vrel**3*dvsg(i,j,k)),50.e-6)
!
              xlb=(pk(i,j,k)+2.)*(rlb/r32(i,j,k))
              xub=(pk(i,j,k)+2.)*(rub/r32(i,j,k))
!
!-----TRUNCATED SMR
!
              q2tr=qtrunc(4,xlb,xub,i,j,k)
              q3tr=qtrunc(6,xlb,xub,i,j,k)
              if(q2tr.gt.q2in*1.e-8.and.q3tr.gt.q3in*1.e-8)then
                r32tr=q3tr/q2tr
                call wqtr(q2tr,i,j,k)
                call wqtr(q3tr,i,j,k)
!
!-----STABLE RADIUS
!
                rstab=12.4*r32tr**0.5
     &            *(dnl(i,j,k)/dng(i,j,k))**0.25
     &            *(dvsl(i,j,k)/(2.*dnl(i,j,k)*vrel))**0.5
!
!-----CHARACTERISTIC TIME SCALE
!
                tscale=cb2*(dnl(i,j,k)*r32tr**3/(2.*st(i,j,k)))**0.5
!
!-----Q1 MOMENT
!
                if(iq.eq.1)then
                  term1=(r32tr**3/rstab**2-r32tr**2)/tscale
!
!-----Q2 MOMENT
!
                elseif(iq.eq.2)then
                  term1=(r32tr**3/rstab-r32tr)/tscale
                endif
              endif
!
!=====STRIPPING BREAK-UP
!
!
!-----INTEGRAL LIMITS
!
              rlb=min(2.*cs1**2*st(i,j,k)**2
     &          /(dng(i,j,k)*vrel**3*dvsg(i,j,k)),50.e-6)
              rub=50.e-6
!
              xlb=(pk(i,j,k)+2.)*(rlb/r32(i,j,k))
              xub=(pk(i,j,k)+2.)*(rub/r32(i,j,k))
!
!-----TRUNCATED SMR
!
              q2tr=qtrunc(4,xlb,xub,i,j,k)
              q3tr=qtrunc(6,xlb,xub,i,j,k)
              if(q2tr.gt.q2in*1.e-10.and.q3tr.gt.q3in*1.e-10)then
                r32tr=q3tr/q2tr
                call wqtr(q2tr,i,j,k)
                call wqtr(q3tr,i,j,k)
!
!-----STABLE RADIUS
!
                rstab=12.4*r32tr**0.5
     &            *(dnl(i,j,k)/dng(i,j,k))**0.25
     &            *(dvsl(i,j,k)/(2.*dnl(i,j,k)*vrel))**0.5
!
!-----CHARACTERISTIC TIME SCALE
!
                tscale=cs2*r32tr*(dnl(i,j,k)/dng(i,j,k))**0.5/vrel
!
!-----Q1 MOMENT
!
                if(iq.eq.1)then
                  term2=(r32tr**3/rstab**2-r32tr**1)/tscale
!
!-----Q2 MOMENT
!
                elseif(iq.eq.2)then
                  term2=(r32tr**3/rstab-r32tr**2)/tscale
                endif
              endif
!
              term3=term1+term2
              if(i.gt.imn+1.and.term3.gt.tiny)then
                su(i,j,k)=su(i,j,k)+term3
                if(iq.eq.1)bq1(i,j,k)=term3
                if(iq.eq.2)bq2(i,j,k)=term3
              endif
            endif
          enddo
        enddo
      enddo

      return
      end
!
!=====BREAK-UP SOURCE TERM  [PILCH AND ERDMAN MODEL] [2]
!
      subroutine breakup2
      include "common.inc"
!
      if(iq.eq.1)call zerovar(bq1)
      if(iq.eq.2)call zerovar(bq2)
      order=1.e-5
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            term1=0;term2=0;term3=0;term4=0;term5=0;term6=0
            if(q1(i,j,k).gt.q1in*order
     &        .and.q2(i,j,k).gt.q2in*order
     &        .and.q3(i,j,k).gt.q3in*order
     &        .and.pk(i,j,k).gt.pkmin
     &        .and.pk(i,j,k).lt.pkmax)then
!
!-----RELATIVE VELOCITY
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
              vrel=(velu**2+velv**2+velw**2)**0.5
!
!-----VIBRATIONAL BREAK-UP
!
              term1=dqdt(1,vrel,i,j,k)
!
!-----BAG BREAK-UP
!
              term2=dqdt(2,vrel,i,j,k)
!
!-----BAG AND STAMEN BREAK-UP
!
              term3=dqdt(3,vrel,i,j,k)
!
!-----SHEET STRIPPING
!
              term4=dqdt(4,vrel,i,j,k)
!
!-----WAVE CREST STRIPPING
!
              term5=dqdt(5,vrel,i,j,k)
!
              term6=term1+term2+term3+term4+term5
              if(i.gt.imn+1.and.term6.gt.tiny)then
                su(i,j,k)=su(i,j,k)+term6
                if(iq.eq.1)bq1(i,j,k)=term6
                if(iq.eq.2)bq2(i,j,k)=term6
              endif
            endif
          enddo
        enddo
      enddo

      return
      end
!
!=====DQDT
!
      function dqdt(ireg,vrel,i,j,k)
      include "common.inc"
      dqdt=0
!
!-----LIMITS
!
      if(ireg.eq.1)coeflb=12.;coefub=18.
      if(ireg.eq.2)coeflb=18.;coefub=45.
      if(ireg.eq.3)coeflb=45.;coefub=350.
      if(ireg.eq.4)coeflb=350.;coefub=2670.
      if(ireg.eq.5)coeflb=2670.;coefub=5000.
!
      rlb=min(coeflb*st(i,j,k)/(dng(i,j,k)*vrel**2),50.e-6)
      rub=min(coefub*st(i,j,k)/(dng(i,j,k)*vrel**2),50.e-6)
!
      xlb=(pk(i,j,k)+2.)*(rlb/r32(i,j,k))
      xub=(pk(i,j,k)+2.)*(rub/r32(i,j,k))
!
!-----TRUNCATED MOMENTS
!
      q2tr=qtrunc(4,xlb,xub,i,j,k)
      q3tr=qtrunc(6,xlb,xub,i,j,k)
      if(q2tr.gt.q2in*1.e-10.and.q3tr.gt.q3in*1.e-10)then
        r32tr=q3tr/q2tr
        call wqtr(q2tr,i,j,k)
        call wqtr(q3tr,i,j,k)
!
!-----CHARACTERISTIC TIME SCALE
!
        if(we(vrel,rlb,i,j,k).gt.wecrit(rlb,i,j,k)
     &    .and.we(vrel,rub,i,j,k).gt.wecrit(rub,i,j,k)
     &    .and.we(vrel,r32tr,i,j,k).gt.12.)then
          if(ireg.eq.1)tcoef=6.*(we(vrel,r32tr,i,j,k)-12.)**(-0.25)
          if(ireg.eq.2)tcoef=2.45*(we(vrel,r32tr,i,j,k)-12.)**0.25
          if(ireg.eq.3)tcoef=14.1*(we(vrel,r32tr,i,j,k)-12.)**(-0.25)
          if(ireg.eq.4)tcoef=.766*(we(vrel,r32tr,i,j,k)-12.)**0.25
          if(ireg.eq.5)tcoef=5.5
          tscale=2.*tcoef*r32tr*(dnl(i,j,k)/dng(i,j,k))**0.5/vrel
!
!-----STABLE RADIUS
!
          rcoef=vrel*(dng(i,j,k)/dnl(i,j,k))**0.5
     &      *(0.375*tcoef+0.2274*tcoef**2)
          rstab=wecrit(r32tr,i,j,k)*st(i,j,k)
     &      /(dng(i,j,k)*vrel**2*(1.-rcoef/vrel)**2)
!
          if(iq.eq.1)then
            dqdt=(r32tr**3/rstab**2-r32tr)/tscale
          elseif(iq.eq.2)then
            dqdt=(r32tr**3/rstab-r32tr**2)/tscale
          endif
        endif
      endif

      return
      end
!
!=====BREAK-UP SOURCE TERM  [HSIANG AND FAETH MODEL] [3]
!
      subroutine breakup
      include "common.inc"
!
      if(iq.eq.1)call zerovar(bq1)
      if(iq.eq.2)call zerovar(bq2)
      order=1.e-6
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            term1=0
            if(q1(i,j,k).gt.q1in*order
     &        .and.q2(i,j,k).gt.q2in*order
     &        .and.q3(i,j,k).gt.q3in*order
     &        .and.pk(i,j,k).gt.pkmin
     &        .and.pk(i,j,k).lt.pkmax)then
!
!-----RELATIVE VELOCITY
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
              vrel=(velu**2+velv**2+velw**2)**0.5
!
!-----INTEGRAL LIMITS
!
              rlb=min(6.*st(i,j,k)/(dng(i,j,k)*vrel**2),50.e-6)
              rub=50.e-6
!
              xlb=(pk(i,j,k)+2.)*(rlb/r32(i,j,k))
              xub=(pk(i,j,k)+2)*(rub/r32(i,j,k))
!
!-----TRUNCATED MOMENTS
!
              q2tr=qtrunc(4,xlb,xub,i,j,k)
              q3tr=qtrunc(6,xlb,xub,i,j,k)
              if(q2tr.gt.q2in*1.e-10.and.q3tr.gt.q3in*1.e-10)then
                r32tr=q3tr/q2tr
!
!-----STABLE RADIUS
!
                rstab=12.4*r32tr**0.5
     &            *(dnl(i,j,k)/dng(i,j,k))**0.25
     &            *(dvsl(i,j,k)/(2.*dnl(i,j,k)*vrel))**0.5
!
!-----CHARACTERISTIC TIME SCALE
!
                tscale=2*5.*r32tr*(dnl(i,j,k)/dng(i,j,k))**0.5
     &            /(vrel*(1.-oh(r32tr,i,j,k)/7.))
!
!-----Q1 MOMENT SOURCE TERM
!
                if(iq.eq.1)then
                  term1=(r32tr**3/rstab**2-r32tr**1)/tscale
!
!-----Q2 MOMENT SOURCE TERM
!
                elseif(iq.eq.2)then
                  term1=(r32tr**3/rstab**1-r32tr**2)/tscale
                endif
              endif
!
!-----NO SOURCE AT INJECTOR
!
              if(i.eq.imn.and.j.eq.jinj.and.k.eq.kinj)term1=0.
              if(i.eq.imn.and.j.eq.jinj+1.and.k.eq.kinj)term1=0.
              if(i.eq.imn.and.j.eq.jinj.and.k.eq.kinj+1)term1=0.
              if(i.eq.imn.and.j.eq.jinj+1.and.k.eq.kinj+1)term1=0.
!
              if(term1.lt.tiny)term1=0.
              su(i,j,k)=su(i,j,k)+term1
              if(iq.eq.1)bq1(i,j,k)=term1
              if(iq.eq.2)bq2(i,j,k)=term1
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

      order=1.e-6
!
!-----NOTE: FCOLL IS DIFFERENT FROM TESTED MODEL BY BECK
!
      fcoll=0.15
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            term1=0
            if(q1(i,j,k).gt.q1in*order
     &        .and.q2(i,j,k).gt.q2in*order
     &        .and.q3(i,j,k).gt.q3in*order
     &        .and.pk(i,j,k).gt.pkmin
     &        .and.pk(i,j,k).lt.pkmax)then
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
!-----NO SOURCE AT INJECTOR
!
                if(i.eq.imn.and.j.eq.jinj.and.k.eq.kinj)term1=0.
                if(i.eq.imn.and.j.eq.jinj+1.and.k.eq.kinj)term1=0.
                if(i.eq.imn.and.j.eq.jinj.and.k.eq.kinj+1)term1=0.
                if(i.eq.imn.and.j.eq.jinj+1.and.k.eq.kinj+1)term1=0.
!
                if(term1.lt.tiny)term1=0.
                su(i,j,k)=su(i,j,k)+term1
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
!
!=====VISCOUS STRESS SOURCE TERM
!
      subroutine vstress
      include "common.inc"
      real uvel(ia,ja,ka),vvel(ia,ja,ka),wvel(ia,ja,ka)
      real temp(ia,ja,ka)

      ivst=0
      if(ivst.eq.0)goto 100

      call zerovar(temp)
!
!-----PHASE
!
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
!
!-----GAS
!
            if(iph.eq.1)then
              uvel(i,j,k)=ug(i,j,k)
              vvel(i,j,k)=vg(i,j,k)
              wvel(i,j,k)=wg(i,j,k)
!
!-----LIQUID
!
            elseif(iph.eq.2)then
              if(iq.eq.1)then
                uvel(i,j,k)=ul1(i,j,k)
                vvel(i,j,k)=vl1(i,j,k)
                wvel(i,j,k)=wl1(i,j,k)
              elseif(iq.eq.2)then
                uvel(i,j,k)=ul2(i,j,k)
                vvel(i,j,k)=vl2(i,j,k)
                wvel(i,j,k)=wl2(i,j,k)
              elseif(iq.eq.3)then
                uvel(i,j,k)=ul3(i,j,k)
                vvel(i,j,k)=vl3(i,j,k)
                wvel(i,j,k)=wl3(i,j,k)
              endif
            endif
          enddo
        enddo
      enddo
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            term1=0;term2=0;term3=0;term4=0
!
!-----X MOMENTUM
!
            if(ivel.eq.1)then
!
!-----d/dx(mu.du/dx)
!
              term1a=f(1,i,j,k,dd)
              term1b=f(2,i,j,k,dd)
              
              term1c=uvel(i,j,k)-uvel(i-1,j,k)
              term1d=uvel(i+1,j,k)-uvel(i,j,k)
              
              term1e=xc(i)-xc(i-1)
              term1f=xc(i+1)-xc(i)
              
              term1g=dx(i)
              
              term1=(term1b*(term1d/term1f)-term1a*(term1c/term1e))
     &          /term1g
!
!-----d/dy(mu.dv/dx)
!
              term1a=f(3,i,j,k,dd)
              term1b=f(4,i,j,k,dd)
              
              temp(i-1,j,k)=f(3,i-1,j,k,vvel)
              temp(i,j,k)=f(3,i,j,k,vvel)
              temp(i+1,j,k)=f(3,i+1,j,k,vvel)
              vel1a=f(1,i,j,k,temp)
              vel1b=f(2,i,j,k,temp)
              term1c=vel1b-vel1a
              temp(i-1,j,k)=f(4,i-1,j,k,vvel)
              temp(i,j,k)=f(4,i,j,k,vvel)
              temp(i+1,j,k)=f(4,i+1,j,k,vvel)
              vel2a=f(1,i,j,k,temp)
              vel2b=f(2,i,j,k,temp)
              term1d=vel2b-vel2a
              
              term1e=dx(i)
              term1f=dx(i)
              
              term1g=dy(j)
              
              term2=(term1b*(term1d/term1f)-term1a*(term1c/term1e))
     &          /term1g
!
!-----d/dz(mu.dw/dx)
!
              term1a=f(5,i,j,k,dd)
              term1b=f(6,i,j,k,dd)
              
              temp(i-1,j,k)=f(5,i-1,j,k,wvel)
              temp(i,j,k)=f(5,i,j,k,wvel)
              temp(i+1,j,k)=f(5,i+1,j,k,wvel)
              vel1a=f(1,i,j,k,temp)
              vel1b=f(2,i,j,k,temp)
              term1c=vel1b-vel1a
              temp(i-1,j,k)=f(6,i-1,j,k,wvel)
              temp(i,j,k)=f(6,i,j,k,wvel)
              temp(i+1,j,k)=f(6,i+1,j,k,wvel)
              vel2a=f(1,i,j,k,temp)
              vel2b=f(2,i,j,k,temp)
              term1d=vel2b-vel2a
              
              term1e=dx(i)
              term1f=dx(i)
              
              term1g=dz(k)
              
              term3=(term1b*(term1d/term1f)-term1a*(term1c/term1e))
     &          /term1g
!
!-----Y MOMENTUM
!
            elseif(ivel.eq.2)then
!
!-----d/dx(mu.du/dy)
!
              term1a=f(1,i,j,k,dd)
              term1b=f(2,i,j,k,dd)
              
              temp(i,j-1,k)=f(1,i,j-1,k,uvel)
              temp(i,j,k)=f(1,i,j,k,uvel)
              temp(i,j+1,k)=f(1,i,j+1,k,uvel)
              vel1a=f(3,i,j,k,temp)
              vel1b=f(4,i,j,k,temp)
              term1c=vel1b-vel1a
              temp(i,j-1,k)=f(2,i,j-1,k,uvel)
              temp(i,j,k)=f(2,i,j,k,uvel)
              temp(i,j+1,k)=f(2,i,j+1,k,uvel)
              vel2a=f(3,i,j,k,temp)
              vel2b=f(4,i,j,k,temp)
              term1d=vel2b-vel2a
              
              term1e=dy(j)
              term1f=dy(j)
              
              term1g=dx(i)
              
              term1=(term1b*(term1d/term1f)-term1a*(term1c/term1e))
     &          /term1g
!
!-----d/dy(mu.dv/dy)
!
              term1a=f(3,i,j,k,dd)
              term1b=f(4,i,j,k,dd)
              
              term1c=vvel(i,j,k)-vvel(i,j-1,k)
              term1d=vvel(i,j+1,k)-vvel(i,j,k)
              
              term1e=yc(j)-yc(j-1)
              term1f=yc(j+1)-yc(j)
              
              term1g=dy(j)
              
              term2=(term1b*(term1d/term1f)-term1a*(term1c/term1e))
     &          /term1g
!
!-----d/dz(mu.dw/dy)
!
              term1a=f(5,i,j,k,dd)
              term1b=f(6,i,j,k,dd)
              
              temp(i,j-1,k)=f(5,i,j-1,k,wvel)
              temp(i,j,k)=f(5,i,j,k,wvel)
              temp(i,j+1,k)=f(5,i,j+1,k,wvel)
              vel1a=f(3,i,j,k,temp)
              vel1b=f(4,i,j,k,temp)
              term1c=vel1b-vel1a
              temp(i,j-1,k)=f(6,i,j-1,k,wvel)
              temp(i,j,k)=f(6,i,j,k,wvel)
              temp(i,j+1,k)=f(6,i,j+1,k,wvel)
              vel2a=f(3,i,j,k,temp)
              vel2b=f(4,i,j,k,temp)
              term1d=vel2b-vel2a
              
              term1e=dy(j)
              term1f=dy(j)
              
              term1g=dz(k)
              
              term3=(term1b*(term1d/term1f)-term1a*(term1c/term1e))
     &          /term1g
!
!-----Z MOMENTUM
!
            elseif(ivel.eq.3)then
!
!-----d/dx(mu.du/dz)
!
              term1a=f(1,i,j,k,dd)
              term1b=f(2,i,j,k,dd)
              
              temp(i,j,k-1)=f(1,i,j,k-1,uvel)
              temp(i,j,k)=f(1,i,j,k,uvel)
              temp(i,j,k+1)=f(1,i,j,k+1,uvel)
              vel1a=f(5,i,j,k,temp)
              vel1b=f(6,i,j,k,temp)
              term1c=vel1b-vel1a
              temp(i,j,k-1)=f(2,i,j,k-1,uvel)
              temp(i,j,k)=f(2,i,j,k,uvel)
              temp(i,j,k+1)=f(2,i,j,k+1,uvel)
              vel2a=f(5,i,j,k,temp)
              vel2b=f(6,i,j,k,temp)
              term1d=vel2b-vel2a
              
              term1e=dz(k)
              term1f=dz(k)
              
              term1g=dx(i)
              
              term1=(term1b*(term1d/term1f)-term1a*(term1c/term1e))
     &          /term1g
!
!-----d/dy(mu.dv/dz)
!
              term1a=f(3,i,j,k,dd)
              term1b=f(4,i,j,k,dd)
              
              temp(i,j,k-1)=f(3,i,j,k-1,vvel)
              temp(i,j,k)=f(3,i,j,k,vvel)
              temp(i,j,k+1)=f(3,i,j,k+1,vvel)
              vel1a=f(5,i,j,k,temp)
              vel1b=f(6,i,j,k,temp)
              term1c=vel1b-vel1a
              temp(i,j,k-1)=f(4,i,j,k-1,vvel)
              temp(i,j,k)=f(4,i,j,k,vvel)
              temp(i,j,k+1)=f(4,i,j,k+1,vvel)
              vel2a=f(5,i,j,k,temp)
              vel2b=f(6,i,j,k,temp)
              term1d=vel2b-vel2a
              
              term1e=dz(k)
              term1f=dz(k)
              
              term1g=dy(j)
              
              term2=(term1b*(term1d/term1f)-term1a*(term1c/term1e))
     &          /term1g
!
!-----d/dz(mu.dw/dz)
!
              term1a=f(5,i,j,k,dd)
              term1b=f(6,i,j,k,dd)
              
              term1c=wvel(i,j,k)-wvel(i,j,k-1)
              term1d=wvel(i,j,k+1)-wvel(i,j,k)
              
              term1e=zc(k)-zc(k-1)
              term1f=zc(k+1)-zc(k)
              
              term1g=dz(k)
              
              term3=(term1b*(term1d/term1f)-term1a*(term1c/term1e))
     &          /term1g
            endif
!
            term4=(term1+term2+term3)*vol(i,j,k)
!
!-----NO SOURCE AT INJECTOR
!
            if(i.le.imn+1)term4=0.
            if(i.eq.imn.and.j.eq.jinj.and.k.eq.kinj)term4=0.
            if(i.eq.imn.and.j.eq.jinj+1.and.k.eq.kinj)term4=0.
            if(i.eq.imn.and.j.eq.jinj.and.k.eq.kinj+1)term4=0.
            if(i.eq.imn.and.j.eq.jinj+1.and.k.eq.kinj+1)term4=0.
!
            su(i,j,k)=su(i,j,k)+term4
          enddo
        enddo
      enddo

100   continue
      return
      end