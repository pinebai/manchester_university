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
                su(i,j,k)=su(i,j,k)&
                  -(1.-cnt*q3(i,j,k))*dpx(i,j,k)*vol(i,j,k)
              elseif(ivel.eq.2)then
                su(i,j,k)=su(i,j,k)&
                  -(1.-cnt*q3(i,j,k))*dpy(i,j,k)*vol(i,j,k)
              elseif(ivel.eq.3)then
                su(i,j,k)=su(i,j,k)&
                  -(1.-cnt*q3(i,j,k))*dpz(i,j,k)*vol(i,j,k)
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
            if((iqmax.eq.4&
              .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order&
              .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order)&
              .or.(iqmax.eq.5&
              .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order&
              .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order&
              .and.q4(i,j,k).gt.q4in*order)&
              .or.(iqmax.eq.6&
              .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order&
              .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order&
              .and.q4(i,j,k).gt.q4in*order.and.q5(i,j,k).gt.q5in*order)&
              .or.(iqmax.eq.7&
              .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order&
              .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order&
              .and.q4(i,j,k).gt.q4in*order.and.q5(i,j,k).gt.q5in*order&
              .and.q6(i,j,k).gt.q6in*order)&
              .or.(iqmax.eq.8&
              .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order&
              .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order&
              .and.q4(i,j,k).gt.q4in*order.and.q5(i,j,k).gt.q5in*order&
              .and.q6(i,j,k).gt.q6in*order.and.q7(i,j,k).gt.q7in*order))then
!
!-----NORMALISATION RADIUS
!
              if(rexpd.gt.tiny)then
                call moment(qa,3.+rexpd,0.,0.,i,j,k)
                call moment(qb,3.,0.,0.,i,j,k)
                rdrag=qa/(qb+tiny)
                r43=q4(i,j,k)/q3(i,j,k)
                if(rdrag.gt.1.or.rdrag.lt.r43)goto 100
              endif
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
                if(qa.gt.tiny.and.qb.gt.tiny)then
                  terma=4.5*dvsg(i,j,k)*urel*qa
                  termb=1.35*urel*(dng(i,j,k)*urelh*qb)**0.687&
                    *((dvsg(i,j,k)*qa)/2.)**0.313
                endif
!
!-----U = U(r)
!
              elseif(rexpd.gt.tiny.and.idrag.eq.0)then
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
                if(qa.gt.tiny.and.qb.gt.tiny)then
                  terma=4.5*dvsg(i,j,k)*urel*qa
                  termb=1.35*urel*(dng(i,j,k)*urelh)**0.687&
                    *qb*(dvsg(i,j,k)/2.)**0.313
                endif
!
!-----U = U(r) ALTERNATIVE
!
              elseif(rexpd.gt.tiny.and.idrag.eq.1)then
                power=iq-2.+rexpd
                ilb=int(abs(power))
                iub=nint(abs(power))
                if(iub.eq.ilb)iub=iub+1
                ilb=sign(real(ilb),power)
                iub=sign(real(iub),power)
                d=abs(power-ilb)
                if(ilb.eq.-2)qa=qm2(i,j,k)
                if(ilb.eq.-1)qa=qm1(i,j,k)
                if(ilb.eq.0)qa=q0(i,j,k)
                if(ilb.eq.1)qa=q1(i,j,k)
                if(ilb.eq.2)qa=q2(i,j,k)
                if(ilb.eq.3)qa=q3(i,j,k)
                if(ilb.eq.4)qa=q4(i,j,k)
                if(ilb.eq.5)qa=q5(i,j,k)
                if(ilb.eq.6)qa=q6(i,j,k)
                if(ilb.eq.7)qa=q7(i,j,k)
                if(iub.eq.-2)qb=qm2(i,j,k)
                if(iub.eq.-1)qb=qm1(i,j,k)
                if(iub.eq.0)qb=q0(i,j,k)
                if(iub.eq.1)qb=q1(i,j,k)
                if(iub.eq.2)qb=q2(i,j,k)
                if(iub.eq.3)qb=q3(i,j,k)
                if(iub.eq.4)qb=q4(i,j,k)
                if(iub.eq.5)qb=q5(i,j,k)
                if(iub.eq.6)qb=q6(i,j,k)
                if(iub.eq.7)qb=q7(i,j,k)
!
                power=iq-1.313+1.687*rexpd
                ilb=int(abs(power))
                iub=nint(abs(power))
                if(iub.eq.ilb)iub=iub+1
                ilb=sign(real(ilb),power)
                iub=sign(real(iub),power)
                d=abs(power-ilb)
                if(ilb.eq.-2)qc=qm2(i,j,k)
                if(ilb.eq.-1)qc=qm1(i,j,k)
                if(ilb.eq.0)qc=q0(i,j,k)
                if(ilb.eq.1)qc=q1(i,j,k)
                if(ilb.eq.2)qc=q2(i,j,k)
                if(ilb.eq.3)qc=q3(i,j,k)
                if(ilb.eq.4)qc=q4(i,j,k)
                if(ilb.eq.5)qc=q5(i,j,k)
                if(ilb.eq.6)qc=q6(i,j,k)
                if(ilb.eq.7)qc=q7(i,j,k)
                if(iub.eq.-2)qd=qm2(i,j,k)
                if(iub.eq.-1)qd=qm1(i,j,k)
                if(iub.eq.0)qd=q0(i,j,k)
                if(iub.eq.1)qd=q1(i,j,k)
                if(iub.eq.2)qd=q2(i,j,k)
                if(iub.eq.3)qd=q3(i,j,k)
                if(iub.eq.4)qd=q4(i,j,k)
                if(iub.eq.5)qd=q5(i,j,k)
                if(iub.eq.6)qd=q6(i,j,k)
                if(iub.eq.7)qd=q7(i,j,k)
!
                if(qa.gt.tiny.and.qb.gt.tiny)then
                  terma=4.5*dvsg(i,j,k)*urel*(qa)**(1-d)*(qb)**(d)
                  termb=1.35*urel*(dng(i,j,k)*urelh)**0.687&
                    *(qc)**(1-d)*(qd)**(d)*(dvsg(i,j,k)/2.)**0.313
                endif
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
100         continue
          enddo
        enddo
      enddo
      
      return
      end
!
!=====BREAK-UP  [HSIANG AND FAETH MODEL]
!
      subroutine breakup
      include "common.inc"
!
      if(iq.eq.0)call zerovar(bq0)
      if(iq.eq.1)call zerovar(bq1)
      if(iq.eq.2)call zerovar(bq2)
      if(iq.eq.3)call zerovar(bq3)
      if(iq.eq.4)call zerovar(bq4)
      if(iq.eq.5)call zerovar(bq5)
      if(iq.eq.6)call zerovar(bq6)
      if(iq.eq.7)call zerovar(bq7)
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            term=0
            if((iqmax.eq.4&
              .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order&
              .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order)&
              .or.(iqmax.eq.5&
              .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order&
              .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order&
              .and.q4(i,j,k).gt.q4in*order)&
              .or.(iqmax.eq.6&
              .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order&
              .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order&
              .and.q4(i,j,k).gt.q4in*order.and.q5(i,j,k).gt.q5in*order)&
              .or.(iqmax.eq.7&
              .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order&
              .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order&
              .and.q4(i,j,k).gt.q4in*order.and.q5(i,j,k).gt.q5in*order&
              .and.q6(i,j,k).gt.q6in*order)&
              .or.(iqmax.eq.8&
              .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order&
              .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order&
              .and.q4(i,j,k).gt.q4in*order.and.q5(i,j,k).gt.q5in*order&
              .and.q6(i,j,k).gt.q6in*order.and.q7(i,j,k).gt.q7in*order))then
!
!-----RELATIVE VELOCITY
!
              if(iq.eq.0)then
                uvel=uq0(i,j,k)-ug(i,j,k)
                vvel=vq0(i,j,k)-vg(i,j,k)
                wvel=wq0(i,j,k)-wg(i,j,k)
              elseif(iq.eq.1)then
                uvel=uq1(i,j,k)-ug(i,j,k)
                vvel=vq1(i,j,k)-vg(i,j,k)
                wvel=wq1(i,j,k)-wg(i,j,k)
              elseif(iq.eq.2)then
                uvel=uq2(i,j,k)-ug(i,j,k)
                vvel=vq2(i,j,k)-vg(i,j,k)
                wvel=wq2(i,j,k)-wg(i,j,k)
              elseif(iq.eq.3)then
                uvel=uq3(i,j,k)-ug(i,j,k)
                vvel=vq3(i,j,k)-vg(i,j,k)
                wvel=wq3(i,j,k)-wg(i,j,k)
              elseif(iq.eq.4)then
                uvel=uq4(i,j,k)-ug(i,j,k)
                vvel=vq4(i,j,k)-vg(i,j,k)
                wvel=wq4(i,j,k)-wg(i,j,k)
              elseif(iq.eq.5)then
                uvel=uq5(i,j,k)-ug(i,j,k)
                vvel=vq5(i,j,k)-vg(i,j,k)
                wvel=wq5(i,j,k)-wg(i,j,k)
              elseif(iq.eq.6)then
                uvel=uq6(i,j,k)-ug(i,j,k)
                vvel=vq6(i,j,k)-vg(i,j,k)
                wvel=wq6(i,j,k)-wg(i,j,k)
              elseif(iq.eq.7)then
                uvel=uq7(i,j,k)-ug(i,j,k)
                vvel=vq7(i,j,k)-vg(i,j,k)
                wvel=wq7(i,j,k)-wg(i,j,k)
              endif
              vrel=sqrt((uvel**2+vvel**2+wvel**2))
!
!-----INTEGRAL LIMITS
!
              rlb=min(6.*st(i,j,k)/(dng(i,j,k)*vrel**2),50.e-6)
!
!-----TRUNCATED MOMENTS
!
              call moment(q2tr,2.,rlb,0.,i,j,k)
              call moment(q3tr,3.,rlb,0.,i,j,k)
              if(q2tr.gt.q2in*1.e-10.and.q3tr.gt.q3in*1.e-10)then
                r32tr=q3tr/q2tr
!
!-----STABLE RADIUS
!
                rstab=12.4*r32tr**0.5&
                  *(dnl(i,j,k)/dng(i,j,k))**0.25&
                  *(dvsl(i,j,k)/(2.*dnl(i,j,k)*vrel))**0.5
!
!-----CHARACTERISTIC TIME SCALE
!
                tscale=2*5.*r32tr*(dnl(i,j,k)/dng(i,j,k))**0.5&
                  /(vrel*(1.-oh(r32tr,i,j,k)/7.))
!
!-----Q0 MOMENT SOURCE TERM
!
                if(iq.eq.0)then
                  term=(r32tr**3/rstab**3-r32tr**0)/tscale
!
!-----Q1 MOMENT SOURCE TERM
!
                elseif(iq.eq.1)then
                  term=(r32tr**3/rstab**2-r32tr**1)/tscale
!
!-----Q2 MOMENT SOURCE TERM
!
                elseif(iq.eq.2)then
                  term=(r32tr**3/rstab**1-r32tr**2)/tscale
!
!-----Q3 MOMENT SOURCE TERM
!
                elseif(iq.eq.3)then
                  term=(r32tr**3/rstab**0-r32tr**3)/tscale ! which equals zero.
!
!-----Q4 MOMENT SOURCE TERM
!
                elseif(iq.eq.4)then
                  term=(r32tr**3/rstab**(-1)-r32tr**4)/tscale
!
!-----Q5 MOMENT SOURCE TERM
!
                elseif(iq.eq.5)then
                  term=(r32tr**3/rstab**(-2)-r32tr**5)/tscale
!
!-----Q6 MOMENT SOURCE TERM
!
                elseif(iq.eq.6)then
                  term=(r32tr**3/rstab**(-3)-r32tr**6)/tscale
!
!-----Q7 MOMENT SOURCE TERM
!
                elseif(iq.eq.7)then
                  term=(r32tr**3/rstab**(-4)-r32tr**7)/tscale
                endif
              endif
!
              if(term.lt.tiny)term=0.
              su(i,j,k)=su(i,j,k)+term*vol(i,j,k)
              if(iq.eq.0)bq0(i,j,k)=term
              if(iq.eq.1)bq1(i,j,k)=term
              if(iq.eq.2)bq2(i,j,k)=term
              if(iq.eq.3)bq3(i,j,k)=term
              if(iq.eq.4)bq4(i,j,k)=term
              if(iq.eq.5)bq5(i,j,k)=term
              if(iq.eq.6)bq6(i,j,k)=term
              if(iq.eq.7)bq7(i,j,k)=term
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
      if(iq.eq.0)call zerovar(cq0)
      if(iq.eq.1)call zerovar(cq1)
      if(iq.eq.2)call zerovar(cq2)
      if(iq.eq.3)call zerovar(cq3)
      if(iq.eq.4)call zerovar(cq4)
      if(iq.eq.5)call zerovar(cq5)
      if(iq.eq.6)call zerovar(cq6)
      if(iq.eq.7)call zerovar(cq7)
!
!-----NOTE: FCOLL IS DIFFERENT FROM TESTED MODEL BY BECK
!
      fcoll=0.15
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            term=0
            if((iqmax.eq.4&
              .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order&
              .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order)&
              .or.(iqmax.eq.5&
              .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order&
              .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order&
              .and.q4(i,j,k).gt.q4in*order)&
              .or.(iqmax.eq.6&
              .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order&
              .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order&
              .and.q4(i,j,k).gt.q4in*order.and.q5(i,j,k).gt.q5in*order)&
              .or.(iqmax.eq.7&
              .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order&
              .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order&
              .and.q4(i,j,k).gt.q4in*order.and.q5(i,j,k).gt.q5in*order&
              .and.q6(i,j,k).gt.q6in*order)&
              .or.(iqmax.eq.8&
              .and.q0(i,j,k).gt.q0in*order.and.q1(i,j,k).gt.q1in*order&
              .and.q2(i,j,k).gt.q2in*order.and.q3(i,j,k).gt.q3in*order&
              .and.q4(i,j,k).gt.q4in*order.and.q5(i,j,k).gt.q5in*order&
              .and.q6(i,j,k).gt.q6in*order.and.q7(i,j,k).gt.q7in*order))then
!
!-----SET AVERAGE RELATIVE VELOCITY BETWEEN DROPLETS
!
              uvel=uq3(i,j,k)-uq0(i,j,k)
              vvel=vq3(i,j,k)-vq0(i,j,k)
              wvel=wq3(i,j,k)-wq0(i,j,k)
              vrel=sqrt((uvel**2+vvel**2+wvel**2))
!
!-----NUMBER OF COLLISIONS PER UNIT VOLUME PERUNIT TIME (*DENSITY)
!
              collno=fcoll*pi*vrel*(q0(i,j,k)*q2(i,j,k)+q1(i,j,k)**2)
!
!-----CRITICAL RADII GIVEN BY CRITICAL WEBER NUMBER
!     ASSUMING THE CRITICAL WEBER NUMBERS ARE THOSE FOR DODECANE
!     NAMELY WEA=2.5, WEB=10, WEC=30
!
              call moment(qq0,0.,0.,0.,i,j,k)
              if(vrel.gt.tiny.and.collno.gt.tiny.and.qq0.gt.tiny)then
                coeff=st(i,j,k)/(2*dng(i,j,k)*vrel**2)
                rlba=min(2.5*coeff,50.e-6)
                rlbb=min(10.*coeff,50.e-6)
                rlbc=min(30.*coeff,50.e-6)
!
!-----PROBABILITY OF R<RCRITA
!
                call moment(proba,0.,rlba,0.,i,j,k)
                proba=proba/(qq0+tiny)
!
!-----PROBABILITY OF R<RCRITB
!
                call moment(probb,0.,rlbb,0.,i,j,k)
                probb=probb/(qq0+tiny)
!
!-----PROBABILITY OF R<RCRITC
!
                call moment(probc,0.,rlbc,0.,i,j,k)
                probc=probc/(qq0+tiny)
                print*,'proba,probb,probc',proba,probb,probc
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
                ra=rlba
                rb=(rlbb+rlbc)/2
                rc=rlbc
                rd=rlbc
!
!-----Q0 MOMENT
!
                if(iq.eq.0)then
                  dqa=(2**(iq/3.)-2)*ra**(iq)
                  dqb=(2**(iq/3.)-2)*rb**(iq)
                  dqc=(2**(iq/3.)-2)*rc**(iq)
                  dqd=(5*0.4**(iq/3.)-2)*rd**(iq)
!
!-----Q1 MOMENT
!
                elseif(iq.eq.1)then
                  dqa=(2**(iq/3.)-2)*ra**(iq)
                  dqb=(2**(iq/3.)-2)*rb**(iq)
                  dqc=(2**(iq/3.)-2)*rc**(iq)
                  dqd=(5*0.4**(iq/3.)-2)*rd**(iq)
!
!-----Q2 MOMENT
!
                elseif(iq.eq.2)then
                  dqa=(2**(iq/3.)-2)*ra**(iq)
                  dqb=(2**(iq/3.)-2)*rb**(iq)
                  dqc=(2**(iq/3.)-2)*rc**(iq)
                  dqd=(5*0.4**(iq/3.)-2)*rd**(iq)
!
!-----Q3 MOMENT
!
                elseif(iq.eq.3)then
                  dqa=(2**(iq/3.)-2)*ra**(iq)
                  dqb=(2**(iq/3.)-2)*rb**(iq)
                  dqc=(2**(iq/3.)-2)*rc**(iq)
                  dqd=(5*0.4**(iq/3.)-2)*rd**(iq) ! which all equal zero.
!
!-----Q4 MOMENT
!
                elseif(iq.eq.4)then
                  dqa=(2**(iq/3.)-2)*ra**(iq)
                  dqb=(2**(iq/3.)-2)*rb**(iq)
                  dqc=(2**(iq/3.)-2)*rc**(iq)
                  dqd=(5*0.4**(iq/3.)-2)*rd**(iq)
!
!-----Q5 MOMENT
!
                elseif(iq.eq.5)then
                  dqa=(2**(iq/3.)-2)*ra**(iq)
                  dqb=(2**(iq/3.)-2)*rb**(iq)
                  dqc=(2**(iq/3.)-2)*rc**(iq)
                  dqd=(5*0.4**(iq/3.)-2)*rd**(iq)
!
!-----Q6 MOMENT
!
                elseif(iq.eq.6)then
                  dqa=(2**(iq/3.)-2)*ra**(iq)
                  dqb=(2**(iq/3.)-2)*rb**(iq)
                  dqc=(2**(iq/3.)-2)*rc**(iq)
                  dqd=(5*0.4**(iq/3.)-2)*rd**(iq)
!
!-----Q7 MOMENT
!
                elseif(iq.eq.7)then
                  dqa=(2**(iq/3.)-2)*ra**(iq)
                  dqb=(2**(iq/3.)-2)*rb**(iq)
                  dqc=(2**(iq/3.)-2)*rc**(iq)
                  dqd=(5*0.4**(iq/3.)-2)*rd**(iq)
                endif
!
!-----COLLISION TERM
!
                term=collno*(pcoala*dqa+pcoalb*dqb+pcoalc*dqc+psep*dqd)
!
                su(i,j,k)=su(i,j,k)+term*vol(i,j,k)
                if(iq.eq.0)cq0(i,j,k)=term
                if(iq.eq.1)cq1(i,j,k)=term
                if(iq.eq.2)cq2(i,j,k)=term
                if(iq.eq.3)cq3(i,j,k)=term
                if(iq.eq.4)cq4(i,j,k)=term
                if(iq.eq.5)cq5(i,j,k)=term
                if(iq.eq.6)cq6(i,j,k)=term
                if(iq.eq.7)cq7(i,j,k)=term
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
            su(i,j,k)=su(i,j,k)+phi1(i,j,k)*phi2(i,j,k)*vol(i,j,k)
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
!
      call zerovar(temp)
      call zerovar(uvel)
      call zerovar(vvel)
      call zerovar(wvel)
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
!             elseif(iph.eq.2)then
!               if(iq.eq.0)then
!                 uvel(i,j,k)=uq0(i,j,k)
!                 vvel(i,j,k)=vq0(i,j,k)
!                 wvel(i,j,k)=wq0(i,j,k)
!               elseif(iq.eq.1)then
!                 uvel(i,j,k)=uq1(i,j,k)
!                 vvel(i,j,k)=vq1(i,j,k)
!                 wvel(i,j,k)=wq1(i,j,k)
!               elseif(iq.eq.2)then
!                 uvel(i,j,k)=uq2(i,j,k)
!                 vvel(i,j,k)=vq2(i,j,k)
!                 wvel(i,j,k)=wq2(i,j,k)
!               elseif(iq.eq.3)then
!                 uvel(i,j,k)=uq3(i,j,k)
!                 vvel(i,j,k)=vq3(i,j,k)
!                 wvel(i,j,k)=wq3(i,j,k)
!               elseif(iq.eq.4)then
!                 uvel(i,j,k)=uq4(i,j,k)
!                 vvel(i,j,k)=vq4(i,j,k)
!                 wvel(i,j,k)=wq4(i,j,k)
!               elseif(iq.eq.5)then
!                 uvel(i,j,k)=uq5(i,j,k)
!                 vvel(i,j,k)=vq5(i,j,k)
!                 wvel(i,j,k)=wq5(i,j,k)
!               elseif(iq.eq.6)then
!                 uvel(i,j,k)=uq6(i,j,k)
!                 vvel(i,j,k)=vq6(i,j,k)
!                 wvel(i,j,k)=wq6(i,j,k)
!               elseif(iq.eq.7)then
!                 uvel(i,j,k)=uq7(i,j,k)
!                 vvel(i,j,k)=vq7(i,j,k)
!                 wvel(i,j,k)=wq7(i,j,k)
!               endif
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
              
              term1=(term1b*(term1d/term1f)-term1a*(term1c/term1e))&
                /term1g
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
              
              term2=(term1b*(term1d/term1f)-term1a*(term1c/term1e))&
                /term1g
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
              
              term3=(term1b*(term1d/term1f)-term1a*(term1c/term1e))&
                /term1g
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
              
              term1=(term1b*(term1d/term1f)-term1a*(term1c/term1e))&
                /term1g
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
              
              term2=(term1b*(term1d/term1f)-term1a*(term1c/term1e))&
                /term1g
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
              
              term3=(term1b*(term1d/term1f)-term1a*(term1c/term1e))&
                /term1g
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
              
              term1=(term1b*(term1d/term1f)-term1a*(term1c/term1e))&
                /term1g
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
              
              term2=(term1b*(term1d/term1f)-term1a*(term1c/term1e))&
                /term1g
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
              
              term3=(term1b*(term1d/term1f)-term1a*(term1c/term1e))&
                /term1g
            endif
!
            term4=(term1+term2+term3)*vol(i,j,k)
            su(i,j,k)=su(i,j,k)+term4
          enddo
        enddo
      enddo
!
      return
      end