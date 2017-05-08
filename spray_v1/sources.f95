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
            qa=0;qb=0;terma=0;termb=0;term=0
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
                q3beta=(q3(i,j,k))**(1-rexpd)*(q4(i,j,k))**(rexpd)
                rdrag=q3beta/q3(i,j,k)
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
                if(qa.gt.tiny)terma=4.5*dvsg(i,j,k)*urel*qa
                if(qb.gt.tiny)termb=1.35*urel*(dng(i,j,k)*urelh*qb)**0.687 &
                  *((dvsg(i,j,k)*qa)/2.)**0.313
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
                if(qa.gt.tiny)terma=4.5*dvsg(i,j,k)*urel*qa
                if(qb.gt.tiny)termb=1.35*urel*(dng(i,j,k)*urelh)**0.687 &
                  *qb*(dvsg(i,j,k)/2.)**0.313
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
                if(qa.gt.tiny.and.qb.gt.tiny)terma=4.5*dvsg(i,j,k)*urel*(qa)**(1-d)*(qb)**(d)
                if(qc.gt.tiny.and.qd.gt.tiny)termb=1.35*urel*(dng(i,j,k)*urelh)**0.687 &
                  *(qc)**(1-d)*(qd)**(d)*(dvsg(i,j,k)/2.)**0.313
              endif
!
!-----SOURCE TERM
!
              term=(terma+termb)/dnl(i,j,k)
              if(iq.eq.3)term=term*(4/3.)*pi*dnl(i,j,k)
!
              if(iph.eq.1)su(i,j,k)=su(i,j,k)+term*vol(i,j,k)
              if(iph.eq.2)su(i,j,k)=su(i,j,k)-term*vol(i,j,k)
            endif
          enddo
        enddo
      enddo
      
      return
      end
!
!=====BREAK-UP  [HSIANG AND FAETH MODEL]
!
      subroutine breakup(rub,i,j,k)
      include "common.inc"
!
      bq0(i,j,k)=0;bq1(i,j,k)=0;bq2(i,j,k)=0
      bq4(i,j,k)=0;bq5(i,j,k)=0;bq6(i,j,k)=0;bq7(i,j,k)=0
!
!-----RELATIVE VELOCITY
!
      termu=uq3(i,j,k)-ug(i,j,k)
      termv=vq3(i,j,k)-vg(i,j,k)
      termw=wq3(i,j,k)-wg(i,j,k)
      urelh=sqrt((termu**2+termv**2+termw**2))
!
!-----INTEGRAL LIMITS
!
      if(urelh.gt.1.e-2)then
        rlb=6.*st(i,j,k)/(dng(i,j,k)*urelh**2)
!
!-----TRUNCATED MOMENTS
!
        call moment(q2t,1.,rlb,rub,i,j,k)
        call moment(q3t,2.,rlb,rub,i,j,k)
        if(q2t.gt.q2in*1.e-10.and.q3t.gt.q3in*1.e-10)then
          rin=q3t/q2t
        endif
        if(rin.gt.rlb.and.rin.lt.rub)then
!
!-----STABLE RADIUS
!
          rstab=(6.2/2)*(2*rin) &
            *(dnl(i,j,k)/dng(i,j,k))**0.25 &
            *sqrt(dvsl(i,j,k)/(dnl(i,j,k)*(2*rin)*urelh))
          rstabc=(6.2/2)*(2) &
            *(dnl(i,j,k)/dng(i,j,k))**0.25 &
            *sqrt(dvsl(i,j,k)/(dnl(i,j,k)*(2)*urelh))
!
!-----CHARACTERISTIC TIME SCALE
!
          tscale=5*(2*rin) &
            *sqrt(dnl(i,j,k)/dng(i,j,k)) &
            /(urelh*(1.-oh(rin,i,j,k)/7.))
          tscalec=5*(2) &
            *sqrt(dnl(i,j,k)/dng(i,j,k)) &
            /(urelh*(1.-oh(rin,i,j,k)/7.))
!
!-----SOURCE
!
          do iqq=0,iqmax-1
            if(rin.gt.rstab)then
!               power=(1+real(iqq))/2.
!               call moment(qa,power,rlb,rub,i,j,k)
!               power=real(iqq)-1.
!               call moment(qb,power,rlb,rub,i,j,k)
!               term=qa/(tscalec*rstabc**(3-iqq))-qb/tscalec
              term=(rin**3/rstab**(3-iqq)-rin**(iqq))/tscale
!               print*,' b:',iqq,rlb/rin,rin/rub,term
!
              if(iqq.eq.0)bq0(i,j,k)=term
              if(iqq.eq.1)bq1(i,j,k)=term
              if(iqq.eq.2)bq2(i,j,k)=term
              if(iqq.eq.4)bq4(i,j,k)=term
              if(iqq.eq.5)bq5(i,j,k)=term
              if(iqq.eq.6)bq6(i,j,k)=term
              if(iqq.eq.7)bq7(i,j,k)=term
            endif
          enddo
        endif
      endif
!
      return
      end
!
!=====COLLISION SOURCE TERM
!
      subroutine collisions(rub,qq0,i,j,k)
      include "common.inc"
!
      cq0(i,j,k)=0;cq1(i,j,k)=0;cq2(i,j,k)=0
      cq4(i,j,k)=0;cq5(i,j,k)=0;cq6(i,j,k)=0;cq7(i,j,k)=0
      rexpdc=0.!rexpd
!
!-----RELATIVE VELOCITY COEFFICIENTS
!
      q0beta=(q0(i,j,k))**(1-rexpdc)*(q1(i,j,k))**(rexpdc)
      termu=(uq0(i,j,k)-ug(i,j,k))*(q0(i,j,k)/q0beta)
      termv=(vq0(i,j,k)-vg(i,j,k))*(q0(i,j,k)/q0beta)
      termw=(wq0(i,j,k)-wg(i,j,k))*(q0(i,j,k)/q0beta)
      alpha0=sqrt(termu**2+termv**2+termw**2)
!
      q3beta=(q3(i,j,k))**(1-rexpdc)*(q4(i,j,k))**(rexpdc)
      termu=(uq3(i,j,k)-ug(i,j,k))*(q3(i,j,k)/q3beta)
      termv=(vq3(i,j,k)-vg(i,j,k))*(q3(i,j,k)/q3beta)
      termw=(wq3(i,j,k)-wg(i,j,k))*(q3(i,j,k)/q3beta)
      alpha3=sqrt(termu**2+termv**2+termw**2)
!
      alpha=alpha3-alpha0
!
      termu=uq3(i,j,k)-uq0(i,j,k)
      termv=vq3(i,j,k)-vq0(i,j,k)
      termw=wq3(i,j,k)-wq0(i,j,k)
      alpha=sqrt(termu**2+termv**2+termw**2)
!
!-----NUMBER OF COLLISIONS PER UNIT VOLUME PER UNIT TIME (*DENSITY)
!
      q1beta=q1(i,j,k)**(1-rexpdc)*q2(i,j,k)**(rexpdc)
      q2beta=q2(i,j,k)**(1-rexpdc)*q3(i,j,k)**(rexpdc)
      collno=0.15*pi*alpha*(q0(i,j,k)*q2beta+q1(i,j,k)*q1beta)
!
!-----CRITICAL RADII GIVEN BY CRITICAL WEBER NUMBER
!     ASSUMING THE CRITICAL WEBER NUMBERS ARE THOSE FOR DODECANE
!     NAMELY WEA=2.5, WEB=10, WEC=30
!
      if(collno.gt.tiny)then
        coeff=st(i,j,k)/(2*dng(i,j,k)*alpha**2)
        rlba=(2.5*coeff)**(1/(1+2.*rexpdc))
        rlbb=(10.*coeff)**(1/(1+2.*rexpdc))
        rlbc=(30.*coeff)**(1/(1+2.*rexpdc))
!
!-----AVERAGE RADII FOR EACH REGIME
!
        ra=rlba; rb=(rlbb+rlbc)/2; rc=rlbc; rd=rlbc
!
!-----PROBABILITY OF R<RCRITA
!
        call moment(ppa,0.,rlba,rub,i,j,k)
        proba=ppa/qq0
!
!-----PROBABILITY OF R<RCRITB
!
        call moment(ppb,0.,rlbb,rub,i,j,k)
        probb=ppb/qq0
!
!-----PROBABILITY OF R<RCRITC
!
        call moment(ppc,0.,rlbc,rub,i,j,k)
        probc=ppc/qq0
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
!-----SOURCE
!
        do iqq=0,iqmax-1
          dqa=0.;dqb=0.;dqc=0.;dqd=0.
          if(ra.gt.1.e-6.and.ra.lt.rub)dqa=(2.**(real(iqq)/3.)-2)*ra**(iqq)
          if(rb.gt.1.e-6.and.rb.lt.rub)dqb=(2.**(real(iqq)/3.)-2)*rb**(iqq)
          if(rc.gt.1.e-6.and.rc.lt.rub)dqc=(2.**(real(iqq)/3.)-2)*rc**(iqq)
          if(rd.gt.1.e-6.and.rd.lt.rub)dqd=(5*0.4**(iqq/3.)-2)*rd**(iqq)
          term=collno*(pcoala*dqa+pcoalb*dqb+pcoalc*dqc+psep*dqd)
!           if(iqq.eq.0.and.abs(term).gt.tiny)print*,' collisions'
!
          if(iqq.eq.0)cq0(i,j,k)=term
          if(iqq.eq.1)cq1(i,j,k)=term
          if(iqq.eq.2)cq2(i,j,k)=term
          if(iqq.eq.4)cq4(i,j,k)=term
          if(iqq.eq.5)cq5(i,j,k)=term
          if(iqq.eq.6)cq6(i,j,k)=term
          if(iqq.eq.7)cq7(i,j,k)=term
        enddo
      endif
!
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