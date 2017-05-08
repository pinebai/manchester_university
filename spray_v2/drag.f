!
!=====INTER-PHASE DRAG
!
      subroutine drag(velg,vell)
      include "common.inc"
      real velg(ia,ja,ka),vell(ia,ja,ka)
      
      if(idrag.eq.0)goto 100
      order=1.e-8
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            term1=0;term2=0;term3=0;r32=0
            call dmom1(i,j,k)
            call dmom2(i,j,k)
            if(pk.gt.pkmax)pk=pkmax-tiny
            call dmom3(i,j,k)
!
!-----RELATIVE VELOCITY
!
            velu=ul3(i,j,k)-ug(i,j,k)
            velv=vl3(i,j,k)-vg(i,j,k)
            velw=wl3(i,j,k)-wg(i,j,k)
            vrel=(velu**2+velv**2+velw**2)**0.5
            vel=vell(i,j,k)-velg(i,j,k)
!
!=====FOUR MOMENT SCHEME
!
C             if(pp.gt.ppmin.and.pp.lt.ppmax
C      &        .and.pq.gt.pp.and.pq.lt.pqmax)then
C               qm1d=q0(i,j,k)*(pp+pq-1)/(rmax*(pp-1))
C               qm2d=qm1d*(pp+pq-2)/(rmax*(pp-2))
C !
C !-----Q0 MOMENTUM
C !
C               if(iq.eq.0)then
C                 term1=4.5*dvsg(i,j,k)*qm2d
C                 term2=1.35*(dvsg(i,j,k)/2)**0.313
C      &                *(dng(i,j,k)*vrel)**0.687
C      &                *(btf(pp-1.313,pq)/btf(pp,pq))
C      &                *rmax**(-1.313)
C !
C !-----Q1 MOMENTUM
C !
C               elseif(iq.eq.1)then
C                 term1=4.5*dvsg(i,j,k)*qm1d
C                 term2=1.35*(dvsg(i,j,k)/2)**0.313
C      &                *(dng(i,j,k)*vrel)**0.687
C      &                *(btf(pp-0.313,pq)/btf(pp,pq))
C      &                *rmax**(-0.313)
C !
C !-----Q2 MOMENTUM
C !
C               elseif(iq.eq.2)then
C                 term1=4.5*dvsg(i,j,k)*q0(i,j,k)
C                 term2=1.35*(dvsg(i,j,k)/2)**0.313
C      &                *(dng(i,j,k)*vrel)**0.687
C      &                *(btf(pp+0.687,pq)/btf(pp,pq))
C      &                *rmax**0.687
C !
C !-----Q3 MOMENTUM
C !
C               elseif(iq.eq.3)then
C                 term1=4.5*dvsg(i,j,k)*q1(i,j,k)
C                 term2=1.35*(dvsg(i,j,k)/2)**0.313
C      &                *(dng(i,j,k)*vrel)**0.687
C      &                *(btf(pp+1.687,pq)/btf(pp,pq))
C      &                *rmax**1.687
C               endif
!
!=====THREE MOMENT SCHEME
!
            if(pk.gt.pkmin.and.pk.lt.pkmax)then
              q0d=q1(i,j,k)*(pk+2)/(r32*pk)
              qm1d=q0d*(pk+2)/(r32*(pk-1))
!
!-----Q1 MOMENTUM
!
              if(iq.eq.1)then
                term1=4.5*dvsg(i,j,k)*qm1d
                term2=1.35*(dvsg(i,j,k)/2)**0.313
     &                *(dng(i,j,k)*vrel)**0.687
     &                *(gmf(pk-0.313)/gmf(pk))
     &                *(r32/(pk+2))**(-0.313)
!
!-----Q2 MOMENTUM
!
              elseif(iq.eq.2)then
                term1=4.5*dvsg(i,j,k)*q0d
                term2=1.35*(dvsg(i,j,k)/2)**0.313
     &                *(dng(i,j,k)*vrel)**0.687
     &                *(gmf(pk+0.687)/gmf(pk))
     &                *(r32/(pk+2))**0.687
!
!-----Q3 MOMENTUM
!
              elseif(iq.eq.3)then
                term1=4.5*dvsg(i,j,k)*q1(i,j,k)
                term2=1.35*(dvsg(i,j,k)/2)**0.313
     &                *(dng(i,j,k)*vrel)**0.687
     &                *(gmf(pk+1.687)/gmf(pk))
     &                *(r32/(pk+2))**1.687
              endif
!
!=====TWO MOMENT SCHEME
!
C             elseif(r32.gt.0.1*r32in.and.r32.lt.10.*r32in)then
C               q0d=q2(i,j,k)/r2
C               q1d=r1*q0d
C !
C !-----Q2 MOMENTUM
C !
C               if(iq.eq.2)then
C                 term1=4.5*dvsg(i,j,k)*q0d
C                 term2=1.35*(dvsg(i,j,k)*q0d/2)**0.313
C      &                *(dng(i,j,k)*q1d*vrel)**0.687
C !
C !-----Q3 MOMENTUM
C !
C               elseif(iq.eq.3)then
C                 term1=4.5*dvsg(i,j,k)*q1d
C                 term2=1.35*(dvsg(i,j,k)*q1d/2)**0.313
C      &                *(dng(i,j,k)*q2(i,j,k)*vrel)**0.687
C               endif
            endif
!
!-----INJECTION CELLS TREATMENT
!
            if(i.eq.imn.or.i.eq.imn+1)then
              if(j.eq.jinj.or.j.eq.jinj+1)then
                if(k.eq.kinj.or.k.eq.kinj+1)then
                  vel=0
                endif
              endif
            endif
!
!-----DRAG SOURCE TERM
!
            term3=(term1+term2)*vel*vol(i,j,k)
            if(ivel.eq.1.and.iq.eq.3)drgu(i,j,k)=term3
            if(ivel.eq.2.and.iq.eq.3)drgv(i,j,k)=term3
            if(ivel.eq.3.and.iq.eq.3)drgw(i,j,k)=term3
            if(iq.eq.3)term3=(4/3.)*pi*term3
            if(iph.eq.1)then
              su(i,j,k)=su(i,j,k)+term3
            elseif(iph.eq.2)then
              su(i,j,k)=su(i,j,k)-term3
            endif
          enddo
        enddo
      enddo
      
100   return
      end