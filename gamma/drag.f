!
!=====INTER-PHASE DRAG
!
      subroutine drag(velg,vell)
      include "common.inc"
      real velg(ia,ja,ka),vell(ia,ja,ka)
      
      if(idrag.eq.0)goto 100
      order=1.e-4
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            term1=0;term2=0;term3=0
            call dmom1(i,j,k)
            call dmom2(i,j,k)
            call dmom3(i,j,k)
!
!-----RELATIVE VELOCITY
!
            if(iq.eq.0)then
              velu=ul0(i,j,k)-ug(i,j,k)
              velv=vl0(i,j,k)-vg(i,j,k)
              velw=wl0(i,j,k)-wg(i,j,k)
            elseif(iq.eq.1)then
              velu=ul1(i,j,k)-ug(i,j,k)
              velv=vl1(i,j,k)-vg(i,j,k)
              velw=wl1(i,j,k)-wg(i,j,k)
            elseif(iq.eq.2)then
              velu=ul2(i,j,k)-ug(i,j,k)
              velv=vl2(i,j,k)-vg(i,j,k)
              velw=wl2(i,j,k)-wg(i,j,k)
            elseif(iq.eq.3)then
              velu=ul3(i,j,k)-ug(i,j,k)
              velv=vl3(i,j,k)-vg(i,j,k)
              velw=wl3(i,j,k)-wg(i,j,k)
            endif
            velu=ul3(i,j,k)-ug(i,j,k)
            velv=vl3(i,j,k)-vg(i,j,k)
            velw=wl3(i,j,k)-wg(i,j,k)
            vrel=(velu**2+velv**2+velw**2)**0.5
            vel=vell(i,j,k)-velg(i,j,k)
!
!-----REYONLDS NUMBER
!
            re=2*dng(i,j,k)*vrel*r32/dvsg(i,j,k)
!
!=====THREE MOMENT SCHEME
!
            if(pk.gt.pkmin.and.pk.lt.pkmax)then
!
              if(re.lt.1000)then
                if(iq.eq.1)then
                  term1=4.5*dvsg(i,j,k)*qm1d
                elseif(iq.eq.2)then
                  term1=4.5*dvsg(i,j,k)*q0d
                elseif(iq.eq.3)then
                 term1=4.5*dvsg(i,j,k)*q1(i,j,k)
                endif
C                 term2=1.35*(dvsg(i,j,k)/2)**0.313
C      &                *(dng(i,j,k)*vrel)**0.687
C      &                *(gmf(pk+iq-1.313)/gmf(pk))
C      &                *(r32/(pk+2))**(iq-1.313)
                term2=1.35*(dvsg(i,j,k)/2)**0.313
     &                *(dng(i,j,k)*vrel)**0.687
     &                *exp(gmfln(pk+iq-1.313)-gmfln(pk))
     &                *(r32/(pk+2))**(iq-1.313)
!
              elseif(re.gt.1000)then
                if(iq.eq.1)then
                  term1=0.159*dng(i,j,k)*qm1d
                elseif(iq.eq.2)then
                  term1=0.159*dng(i,j,k)*q0d
                elseif(iq.eq.3)then
                 term1=0.159*dng(i,j,k)*q1(i,j,k)
                endif
                term2=0
              endif
            endif
!
!-----INJECTION CELLS TREATMENT
!
            if(i.eq.imn)then
              vel=0
            endif
!
!=====DRAG SOURCE TERM
!
            term3=(term1+term2)*vel*vol(i,j,k)
            if(iq.eq.3)term3=(4/3.)*pi*term3
            if(ivel.eq.1.and.iq.eq.3)drgu(i,j,k)=term3
            if(ivel.eq.2.and.iq.eq.3)drgv(i,j,k)=term3
            if(ivel.eq.3.and.iq.eq.3)drgw(i,j,k)=term3
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