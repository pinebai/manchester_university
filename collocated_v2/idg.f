!
!=====INTER-PHASE DRAG
!
      subroutine idg(velg,vell)
      include "common.inc"
      
      real velg(ia,ja),vell(ia,ja)
      
      do i=imn,imx
        do j=jmn,jmx
!
!-----SOLVE DRAG ONLY IF K IS GREATER THAN CUT-OFF
!
          if(pk(i,j).gt.pkmin)then
!
!-----SOLVE VELOCITY
!
            if(iq.eq.1)then
              vrel=((ul1(i,j)-ug(i,j))**2
     &             +(vl1(i,j)-vg(i,j))**2)**0.5
            elseif(iq.eq.2)then
              vrel=((ul2(i,j)-ug(i,j))**2
     &             +(vl2(i,j)-vg(i,j))**2)**0.5
            elseif(iq.eq.3)then
              vrel=((ul3(i,j)-ug(i,j))**2
     &             +(vl3(i,j)-vg(i,j))**2)**0.5
            endif
!
!-----Q1 MOMENTUM
!
            if(iq.eq.1)then
              drg1a(i,j)=4.5*dvsg(i,j)*qm1(i,j)
              drg1b(i,j)=1.35*(dvsg(i,j)/2)**0.313
     &                   *(dng(i,j)*vrel)**0.687*q0(i,j)
     &                   *(gmf(pk(i,j)-0.313)/gmf(pk(i,j)))
     &                   *(r32(i,j)/(pk(i,j)+2))**(-0.313)
              drg(i,j)=(drg1a(i,j)+drg1b(i,j))
     &                 *(vell(i,j)-velg(i,j))*dvs(i,j)
              if(ivel.eq.1)drgu1(i,j)=drg(i,j)
              if(ivel.eq.2)drgv1(i,j)=drg(i,j)
!
!-----Q2 MOMENTUM
!
            elseif(iq.eq.2)then
              drg2a(i,j)=4.5*dvsg(i,j)*q0(i,j)
              drg2b(i,j)=1.35*(dvsg(i,j)/2)**0.313
     &                   *(dng(i,j)*vrel)**0.687*q0(i,j)
     &                   *(gmf(pk(i,j)+0.687)/gmf(pk(i,j)))
     &                   *(r32(i,j)/(pk(i,j)+2))**0.687
              drg(i,j)=(drg2a(i,j)+drg2b(i,j))
     &                 *(vell(i,j)-velg(i,j))*dvs(i,j)
              if(ivel.eq.1)drgu2(i,j)=drg(i,j)
              if(ivel.eq.2)drgv2(i,j)=drg(i,j)
!
!-----Q3 MOMENTUM
!
            elseif(iq.eq.3)then
              drg3a(i,j)=6*pi*dvsg(i,j)*q1(i,j)
              drg3b(i,j)=1.8*pi*(dvsg(i,j)/2)**0.313
     &                   *(dng(i,j)*vrel)**0.687*q0(i,j)
     &                   *(gmf(pk(i,j)+1.687)/gmf(pk(i,j)))
     &                   *(r32(i,j)/(pk(i,j)+2))**1.687
              drg(i,j)=(drg3a(i,j)+drg3b(i,j))
     &                 *(vell(i,j)-velg(i,j))*dvs(i,j)
              if(ivel.eq.1)drgu3(i,j)=drg(i,j)
              if(ivel.eq.2)drgv3(i,j)=drg(i,j)
            endif
          else
            drg(i,j)=0.
          endif
        enddo
      enddo
      
      return
      end