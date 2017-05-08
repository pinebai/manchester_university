!
!=====INJECTION CONDITIONS
!
      subroutine ics(phi)
      include "common.inc"
      
      real phi(ia,ja)  ! NOT USED
      
      i=iinj
      do j=jinj,jinj+ninjc-1
!
!-----MOMENTS
!
        if(ivel.eq.0)then
          if(iq.eq.1)then
            qin=q1in
            su(i,j)=great*qin/1.e+20
            sp(i,j)=-great/1.e+20
          elseif(iq.eq.2)then
            qin=q2in
            su(i,j)=great*qin/1.e+25
            sp(i,j)=-great/1.e+25
          elseif(iq.eq.3)then
            qin=q3in
            su(i,j)=great*qin
            sp(i,j)=-great
          endif
!
!-----U - VELOCITIES
!
        elseif(ivel.eq.1)then
          if(iph.eq.1)uin=0.9*ulin
          if(iph.eq.2)uin=ulin
          su(i,j)=great*uin
          sp(i,j)=-great
!
!-----V - VELOCITIES
!
        elseif(ivel.eq.2)then
          if(iph.eq.1)vin=0.9*vlin
          if(iph.eq.2)vin=vlin
          su(i,j)=great*vin
          sp(i,j)=-great
        endif
      enddo
      
      return
      end