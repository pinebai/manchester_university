!
!=====MOMENT Q0
!
      subroutine calcq0
      include "common.inc"
!
      ivel=0;iph=2;iq=0
!
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=1.
            ff(i,j,k)=1.
          enddo
        enddo
      enddo
!
      call zerocoef
      call convection(uq0,vq0,wq0)
      call hybrid
      if(iquick.eq.1)call quick(q0)
      do ib=2,6
        call nbc(ib,q0)
      enddo
      call dbc(1,q0,0.)
      call injbc(q0)
      call spbc
      call temporal(q00)
      call solve(q0,q0in)
!
      return
      end
!
!=====MOMENT Q1
!
      subroutine calcq1
      include "common.inc"
!
      ivel=0;iph=2;iq=1
!
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=1.
            ff(i,j,k)=1.
          enddo
        enddo
      enddo
!
      call zerocoef
      call convection(uq1,vq1,wq1)
      call hybrid
      if(iquick.eq.1)call quick(q1)
      do ib=2,6
        call nbc(ib,q1)
      enddo
      call dbc(1,q1,0.)
      call injbc(q1)
      call spbc
      call temporal(q10)
      call solve(q1,q1in)
!
      return
      end
!
!=====MOMENT Q2
!
      subroutine calcq2
      include "common.inc"
!
      ivel=0;iph=2;iq=2
!
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=1.
            ff(i,j,k)=1.
          enddo
        enddo
      enddo
!
      call zerocoef
      call convection(uq2,vq2,wq2)
      call hybrid
      if(iquick.eq.1)call quick(q2)
      do ib=2,6
        call nbc(ib,q2)
      enddo
      call dbc(1,q2,0.)
      call injbc(q2)
      call spbc
      call temporal(q20)
      call solve(q2,q2in)
!
      return
      end
!
!=====MOMENT Q3
!
      subroutine calcq3
      include "common.inc"
!
      ivel=0;iph=2;iq=3
!
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=dnl0(i,j,k)
            ff(i,j,k)=dnl(i,j,k)
          enddo
        enddo
      enddo
!
      call zerocoef
      call convection(uq3,vq3,wq3)
      call hybrid
      if(iquick.eq.1)call quick(q3)
      do ib=2,6
        call nbc(ib,q3)
      enddo
      call dbc(1,q3,0.)
      call injbc(q3)
      call spbc
      call temporal(q30)
      call solve(q3,q3in)
!
      return
      end
!
!=====MOMENT Q4
!
      subroutine calcq4
      include "common.inc"
!
      ivel=0;iph=2;iq=4
!
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=1.
            ff(i,j,k)=1.
          enddo
        enddo
      enddo
!
      call zerocoef
      call convection(uq4,vq4,wq4)
      call hybrid
      if(iquick.eq.1)call quick(q4)
      do ib=2,6
        call nbc(ib,q4)
      enddo
      call dbc(1,q4,0.)
      call injbc(q4)
      call spbc
      call temporal(q40)
      call solve(q4,q4in)
!
      return
      end
!
!=====MOMENT Q5
!
      subroutine calcq5
      include "common.inc"
!
      ivel=0;iph=2;iq=5
!
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=1.
            ff(i,j,k)=1.
          enddo
        enddo
      enddo
!
      call zerocoef
      call convection(uq5,vq5,wq5)
      call hybrid
      if(iquick.eq.1)call quick(q5)
      do ib=2,6
        call nbc(ib,q5)
      enddo
      call dbc(1,q5,0.)
      call injbc(q5)
      call spbc
      call temporal(q50)
      call solve(q5,q5in)
!
      return
      end
!
!=====MOMENT Q6
!
      subroutine calcq6
      include "common.inc"
!
      ivel=0;iph=2;iq=6
!
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=1.
            ff(i,j,k)=1.
          enddo
        enddo
      enddo
!
      call zerocoef
      call convection(uq6,vq6,wq6)
      call hybrid
      if(iquick.eq.1)call quick(q6)
      do ib=2,6
        call nbc(ib,q6)
      enddo
      call dbc(1,q6,0.)
      call injbc(q6)
      call spbc
      call temporal(q60)
      call solve(q6,q6in)
!
      return
      end
!
!=====MOMENT Q7
!
      subroutine calcq7
      include "common.inc"
!
      ivel=0;iph=2;iq=7
!
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=1.
            ff(i,j,k)=1.
          enddo
        enddo
      enddo
!
      call zerocoef
      call convection(uq7,vq7,wq7)
      call hybrid
      if(iquick.eq.1)call quick(q7)
      do ib=2,6
        call nbc(ib,q7)
      enddo
      call dbc(1,q7,0.)
      call injbc(q7)
      call spbc
      call temporal(q70)
      call solve(q7,q7in)
!
      return
      end
!
!=====DERIVED MOMENTS
!
      subroutine calcqd
      include "common.inc"
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
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
!-----CHECK DISTRIBUTION
!
              rexpm=3.
              qcheck=dmom(0.,0.,i,j,k)/q3(i,j,k)
C               print*,'Q3 check =',qcheck
              if(qcheck.gt.tiny)then
!
!-----MOMENT RATIOS
!
                r32(i,j,k)=q3(i,j,k)/q2(i,j,k)
!
!-----MOMENTS REQUIRED FOR STANDARD DRAG MODEL
!
                if(rexpd.lt.tiny)then
                  do iqq=-1,-2,-1
                    rexpm=iqq
                    qa=dmom(0.,0.,i,j,k)
                    if(iqq.eq.-1)qm1(i,j,k)=qa
                    if(iqq.eq.-2)qm2(i,j,k)=qa
                  enddo
!
!-----FIRST MOMENT REQUIRED FOR ADVANCED DRAG MODEL
!
                elseif(rexpd.gt.tiny)then
                  do iqq=0,iqmax
                    rexpm=iqq-2.+rexpd
                    qa=dmom(0.,0.,i,j,k)
                    if(iqq.eq.0)qa0(i,j,k)=qa
                    if(iqq.eq.1)qa1(i,j,k)=qa
                    if(iqq.eq.2)qa2(i,j,k)=qa
                    if(iqq.eq.3)qa3(i,j,k)=qa
                    if(iqq.eq.4)qa4(i,j,k)=qa
                    if(iqq.eq.5)qa5(i,j,k)=qa
                    if(iqq.eq.6)qa6(i,j,k)=qa
                    if(iqq.eq.7)qa7(i,j,k)=qa
!
!-----SECOND MOMENT REQUIRED ADVANCED DRAG MODEL
!
                    rexpm=iqq-1.313+1.687*rexpd
                    qb=dmom(0.,0.,i,j,k)
                    if(iqq.eq.0)qb0(i,j,k)=qb
                    if(iqq.eq.1)qb1(i,j,k)=qb
                    if(iqq.eq.2)qb2(i,j,k)=qb
                    if(iqq.eq.3)qb3(i,j,k)=qb
                    if(iqq.eq.4)qb4(i,j,k)=qb
                    if(iqq.eq.5)qb5(i,j,k)=qb
                    if(iqq.eq.6)qb6(i,j,k)=qb
                    if(iqq.eq.7)qb7(i,j,k)=qb
                  enddo
                endif
              endif
!
!-----ZERO
!
            else
              r32(i,j,k)=0
!
              qa0(i,j,k)=0
              qa1(i,j,k)=0
              qa2(i,j,k)=0
              qa3(i,j,k)=0
              qa4(i,j,k)=0
              qa5(i,j,k)=0
!
              qb0(i,j,k)=0
              qb1(i,j,k)=0
              qb2(i,j,k)=0
              qb3(i,j,k)=0
              qb4(i,j,k)=0
              qb5(i,j,k)=0
!
              qm2(i,j,k)=0
              qm1(i,j,k)=0
              q0(i,j,k)=0
              q1(i,j,k)=0
              q2(i,j,k)=0
              q3(i,j,k)=0
              q4(i,j,k)=0
              q5(i,j,k)=0
              q6(i,j,k)=0
              q7(i,j,k)=0
            endif
          enddo
        enddo
      enddo
      
      return
      end