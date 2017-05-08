!
!=====Q1 MOMENT
!
      subroutine calcq1
      include "common.inc"
      
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=dnl0(i,j,k)
            ff(i,j,k)=dnl(i,j,k)
          enddo
        enddo
      enddo
      ivel=0;iph=2;iq=1;iequ=18
      call zerocoef
      call convection(ul1,vl1,wl1)
      call hybrid
      if(iquick.eq.1)call quick(q1)
      do ib=2,6
        call nbc(ib,q1)
      enddo
      call dbc(1,q1,0.)
      call injbc2(q1,nradic)
      call spbc2(nradic)
      if(ibreak.eq.1)call breakup
      if(icoll.eq.1)call collisions
      call temporal(q10)
      call solve(q1)

      return
      end
!
!=====Q2 MOMENT
!
      subroutine calcq2
      include "common.inc"
      
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=dnl0(i,j,k)
            ff(i,j,k)=dnl(i,j,k)
          enddo
        enddo
      enddo
      ivel=0;iph=2;iq=2;iequ=19
      call zerocoef
      call convection(ul2,vl2,wl2)
      call hybrid
      if(iquick.eq.1)call quick(q2)
      do ib=2,6
        call nbc(ib,q2)
      enddo
      call dbc(1,q2,0.)
      call injbc2(q2,nradic)
      call spbc2(nradic)
      if(ibreak.eq.1)call breakup
      if(icoll.eq.1)call collisions
      call temporal(q20)
      call solve(q2)

      return
      end
!
!=====Q3 MOMENT
!
      subroutine calcq3
      include "common.inc"
      
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=dnl0(i,j,k)
            ff(i,j,k)=dnl(i,j,k)
          enddo
        enddo
      enddo
      ivel=0;iph=2;iq=3;iequ=20
      call zerocoef
      call convection(ul3,vl3,wl3)
      call hybrid
      if(iquick.eq.1)call quick(q3)
      do ib=2,6
        call nbc(ib,q3)
      enddo
      call dbc(1,q3,0.)
      call injbc2(q3,nradic)
      call spbc2(nradic)
      call temporal(q30)
      call solve(q3)
!
!-----DERIVED MOMENTS
!
      order=1.e-5
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            if(q1(i,j,k).gt.q1in*order
     &        .and.q2(i,j,k).gt.q2in*order
     &        .and.q3(i,j,k).gt.q3in*order)then
     
              r21(i,j,k)=q2(i,j,k)/(q1(i,j,k)+tiny)
              r32(i,j,k)=q3(i,j,k)/(q2(i,j,k)+tiny)
              term1=r21(i,j,k)/(r32(i,j,k)+tiny)
              term2=(1.-2.*term1)/(term1-1.+tiny)
              
              if(term2.lt.pkmin)then
                term2=pkmin+0.5
                q2(i,j,k)=(((term2+1)/(term2+2))
     &            *(q3(i,j,k)/(q1(i,j,k)+tiny)))**0.5
              endif
              if(term2.gt.pkmax)then
                term2=pkmin-0.5
                q1(i,j,k)=((term2+2)/(term2+1))
     &            *(q2(i,j,k)**2/(q3(i,j,k)+tiny))
              endif
              
              pk(i,j,k)=term2
              q0(i,j,k)=q1(i,j,k)*(pk(i,j,k)+2.)
     &          /(r32(i,j,k)*pk(i,j,k)+tiny)
              qm1(i,j,k)=q0(i,j,k)*(pk(i,j,k)+2.)
     &          /(r32(i,j,k)*(pk(i,j,k)-1.)+tiny)
            else
              r21(i,j,k)=0.
              r32(i,j,k)=0.
              pk(i,j,k)=0.
              qm1(i,j,k)=0.
              q0(i,j,k)=0.
              q1(i,j,k)=0.
              q2(i,j,k)=0.
              q3(i,j,k)=0.
            endif
          enddo
        enddo
      enddo
      
      return
      end