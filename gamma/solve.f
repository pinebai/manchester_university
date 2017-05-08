!
!=====SOLUTION AND CONVERGENCE
!
      subroutine solve(phi)
      include "common.inc"
      real phi(ia,ja,ka),aa(ia,ja,ka),cc(ia,ja,ka)
!
!-----UNDER RELAXATION
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            if(abs(sp(i,j,k)).lt.great)then
              ap(i,j,k)=ap(i,j,k)/urf(iequ)
              su(i,j,k)=su(i,j,k)+(1-urf(iequ))*ap(i,j,k)*phi(i,j,k)
            endif
          enddo
        enddo
      enddo   
!
!-----RESIDUAL
!
      do iti=1,itimx
        call zerovar(aa)
        call zerovar(cc)
        resl=0
        do i=imn,imx
          do j=jmn,jmx
            do k=kmn,kmx
              sum=-ap(i,j,k)*phi(i,j,k)
     &            +aw(i,j,k)*phi(i-1,j,k)
     &            +ae(i,j,k)*phi(i+1,j,k)
     &            +as(i,j,k)*phi(i,j-1,k)
     &            +an(i,j,k)*phi(i,j+1,k)
     &            +ab(i,j,k)*phi(i,j,k-1)
     &            +at(i,j,k)*phi(i,j,k+1)
     &            +su(i,j,k)
              resl=resl+abs(sum)
            enddo
          enddo
        enddo
        if(iti.eq.1)res(iequ)=resl
!
!-----SOLVE
!
      do k=kmn,kmx
        do j=jmn,jmx
          do i=imn,imx
            term1=ap(i,j,k)-aw(i,j,k)*aa(i-1,j,k)+tiny
            term2=su(i,j,k)
     &            +as(i,j,k)*phi(i,j-1,k)
     &            +an(i,j,k)*phi(i,j+1,k)
     &            +ab(i,j,k)*phi(i,j,k-1)
     &            +at(i,j,k)*phi(i,j,k+1)
            aa(i,j,k)=ae(i,j,k)/term1
            cc(i,j,k)=(aw(i,j,k)*cc(i-1,j,k)+term2)/term1
          enddo
        enddo
      enddo
      do k=kmn,kmx
        do j=jmn,jmx
          do i=imx,imn,-1
            phi(i,j,k)=phi(i+1,j,k)*aa(i,j,k)+cc(i,j,k)
          enddo
        enddo
      enddo
!
      do i=imn,imx
        do k=kmn,kmx
          do j=jmn,jmx
            term1=ap(i,j,k)-as(i,j,k)*aa(i,j-1,k)+tiny
            term2=su(i,j,k)
     &            +aw(i,j,k)*phi(i-1,j,k)
     &            +ae(i,j,k)*phi(i+1,j,k)
     &            +ab(i,j,k)*phi(i,j,k-1)
     &            +at(i,j,k)*phi(i,j,k+1)
            aa(i,j,k)=an(i,j,k)/term1
            cc(i,j,k)=(as(i,j,k)*cc(i,j-1,k)+term2)/term1
          enddo
        enddo
      enddo
      do i=imn,imx
        do k=kmn,kmx
          do j=jmx,jmn,-1
            phi(i,j,k)=phi(i,j+1,k)*aa(i,j,k)+cc(i,j,k)
          enddo
        enddo
      enddo
!
      do j=jmn,jmx
        do i=imn,imx
          do k=kmn,kmx
            term1=ap(i,j,k)-ab(i,j,k)*aa(i,j,k-1)+tiny
            term2=su(i,j,k)
     &            +aw(i,j,k)*phi(i-1,j,k)
     &            +ae(i,j,k)*phi(i+1,j,k)
     &            +as(i,j,k)*phi(i,j-1,k)
     &            +an(i,j,k)*phi(i,j+1,k)
            aa(i,j,k)=at(i,j,k)/term1
            cc(i,j,k)=(ab(i,j,k)*cc(i,j,k-1)+term2)/term1
          enddo
        enddo
      enddo
      do j=jmn,jmx
        do i=imn,imx
          do k=kmx,kmn,-1
            phi(i,j,k)=phi(i,j,k+1)*aa(i,j,k)+cc(i,j,k)
          enddo
        enddo
      enddo
!
!-----INNER CONVERGENCE
!
        rsm=resl/(res(iequ)+tiny)
        if(rsm.lt.0.2.or.iti.eq.itimx)then
          write(6,'('' IEQU ='',I3,''  RSM ='',F6.3,''  ITI ='',I3)')
     &      iequ,rsm,iti
          goto 100
        endif
      enddo
100   continue
      
      return
      end
!
!=====OUTER CONVERGRENCE
!
      subroutine converge(iphase)
      include "common.inc"
10    format(1(e12.3))
20    format(3(e12.3))
!
!-----NORMALISE ALL RESIDUALS
!
      i=imn;j=jinj;k=kinj
      flux1=dng(i,j,k)*pi*rs**2*ugin
      flux2=dnl(i,j,k)*pi*rs**2*ulin
      
      res(1)=res(1)/(flux1*(1-cnt*q3in)*ugin)
      res(3)=res(3)/(flux2*q1in*ulin)
      res(4)=res(4)/(flux2*q2in*ulin)
      res(5)=res(5)/(flux2*cnt*q3in*ulin)
      
      res(6)=res(6)/(flux1*(1-cnt*q3in)*vgin)
      res(8)=res(8)/(flux2*q1in*vlin)
      res(9)=res(9)/(flux2*q2in*vlin)
      res(10)=res(10)/(flux2*cnt*q3in*vlin)
      
      res(11)=res(11)/(flux1*(1-cnt*q3in)*wgin)
      res(13)=res(13)/(flux2*q1in*wlin)
      res(14)=res(14)/(flux2*q2in*wlin)
      res(15)=res(15)/(flux2*cnt*q3in*wlin)
      
      res(16)=res(16)/1
      
      res(18)=res(18)/(flux2*q1in)
      res(19)=res(19)/(flux2*q2in)
      res(20)=res(20)/(flux2*q3in)    
!
!=====GAS PHASE
!
      if(iphase.eq.1)then
!
!-----MONITOR
!
        i=imon;j=jmon;k=kmon
        write(6,*)''
        write(6,*)'MONITOR:'
        write(6,*)'UG/VG/WG/P'
        write(6,*)''
        write(6,10)ug(i,j,k)
        write(6,10)vg(i,j,k)
        write(6,10)wg(i,j,k)
        write(6,10)p(i,j,k)
!
!-----PRINT RESIDUALS
!
        write(6,*)''
        write(6,*)'NORMALISED RESIDUALS:'
        write(6,*)'UG/VG/WG/FLUX'
        write(6,*)''
        write(6,10)res(1)
        write(6,10)res(6)
        write(6,10)res(11)
        write(6,10)res(16)
!
!-----SOURCE
!
        sgas=max(res(1),res(6),res(11),res(16))
        write(6,*)''
        write(6,'('' GAS SOURCE = '',E9.3)')sgas
!
!=====LIQUID PHASE
!
      elseif(iphase.eq.2)then
!
!-----MONITOR
!
        i=imon;j=jmon;k=kmon
        write(6,*)''
        write(6,*)'MONITOR:'
        write(6,*)'UL1,2,3/VL1,2,3/WL1,2,3/Q1,2,3'
        write(6,*)''
        write(6,20)ul1(i,j,k),ul2(i,j,k),ul3(i,j,k)
        write(6,20)vl1(i,j,k),vl2(i,j,k),vl3(i,j,k)
        write(6,20)wl1(i,j,k),wl2(i,j,k),wl3(i,j,k)
        write(6,20)q1(i,j,k),q2(i,j,k),q3(i,j,k)
        write(6,*)''
!
!-----PRINT RESIDUALS
!
        write(6,*)''
        write(6,*)'NORMALISED RESIDUALS:'
        write(6,*)'UL1,2,3/VL1,2,3/WL1,2,3/Q1,2,3'
        write(6,*)''
        write(6,20)res(3),res(4),res(5)
        write(6,20)res(8),res(9),res(10)
        write(6,20)res(13),res(14),res(15)
        write(6,20)res(18),res(19),res(20)
!
!-----SOURCE
!
        sliq=max(res(2),res(3),res(4),res(5)
     &    ,res(7),res(8),res(9),res(10)
     &    ,res(12),res(13),res(14),res(15)
     &    ,res(17),res(18),res(19),res(20))
        write(6,*)''
        write(6,'('' LIQUID SOURCE = '',E9.3)')sliq
!
!=====BOTH PHASES
!
      elseif(iphase.eq.0)then
        gsmax0=gsmax
        gsmax=max(sgas,sliq)
        gdelta=abs(gsmax-gsmax0)
        write(6,*)''
        write(6,'('' TIME = '',E9.3)')time
        write(6,'('' TIME STEP '',I4,'' OF '',I4)')nts,ntsmx
        write(6,'('' ITERATION '',I4,'' OF '',I4)')ito,itomx
        write(6,*)''
        write(6,'(''            DELTA = '',E9.3)')gdelta
      endif
      write(6,*)'__________________________________________________'
      write(6,*)'__________________________________________________'

      return
      end
      