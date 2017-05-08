!
!=====SOLUTION AND CONVERGENCE
!
      subroutine solve(phi,phiin)
      include "common.inc"
      real phi(ia,ja,ka),aa(ia,ja,ka),cc(ia,ja,ka)
      real temp(ia,ja,ka)

      iequ=2
!
!-----UNDER RELAXATION
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            if(abs(sp(i,j,k)).lt.1.e+15)then
              ap(i,j,k)=ap(i,j,k)/urf(iequ)
              su(i,j,k)=su(i,j,k)+(1.-urf(iequ))*ap(i,j,k)*phi(i,j,k)
            endif
          enddo
        enddo
      enddo   
!
!-----ITIALISE
!
      do iti=1,itimx
        call zerovar(aa)
        call zerovar(cc)
        resl=0
        do i=imn,imx
          do j=jmn,jmx
            do k=kmn,kmx
              if(abs(sp(i,j,k)).lt.1.e+15)then
                sum=-ap(i,j,k)*phi(i,j,k)&
                  +aw(i,j,k)*phi(i-1,j,k)&
                  +ae(i,j,k)*phi(i+1,j,k)&
                  +as(i,j,k)*phi(i,j-1,k)&
                  +an(i,j,k)*phi(i,j+1,k)&
                  +ab(i,j,k)*phi(i,j,k-1)&
                  +at(i,j,k)*phi(i,j,k+1)&
                  +su(i,j,k)
                endif
              resl=resl+abs(sum)
            enddo
          enddo
        enddo
        if(iti.eq.1)then
          res(iequ)=resl+tiny
        endif
!
!-----SOLVE
!
        do k=kmn,kmx
          do j=jmn,jmx
            do i=imn,imx
              term1=ap(i,j,k)-aw(i,j,k)*aa(i-1,j,k)+tiny
              term2=su(i,j,k)&
                +as(i,j,k)*phi(i,j-1,k)&
                +an(i,j,k)*phi(i,j+1,k)&
                +ab(i,j,k)*phi(i,j,k-1)&
                +at(i,j,k)*phi(i,j,k+1)
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
              term2=su(i,j,k)&
                +aw(i,j,k)*phi(i-1,j,k)&
                +ae(i,j,k)*phi(i+1,j,k)&
                +ab(i,j,k)*phi(i,j,k-1)&
                +at(i,j,k)*phi(i,j,k+1)
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
              term2=su(i,j,k)&
                +aw(i,j,k)*phi(i-1,j,k)&
                +ae(i,j,k)*phi(i+1,j,k)&
                +as(i,j,k)*phi(i,j-1,k)&
                +an(i,j,k)*phi(i,j+1,k)
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
        if(rsm.lt.0.2 .or.iti.eq.itimx)goto 100
      enddo
100   continue
!
!-----RESIDUAL
!
      call zerovar(temp)
      sum=0
      res(iequ)=0
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            if(abs(sp(i,j,k)).lt.1.e+15)then
              sum=-ap(i,j,k)*phi(i,j,k)&
                +aw(i,j,k)*phi(i-1,j,k)&
                +ae(i,j,k)*phi(i+1,j,k)&
                +as(i,j,k)*phi(i,j-1,k)&
                +an(i,j,k)*phi(i,j+1,k)&
                +ab(i,j,k)*phi(i,j,k-1)&
                +at(i,j,k)*phi(i,j,k+1)&
                +su(i,j,k)
            endif
            temp(i,j,k)=sum
            res(iequ)=res(iequ)+abs(sum)
          enddo
        enddo
      enddo
      res(iequ)=res(iequ)+tiny
!
!-----RESIDUAL NORMALIASTION
!
      i=imn
      j=jinj
      k=kinj
      if(iph.eq.1)fluxn=fluxn2
      if(iph.eq.2)fluxn=fluxn2
      resn(iequ)=fluxn*ff(i,j,k)*phiin+tiny
!
!-----MAGNITUDE, SIGN AND LOCATION OF MAXIMUM RESIDUAL
!
      resmx=0
      inew=0;jnew=0;knew=0
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            resmx0=resmx
            resmx=max(abs(temp(i,j,k)),abs(temp(i,j,k-1)))
!
            if(abs(temp(i,j,k)).gt.abs(temp(i,j,k-1)))then
              sign=temp(i,j,k)/abs(temp(i,j,k)+tiny)
            else
              sign=temp(i,j,k-1)/abs(temp(i,j,k-1)+tiny)
            endif
            resmx=sign*resmx
            if(abs(resmx).lt.abs(resmx0))resmx=resmx0
!
            iold=inew
            jold=jnew
            kold=knew
            if(resmx.gt.resmx0)then
              if(abs(temp(i,j,k)).gt.abs(temp(i,j,k-1)))then
                inew=i
                jnew=j
                knew=k
              endif
              if(abs(temp(i,j,k)).lt.abs(temp(i,j,k-1)))then
                inew=i
                jnew=j
                knew=k-1
              endif
            else
              inew=iold
              jnew=jold
              knew=kold
            endif
          enddo
        enddo
      enddo
!
!-----OUTPUT
!
!       write(6,'('' IEQU='',I3,'' RSM='',F6.3,'' I='',I3,'' J='',I3,
!      &  '' K='',I3,'' RESMX='',E10.3,'' |RESR|='',F6.3)')
!      &  iequ,rsm,inew,jnew,knew,resmx/resn(iequ),abs(resmx)/res(iequ)
      return
      end
!
!=====OUTER CONVERGRENCE
!
      subroutine converge
      include "common.inc"
10    format(1(e12.3))
20    format(3(e12.3))   
!
!=====GAS PHASE
!
      if(iph.eq.1)then
!
!-----PRINT RESIDUALS
!
        write(6,*)''
        write(6,*)'NORMALISED RESIDUALS:'
        write(6,*)'UG/VG/WG/FLUX'
        write(6,*)''
        write(6,10)res(1)/resn(1)
        write(6,10)res(6)/resn(6)
        write(6,10)res(11)/resn(11)
        write(6,10)res(16)/resn(16)
!
!-----SOURCE
!
        sgas=max(res(1)/resn(1),res(6)/resn(6)&
          ,res(11)/resn(11),res(16)/resn(16))
        write(6,*)''
        write(6,'('' GAS SOURCE = '',E9.3)')sgas
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
!=====LIQUID PHASE
!
      elseif(iph.eq.2)then
!
!-----PRINT RESIDUALS
!
        write(6,*)''
        write(6,*)'NORMALISED RESIDUALS:'
        write(6,*)'UL1,2,3/VL1,2,3/WL1,2,3/Q1,2,3'
        write(6,*)''
        write(6,20)res(3)/resn(3),res(4)/resn(4),res(5)/resn(5)
        write(6,20)res(8)/resn(8),res(9)/resn(9),res(10)/resn(10)
        write(6,20)res(13)/resn(13),res(14)/resn(14),res(15)/resn(15)
        write(6,20)res(18)/resn(18),res(19)/resn(19),res(20)/resn(20)
!
!-----SOURCE
!
        sliq=max(res(3)/resn(3),res(4)/resn(4),res(5)/resn(5)&
          ,res(8)/resn(8),res(9)/resn(9),res(10)/resn(10)&
          ,res(13)/resn(13),res(14)/resn(14),res(15)/resn(15)&
          ,res(18)/resn(18),res(19)/resn(19),res(20)/resn(20))
        write(6,*)''
        write(6,'('' LIQUID SOURCE = '',E9.3)')sliq
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
!
!=====BOTH PHASES
!
      elseif(iph.eq.0)then
        gsmax=max(sgas,sliq)
        write(6,*)''
        write(6,'('' TIME = '',E9.3)')time
        write(6,*)''
        write(6,'('' TIME STEP '',I4,'' OF '',I4)')nts,ntsmx
        write(6,'('' ITERATION '',I4,'' OF '',I4)')ito,itomx
        write(6,*)''
        write(6,'('' MAX SOURCE = '',E9.3)')gsmax
      endif
      write(6,*)'-----------------------------------------------'

      return
      end
      