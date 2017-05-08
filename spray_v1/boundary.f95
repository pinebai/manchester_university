!
!=====CELL FACE INTERPOLATED VALUES
!
      function f(ib,i,j,k,phi)
      include "common.inc"
      real phi(ia,ja,ka)
!
!-----WEST
!
      if(ib.eq.1)then
        f=fxw(i)*phi(i-1,j,k)+(1-fxw(i))*phi(i,j,k)
!
!-----EAST
!
      elseif(ib.eq.2)then
        f=(1-fxe(i))*phi(i,j,k)+fxe(i)*phi(i+1,j,k)
!
!-----SOUTH
!
      elseif(ib.eq.3)then
        f=fys(j)*phi(i,j-1,k)+(1-fys(j))*phi(i,j,k)
!
!-----NORTH
!
      elseif(ib.eq.4)then
        f=(1-fyn(j))*phi(i,j,k)+fyn(j)*phi(i,j+1,k)
!
!-----BOTTOM
!
      elseif(ib.eq.5)then
        f=fzb(k)*phi(i,j,k-1)+(1-fzb(k))*phi(i,j,k)
!
!-----TOP
!
      elseif(ib.eq.6)then
        f=(1-fzt(k))*phi(i,j,k)+fzt(k)*phi(i,j,k+1)
      endif

      return
      end
!
!=====CELL FACE TWO INTERPOLATED VALUES
!
      function f2(ib,i,j,k,phi1,phi2)
      include "common.inc"
      real phi1(ia,ja,ka),phi2(ia,ja,ka)
!
!-----WEST
!
      if(ib.eq.1)then
          f2=fxw(i)*phi1(i-1,j,k)*phi2(i-1,j,k)&
            +(1-fxw(i))*phi1(i,j,k)*phi2(i,j,k)
!
!-----EAST
!
      elseif(ib.eq.2)then
          f2=(1-fxe(i))*phi1(i,j,k)*phi2(i,j,k)&
            +fxe(i)*phi1(i+1,j,k)*phi2(i+1,j,k)
!
!-----SOUTH
!
      elseif(ib.eq.3)then
          f2=fys(j)*phi1(i,j-1,k)*phi2(i,j-1,k)&
            +(1-fys(j))*phi1(i,j,k)*phi2(i,j,k)
!
!-----NORTH
!
      elseif(ib.eq.4)then
          f2=(1-fyn(j))*phi1(i,j,k)*phi2(i,j,k)&
            +fyn(j)*phi1(i,j+1,k)*phi2(i,j+1,k)
!
!-----BOTTOM
!
      elseif(ib.eq.5)then
          f2=fzb(k)*phi1(i,j,k-1)*phi2(i,j,k-1)&
            +(1-fzb(k))*phi1(i,j,k)*phi2(i,j,k)
!
!-----TOP
!
      elseif(ib.eq.6)then
          f2=(1-fzt(k))*phi1(i,j,k)*phi2(i,j,k)&
            +fzt(k)*phi1(i,j,k+1)*phi2(i,j,k+1)
      endif

      return
      end
!
!=====CELL FACE THREE INTERPOLATED VALUES
!
      function f3(ib,i,j,k,phi1,phi2,phi3)
      include "common.inc"
      real phi1(ia,ja,ka),phi2(ia,ja,ka),phi3(ia,ja,ka)
!
!-----WEST
!
      if(ib.eq.1)then
        f3=fxw(i)*phi1(i-1,j,k)*phi2(i-1,j,k)*phi3(i-1,j,k)&
          +(1-fxw(i))*phi1(i,j,k)*phi2(i,j,k)*phi3(i,j,k)
!
!-----EAST
!
      elseif(ib.eq.2)then
        f3=(1-fxe(i))*phi1(i,j,k)*phi2(i,j,k)*phi3(i,j,k)&
          +fxe(i)*phi1(i+1,j,k)*phi2(i+1,j,k)*phi3(i+1,j,k)
!
!-----SOUTH
!
      elseif(ib.eq.3)then
        f3=fys(j)*phi1(i,j-1,k)*phi2(i,j-1,k)*phi3(i,j-1,k)&
          +(1-fys(j))*phi1(i,j,k)*phi2(i,j,k)*phi3(i,j,k)
!
!-----NORTH
!
      elseif(ib.eq.4)then
        f3=(1-fyn(j))*phi1(i,j,k)*phi2(i,j,k)*phi3(i,j,k)&
          +fyn(j)*phi1(i,j+1,k)*phi2(i,j+1,k)*phi3(i,j+1,k)
!
!-----BOTTOM
!
      elseif(ib.eq.5)then
        f3=fzb(k)*phi1(i,j,k-1)*phi2(i,j,k-1)*phi3(i,j,k-1)&
          +(1-fzb(k))*phi1(i,j,k)*phi2(i,j,k)*phi3(i,j,k)
!
!-----TOP
!
      elseif(ib.eq.6)then
        f3=(1-fzt(k))*phi1(i,j,k)*phi2(i,j,k)*phi3(i,j,k)&
          +fzt(k)*phi1(i,j,k+1)*phi2(i,j,k+1)*phi3(i,j,k+1)
      endif

      return
      end
!
!=====DIRICHLET CONDITIONS
!
      subroutine dbc(ib,phi,var)
      include "common.inc"
      real phi(ia,ja,ka)
!
!-----WEST
!
      if(ib.eq.1)then
        i=imn
        do j=jmn,jmx
          do k=kmn,kmx
            phi(i-1,j,k)=var
          enddo
        enddo
!
!-----EAST
!
      elseif(ib.eq.2)then
        i=imx
        do j=jmn,jmx
          do k=kmn,kmx
            phi(i+1,j,k)=var
          enddo
        enddo
!
!-----SOUTH
!
      elseif(ib.eq.3)then
        j=jmn
        do i=imn,imx
          do k=kmn,kmx
            phi(i,j-1,k)=var
          enddo
        enddo
!
!-----NORTH
!
      elseif(ib.eq.4)then
        j=jmx
        do i=imn,imx
          do k=kmn,kmx
            phi(i,j+1,k)=var
          enddo
        enddo
!
!-----BOTTOM
!
      elseif(ib.eq.5)then
        k=kmn
        do i=imn,imx
          do j=jmn,jmx
            phi(i,j,k-1)=var
          enddo
        enddo
!
!-----TOP
!
      elseif(ib.eq.6)then
        k=kmx
        do i=imn,imx
          do j=jmn,jmx
            phi(i,j,k+1)=var
          enddo
        enddo
      endif

      return
      end
!
!=====NEUMANN CONDITIONS
!
      subroutine nbc(ib,phi)
      include "common.inc"
      real phi(ia,ja,ka)
!
!-----WEST
!
      if(ib.eq.1)then
        i=imn
        do j=jmn,jmx
          do k=kmn,kmx
            phi(i-1,j,k)=phi(i,j,k)
          enddo
        enddo
!
!-----EAST
!
      elseif(ib.eq.2)then
        i=imx
        do j=jmn,jmx
          do k=kmn,kmx
            phi(i+1,j,k)=phi(i,j,k)
          enddo
        enddo
!
!-----SOUTH
!
      elseif(ib.eq.3)then
        j=jmn
        do i=imn,imx
          do k=kmn,kmx
            phi(i,j-1,k)=phi(i,j,k)
          enddo
        enddo
!
!-----NORTH
!
      elseif(ib.eq.4)then
        j=jmx
        do i=imn,imx
          do k=kmn,kmx
            phi(i,j+1,k)=phi(i,j,k)
          enddo
        enddo
!
!-----BOTTOM
!
      elseif(ib.eq.5)then
        k=kmn
        do i=imn,imx
          do j=jmn,jmx
            phi(i,j,k-1)=phi(i,j,k)
          enddo
        enddo
!
!-----TOP
!
      elseif(ib.eq.6)then
        k=kmx
        do i=imn,imx
          do j=jmn,jmx
            phi(i,j,k+1)=phi(i,j,k)
          enddo
        enddo
      endif

      return
      end
!
!=====INJECTION BOUNDARY CONDITIONS
!
      subroutine injbc(phi)
      include "common.inc"
      real phi(ia,ja,ka)
!
!-----INJECTOR REFERENCE CELL LOCATION
!
      i=imn-1
      j=jinj
      k=kinj
!
!-----INJECTOR VARIABLE
!
      if(ivel.eq.0)then
        if(iq.eq.0)varin=q0in
        if(iq.eq.1)varin=q1in
        if(iq.eq.2)varin=q2in
        if(iq.eq.3)varin=q3in
        if(iq.eq.4)varin=q4in
        if(iq.eq.5)varin=q5in
      elseif(ivel.gt.0)then
        if(ivel.eq.1)then
          if(iph.eq.1)varin=ugin
          if(iph.eq.2)varin=uqin*0.8
        elseif(ivel.eq.2)then
          if(iph.eq.1)varin=rugin
          if(iph.eq.2)varin=ruqin*0.8
        elseif(ivel.eq.3)then
          if(iph.eq.1)varin=rugin
          if(iph.eq.2)varin=ruqin*0.8
        endif
      endif
!
!-----CONTOURS
!
      const1=1;const2=1;const3=1
!
      if(nradic.eq.1)then
        const1=1
      elseif(nradic.eq.2)then
        if(ivel.eq.1)then
          const1=1
          const2=0.5
        elseif(ivel.eq.2.or.ivel.eq.3)then
          const1=0.5
          const2=1
        endif
      elseif(nradic.eq.3)then
        if(ivel.eq.1)then
          const1=0.999
          const2=0.666
          const3=0.333
        elseif(ivel.eq.2.or.ivel.eq.3)then
          const1=0.333
          const2=0.666
          const3=0.999
        endif
      endif
!
!-----1 RADIAL INJECTION CELL
!
      if(nradic.eq.1)then
        phi(i,j,k+1)=varin*comp(j,k+1)*const1
        phi(i,j,k)=varin*comp(j,k)*const1
!
        phi(i,j+1,k+1)=varin*comp(j+1,k+1)*const1
        phi(i,j+1,k)=varin*comp(j+1,k)*const1
!
!-----2 RADIAL INJECTION CELLS
!
      elseif(nradic.eq.2)then
        phi(i,j-1,k+1)=varin*comp(j-1,k+1)*const2
        phi(i,j-1,k)=varin*comp(j-1,k)*const2
!
        phi(i,j,k+2)=varin*comp(j,k+2)*const2
        phi(i,j,k+1)=varin*comp(j,k+1)*const1
        phi(i,j,k)=varin*comp(j,k)*const1
        phi(i,j,k-1)=varin*comp(j,k-1)*const2
!
        phi(i,j+1,k+2)=varin*comp(j+1,k+2)*const2
        phi(i,j+1,k+1)=varin*comp(j+1,k+1)*const1
        phi(i,j+1,k)=varin*comp(j+1,k)*const1
        phi(i,j+1,k-1)=varin*comp(j+1,k-1)*const2
!
        phi(i,j+2,k+1)=varin*comp(j+2,k+1)*const2
        phi(i,j+2,k)=varin*comp(j+2,k)*const2
!
!-----3 RADIAL INJECTION CELLS
!
      elseif(nradic.eq.3)then
        phi(i,j-2,k+1)=varin*comp(j-2,k+1)*const3
        phi(i,j-2,k)=varin*comp(j-2,k)*const3
!
        phi(i,j-1,k+2)=varin*comp(j-1,k+2)*const3
        phi(i,j-1,k+1)=varin*comp(j-1,k+1)*const2
        phi(i,j-1,k)=varin*comp(j-1,k)*const2
        phi(i,j-1,k-1)=varin*comp(j-1,k-1)*const3
!
        phi(i,j,k+3)=varin*comp(j,k+3)*const3
        phi(i,j,k+2)=varin*comp(j,k+2)*const2
        phi(i,j,k+1)=varin*comp(j,k+1)*const1
        phi(i,j,k)=varin*comp(j,k)*const1
        phi(i,j,k-1)=varin*comp(j,k-1)*const2
        phi(i,j,k-2)=varin*comp(j,k-2)*const3
!
        phi(i,j+1,k+3)=varin*comp(j+1,k+3)*const3
        phi(i,j+1,k+2)=varin*comp(j+1,k+2)*const2
        phi(i,j+1,k+1)=varin*comp(j+1,k+1)*const1
        phi(i,j+1,k)=varin*comp(j+1,k)*const1
        phi(i,j+1,k-1)=varin*comp(j+1,k-1)*const2
        phi(i,j+1,k-2)=varin*comp(j+1,k-2)*const3
!
        phi(i,j+2,k+2)=varin*comp(j+2,k+2)*const3
        phi(i,j+2,k+1)=varin*comp(j+2,k+1)*const2
        phi(i,j+2,k)=varin*comp(j+2,k)*const2
        phi(i,j+2,k-1)=varin*comp(j+2,k-1)*const3
!
        phi(i,j+3,k+1)=varin*comp(j+3,k+1)*const3
        phi(i,j+3,k)=varin*comp(j+3,k)*const3
      endif
!
      return
      end
!
!=====COMPONENT
!
      function comp(j,k)
      include "common.inc"
      real adj,opp,comp
!
      if(ivel.eq.0.or.ivel.eq.1)then
        comp=1.
        return
      endif
!
      opp=zc(k)-z(kinj)
      if(ivel.eq.2)opp=abs(opp)
!
      adj=yc(j)-y(jinj)
      if(ivel.eq.3)adj=abs(adj)
!
      comp=atan(opp/adj)
!
      return
      end
!
!=====SPRAY NOZZLE CONDITIONS
!
      subroutine spbc
      include "common.inc"
!
!-----INJECTOR REFERENCE CELL LOCATION
!
      i=imn
      j=jinj
      k=kinj
!
!-----1 RADIAL INJECTION CELL
!
      if(nradic.eq.1)then
        an(i,j-1,k+1)=0
        an(i,j-1,k)=0
!
        ab(i,j,k+2)=0
        at(i,j,k-1)=0
!
        ab(i,j+1,k+2)=0
        at(i,j+1,k-1)=0
!
        as(i,j+2,k+1)=0
        as(i,j+2,k)=0
!
!-----2 RADIAL INJECTION CELLS
!
      elseif(nradic.eq.2)then
        an(i,j-2,k+1)=0
        an(i,j-2,k)=0
!
        an(i,j-1,k+2)=0;ab(i,j-1,k+2)=0
        an(i,j-1,k-1)=0;at(i,j-1,k-1)=0
!
        ab(i,j,k+3)=0
        at(i,j,k-2)=0
!
        ab(i,j+1,k+3)=0
        at(i,j+1,k-2)=0
!
        as(i,j+2,k+2)=0;ab(i,j+2,k+2)=0
        as(i,j+2,k-1)=0;at(i,j+2,k-1)=0
!
        as(i,j+3,k+1)=0
        as(i,j+3,k)=0
!
!-----3 RADIAL INJECTION CELLS
!
      elseif(nradic.eq.3)then
        an(i,j-3,k+1)=0
        an(i,j-3,k)=0
!
        an(i,j-2,k+2)=0;ab(i,j-2,k+2)=0
        an(i,j-2,k-1)=0;at(i,j-2,k-1)=0
!
        an(i,j-1,k+3)=0;ab(i,j-1,k+3)=0
        an(i,j-1,k-2)=0;at(i,j-1,k-2)=0
!
        ab(i,j,k+4)=0
        at(i,j,k-3)=0
!
        ab(i,j+1,k+4)=0
        at(i,j+1,k-3)=0
!
        as(i,j+2,k+3)=0;ab(i,j+2,k+3)=0
        as(i,j+2,k-2)=0;at(i,j+2,k-2)=0
!
        as(i,j+3,k+2)=0;ab(i,j+3,k+2)=0
        as(i,j+3,k-1)=0;at(i,j+3,k-1)=0
!
        as(i,j+4,k+1)=0
        as(i,j+4,k)=0
      endif
!
      return
      end