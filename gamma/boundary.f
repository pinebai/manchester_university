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
          f2=fxw(i)*phi1(i-1,j,k)*phi2(i-1,j,k)
     &     +(1-fxw(i))*phi1(i,j,k)*phi2(i,j,k)
!
!-----EAST
!
      elseif(ib.eq.2)then
          f2=(1-fxe(i))*phi1(i,j,k)*phi2(i,j,k)
     &       +fxe(i)*phi1(i+1,j,k)*phi2(i+1,j,k)
!
!-----SOUTH
!
      elseif(ib.eq.3)then
          f2=fys(j)*phi1(i,j-1,k)*phi2(i,j-1,k)
     &       +(1-fys(j))*phi1(i,j,k)*phi2(i,j,k)
!
!-----NORTH
!
      elseif(ib.eq.4)then
          f2=(1-fyn(j))*phi1(i,j,k)*phi2(i,j,k)
     &       +fyn(j)*phi1(i,j+1,k)*phi2(i,j+1,k)
!
!-----BOTTOM
!
      elseif(ib.eq.5)then
          f2=fzb(k)*phi1(i,j,k-1)*phi2(i,j,k-1)
     &       +(1-fzb(k))*phi1(i,j,k)*phi2(i,j,k)
!
!-----TOP
!
      elseif(ib.eq.6)then
          f2=(1-fzt(k))*phi1(i,j,k)*phi2(i,j,k)
     &       +fzt(k)*phi1(i,j,k+1)*phi2(i,j,k+1)
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
        f3=fxw(i)*phi1(i-1,j,k)*phi2(i-1,j,k)*phi3(i-1,j,k)
     &     +(1-fxw(i))*phi1(i,j,k)*phi2(i,j,k)*phi3(i,j,k)
!
!-----EAST
!
      elseif(ib.eq.2)then
        f3=(1-fxe(i))*phi1(i,j,k)*phi2(i,j,k)*phi3(i,j,k)
     &     +fxe(i)*phi1(i+1,j,k)*phi2(i+1,j,k)*phi3(i+1,j,k)
!
!-----SOUTH
!
      elseif(ib.eq.3)then
        f3=fys(j)*phi1(i,j-1,k)*phi2(i,j-1,k)*phi3(i,j-1,k)
     &     +(1-fys(j))*phi1(i,j,k)*phi2(i,j,k)*phi3(i,j,k)
!
!-----NORTH
!
      elseif(ib.eq.4)then
        f3=(1-fyn(j))*phi1(i,j,k)*phi2(i,j,k)*phi3(i,j,k)
     &     +fyn(j)*phi1(i,j+1,k)*phi2(i,j+1,k)*phi3(i,j+1,k)
!
!-----BOTTOM
!
      elseif(ib.eq.5)then
        f3=fzb(k)*phi1(i,j,k-1)*phi2(i,j,k-1)*phi3(i,j,k-1)
     &     +(1-fzb(k))*phi1(i,j,k)*phi2(i,j,k)*phi3(i,j,k)
!
!-----TOP
!
      elseif(ib.eq.6)then
        f3=(1-fzt(k))*phi1(i,j,k)*phi2(i,j,k)*phi3(i,j,k)
     &     +fzt(k)*phi1(i,j,k+1)*phi2(i,j,k+1)*phi3(i,j,k+1)
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
      
      i=imn-1
      j=jinj
      k=kinj
!
!-----MOMENTS
!
      if(ivel.eq.0)then
        if(iq.eq.1)qin=q1in
        if(iq.eq.2)qin=q2in
        if(iq.eq.3)qin=q3in
        phi(i,j,k)=qin
        phi(i,j+1,k)=qin
        phi(i,j+1,k+1)=qin
        phi(i,j,k+1)=qin
!
!-----U - VELOCITIES
!
      elseif(ivel.eq.1)then
        if(iph.eq.2)uin=ulin
        phi(i,j,k)=uin
        phi(i,j+1,k)=uin
        phi(i,j+1,k+1)=uin
        phi(i,j,k+1)=uin
!
!-----V - VELOCITIES
!
      elseif(ivel.eq.2)then
        if(iph.eq.2)vin=vlin
        phi(i,j,k)=-vin
        phi(i,j+1,k)=vin
        phi(i,j+1,k+1)=vin
        phi(i,j,k+1)=-vin
!
!-----W - VELOCITIES
!
      elseif(ivel.eq.3)then
        if(iph.eq.2)win=wlin
        phi(i,j,k)=-win
        phi(i,j+1,k)=-win
        phi(i,j+1,k+1)=win
        phi(i,j,k+1)=win
      endif
      
      return
      end
!
!=====INJECTION BOUNDARY CONDITIONS
!
      subroutine injbc2(phi,n)
      include "common.inc"
      real phi(ia,ja,ka)
!
!-----INJECTOR REFERENCE CELL LOCATION
!
      i=imn-1
      j=jinj
      k=kinj
!
!=====1 RADIAL INJECTION CELL
!
      if(n.eq.1)then
!
!-----MOMENTS
!
        if(ivel.eq.0)then
          if(iq.eq.1)qin=q1in
          if(iq.eq.2)qin=q2in
          if(iq.eq.3)qin=q3in
          phi(i,j,k+1)=qin
          phi(i,j,k)=qin
      
          phi(i,j+1,k+1)=qin
          phi(i,j+1,k)=qin
!
!-----U - VELOCITIES
!
        elseif(ivel.eq.1)then
          if(iph.eq.1)uin=ugin
          if(iph.eq.2)uin=ulin
          phi(i,j,k+1)=uin
          phi(i,j,k)=uin
      
          phi(i,j+1,k+1)=uin
          phi(i,j+1,k)=uin
!
!-----V - VELOCITIES
!
        elseif(ivel.eq.2)then
          if(iph.eq.1)vin=vgin
          if(iph.eq.2)vin=vlin
          phi(i,j,k+1)=-vin
          phi(i,j,k)=-vin
      
          phi(i,j+1,k+1)=vin
          phi(i,j+1,k)=vin
!
!-----W - VELOCITIES
!
        elseif(ivel.eq.3)then
          if(iph.eq.1)win=wgin
          if(iph.eq.2)win=wlin
          phi(i,j,k+1)=win
          phi(i,j,k)=-win
      
          phi(i,j+1,k+1)=win
          phi(i,j+1,k)=-win
        endif
!
!=====2 RADIAL INJECTION CELLS
!
      elseif(n.eq.2)then
!
!-----MOMENTS
!
        if(ivel.eq.0)then
          if(iq.eq.1)qin=q1in
          if(iq.eq.2)qin=q2in
          if(iq.eq.3)qin=q3in
          phi(i,j-1,k+1)=qin
          phi(i,j-1,k)=qin
      
          phi(i,j,k+2)=qin
          phi(i,j,k+1)=qin
          phi(i,j,k)=qin
          phi(i,j,k-1)=qin
          
          phi(i,j+1,k+2)=qin
          phi(i,j+1,k+1)=qin
          phi(i,j+1,k)=qin
          phi(i,j+1,k-1)=qin
          
          phi(i,j+2,k+1)=qin
          phi(i,j+2,k)=qin
!
!-----U - VELOCITIES
!
        elseif(ivel.eq.1)then
          if(iph.eq.2)uin=ulin
          phi(i,j-1,k+1)=uin
          phi(i,j-1,k)=uin
      
          phi(i,j,k+2)=uin
          phi(i,j,k+1)=uin
          phi(i,j,k)=uin
          phi(i,j,k-1)=uin
          
          phi(i,j+1,k+2)=uin
          phi(i,j+1,k+1)=uin
          phi(i,j+1,k)=uin
          phi(i,j+1,k-1)=uin
          
          phi(i,j+2,k+1)=uin
          phi(i,j+2,k)=uin
!
!-----V - VELOCITIES
!
        elseif(ivel.eq.2)then
          if(iph.eq.2)vin=vlin
          phi(i,j-1,k+1)=-vin
          phi(i,j-1,k)=-vin
      
          phi(i,j,k+2)=-vin
          phi(i,j,k+1)=-vin*0.5
          phi(i,j,k)=-vin*0.5
          phi(i,j,k-1)=-vin
          
          phi(i,j+1,k+2)=vin
          phi(i,j+1,k+1)=vin*0.5
          phi(i,j+1,k)=vin*0.5
          phi(i,j+1,k-1)=vin
          
          phi(i,j+2,k+1)=vin
          phi(i,j+2,k)=vin
!
!-----W - VELOCITIES
!
        elseif(ivel.eq.3)then
          if(iph.eq.2)win=wlin
          phi(i,j-1,k+1)=win
          phi(i,j-1,k)=-win
      
          phi(i,j,k+2)=win
          phi(i,j,k+1)=win*0.5
          phi(i,j,k)=-win*0.5
          phi(i,j,k-1)=-win
          
          phi(i,j+1,k+2)=win
          phi(i,j+1,k+1)=win*0.5
          phi(i,j+1,k)=-win*0.5
          phi(i,j+1,k-1)=-win
          
          phi(i,j+2,k+1)=win
          phi(i,j+2,k)=-win
        endif
!
!=====3 RADIAL INJECTION CELLS
!
      elseif(n.eq.3)then
!
!-----MOMENTS
!
        if(ivel.eq.0)then
          if(iq.eq.1)qin=q1in
          if(iq.eq.2)qin=q2in
          if(iq.eq.3)qin=q3in
          do k=kinj,kinj+1
            phi(i,j-2,k)=qin
            phi(i,j+3,k)=qin
          enddo
          do k=kinj-1,kinj+2
            phi(i,j-1,k)=qin
            phi(i,j+2,k)=qin
          enddo
          do k=kinj-2,kinj+3
            phi(i,j,k)=qin
            phi(i,j+1,k)=qin
          enddo
!
!-----U - VELOCITIES
!
        elseif(ivel.eq.1)then
          if(iph.eq.2)uin=ulin
          do k=kinj,kinj+1
            phi(i,j-2,k)=uin
            phi(i,j+3,k)=uin
          enddo
          do k=kinj-1,kinj+2
            phi(i,j-1,k)=uin
            phi(i,j+2,k)=uin
          enddo
          do k=kinj-2,kinj+3
            phi(i,j,k)=uin
            phi(i,j+1,k)=uin
          enddo
!
!-----V - VELOCITIES
!
        elseif(ivel.eq.2)then
          if(iph.eq.2)vin=vlin
          phi(i,j-2,k+1)=-vin
          phi(i,j-2,k)=-vin
    
          phi(i,j-1,k+2)=-vin
          phi(i,j-1,k+1)=-vin*0.5!*0.666
          phi(i,j-1,k)=-vin*0.5!*0.666
          phi(i,j-1,k-1)=-vin
    
          phi(i,j,k+3)=-vin
          phi(i,j,k+2)=-vin*0.5!*0.666
          phi(i,j,k+1)=-vin*0.!333
          phi(i,j,k)=-vin*0.!333
          phi(i,j,k-1)=-vin*0.5!*0.666
          phi(i,j,k-2)=-vin
          
          phi(i,j+1,k+3)=vin
          phi(i,j+1,k+2)=vin*0.5!*0.666
          phi(i,j+1,k+1)=vin*0.!333
          phi(i,j+1,k)=vin*0.!333
          phi(i,j+1,k-1)=vin*0.5!*0.666
          phi(i,j+1,k-2)=vin
    
          phi(i,j+2,k+2)=vin
          phi(i,j+2,k+1)=vin*0.5!*0.666
          phi(i,j+2,k)=vin*0.5!*0.666
          phi(i,j+2,k-1)=vin
    
          phi(i,j+3,k+1)=vin
          phi(i,j+3,k)=vin
!
!-----W - VELOCITIES
!
        elseif(ivel.eq.3)then
          if(iph.eq.2)win=wlin
          phi(i,j-2,k+1)=win
          phi(i,j-2,k)=-win
    
          phi(i,j-1,k+2)=win
          phi(i,j-1,k+1)=win*0.5!0.666
          phi(i,j-1,k)=-win*0.5!*0.666
          phi(i,j-1,k-1)=-win
    
          phi(i,j,k+3)=win
          phi(i,j,k+2)=win*0.5!*0.666
          phi(i,j,k+1)=win*0.!333
          phi(i,j,k)=-win*0.!333
          phi(i,j,k-1)=-win*0.5!*0.666
          phi(i,j,k-2)=-win
          
          phi(i,j+1,k+3)=win
          phi(i,j+1,k+2)=win*0.5!*0.666
          phi(i,j+1,k+1)=win*0.!333
          phi(i,j+1,k)=-win*0.!333
          phi(i,j+1,k-1)=-win*0.5!*0.666
          phi(i,j+1,k-2)=-win
    
          phi(i,j+2,k+2)=win
          phi(i,j+2,k+1)=win*0.5!*0.666
          phi(i,j+2,k)=-win*0.5!*0.666
          phi(i,j+2,k-1)=-win
    
          phi(i,j+3,k+1)=win
          phi(i,j+3,k)=-win
        endif
      
      endif
      
      return
      end
!
!=====SPRAY NOZZLE CONDITIONS
!
      subroutine spbc2(n)
      include "common.inc"
!
!-----INJECTOR REFERENCE CELL LOCATION
!
      i=imn
      j=jinj
      k=kinj
!
!=====1 RADIAL INJECTION CELL
!
      if(n.eq.1)then
        ae(i,j-1,k+1)=0
        ae(i,j-1,k)=0
    
        as(i,j,k+2)=0
        as(i,j,k-1)=0
        
        an(i,j+1,k+2)=0
        an(i,j+1,k-1)=0
        
        aw(i,j+2,k+1)=0
        aw(i,j+2,k)=0
!
!=====2 RADIAL INJECTION CELLS
!
      elseif(n.eq.2)then
        ae(i,j-2,k+1)=0
        ae(i,j-2,k)=0
  
        ae(i,j-1,k+2)=0
        as(i,j-1,k+2)=0
        ae(i,j-1,k-1)=0
        an(i,j-1,k-1)=0

        as(i,j,k+3)=0
        an(i,j,k-2)=0
        
        as(i,j+1,k+3)=0
        an(i,j+1,k-2)=0
  
        aw(i,j+2,k+2)=0
        as(i,j+2,k+2)=0
        aw(i,j+2,k-1)=0
        an(i,j+2,k-1)=0
  
        aw(i,j+3,k+1)=0
        aw(i,j+3,k)=0
!
!=====3 RADIAL INJECTION CELLS
!
      elseif(n.eq.3)then
        ae(i,j-3,k+1)=0
        ae(i,j-3,k)=0
  
        ae(i,j-2,k+2)=0
        as(i,j-2,k+2)=0
        ae(i,j-2,k-1)=0
        an(i,j-2,k-1)=0
  
        ae(i,j-1,k+3)=0
        as(i,j-1,k+3)=0
        ae(i,j-1,k-2)=0
        an(i,j-1,k-2)=0

        as(i,j,k+4)=0
        an(i,j,k-3)=0
        
        as(i,j+1,k+4)=0
        an(i,j+1,k-3)=0
        
        aw(i,j+2,k+3)=0
        as(i,j+2,k+3)=0
        aw(i,j+2,k-2)=0
        an(i,j+2,k-2)=0
  
        aw(i,j+3,k+2)=0
        as(i,j+3,k+2)=0
        aw(i,j+3,k-1)=0
        an(i,j+3,k-1)=0
  
        aw(i,j+4,k+1)=0
        aw(i,j+4,k)=0
      endif
      
      return
      end