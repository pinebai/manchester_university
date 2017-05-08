!
!=====ANALITICAL GAMMA FUNCTION
!
      function gmf2(parak)
      
      pi=4.*atan(1.)
      gmf2=((parak/exp(1.)*(parak*sinh(1/parak
     &     +1./(810.*parak**6)))**0.5)**parak)
     &     *(2.*pi/parak)**0.5

      return
      end
!
!=====GAMMA FUNCTION
!
	  function gmf(parak)
      include "common.inc"

	  k=int(parak)
	  dparak=parak-k
	  prod=1
	  if(k.ge.2)then
	    do jj=1,k-1
	      prod=prod*(parak-jj)
	    enddo
	  elseif(k.eq.0)then
	    prod=prod/parak
	  endif
	  jdown=int(dparak/0.01+1)
	  dpar=dparak-jdown*0.01
	  jup=jdown+1
	  gmf=prod*(gf(jdown)+dpar*(gf(jup)-gf(jdown)))
      
      return
	  end
!
!=====BETA FUNCTION
!
	  function btf(parap,paraq)

	  gamp=gmf(parap)
      gamq=gmf(paraq)
      gampq=gmf(parap+paraq)
      btf=(gamp*gamq)/gampq
      
      return
	  end
!
!=====DERIVED MOMENTS FOR THE BETA DISTRIBUTION
!
      subroutine dmom1(i,j,k)
      include "common.inc"

      rmax=0;pp=0;pq=0
!
      if(q0(i,j,k).gt.q0in*order
     &  .and.q1(i,j,k).gt.q1in*order
     &  .and.q2(i,j,k).gt.q2in*order
     &  .and.q3(i,j,k).gt.q3in*order)then
!
        r10=q1(i,j,k)/(q0(i,j,k)+tiny)
        r21=q2(i,j,k)/(q1(i,j,k)+tiny)
        r32=q3(i,j,k)/(q2(i,j,k)+tiny)
!
        a=r21/(r32+tiny)
        b=r10/(r32+tiny)
!
        term1=a+a*b-2*b
        term2=2*a-b-1
        if(term1.gt.tiny.and.term2.gt.tiny)then
          pg=(2*a-b-1)/(term1+tiny)
          pp=2*b*(1-a)/(term1+tiny)
          pq=2*(a-b)*(1-a)*(1-b)/(term1*term2+tiny)
          if(pg.gt.tiny)rmax=r32/(pg+tiny)
        endif
      endif

      
      return
      end
!
!=====DERIVED MOMENTS FOR THE GAMMA DISTRIBUTION
!
      subroutine dmom2(i,j,k)
      include "common.inc"

      pk=0
!
      if(q1(i,j,k).gt.q1in*order
     &  .and.q2(i,j,k).gt.q2in*order
     &  .and.q3(i,j,k).gt.q3in*order)then
!
        r21=q2(i,j,k)/(q1(i,j,k)+tiny)
        r32=q3(i,j,k)/(q2(i,j,k)+tiny)
!
        ratio=r21/(r32+tiny)
        if(ratio.lt.0.999)pk=(1-2*ratio)/(ratio-1)
      endif
      
      return
      end
!
!=====DERIVED MOMENTS FOR BECK'S DISTRIBUTION
!
      subroutine dmom3(i,j,k)
      include "common.inc"

      r1=0;r2=0
!
      if(q2(i,j,k).gt.q2in*order
     &  .and.q3(i,j,k).gt.q3in*order)then
!
        r32=q3(i,j,k)/(q2(i,j,k)+tiny)
        rmx=100.e-6
        if(r32.lt.rmx)then
          if(r32.lt.0.01*r32in)r32=0.01*r32in
          a=16./r32**2
          b=-4./r32
!
!-----Q1/Q0 = R1
!
          term1=a*rmx**2*exp(b*rmx)/b
     &          -(2*a*rmx*exp(b*rmx)/b**2
     &          -2*a*exp(b*rmx)/b**3)
          term2=2*a/b**3
          r1=term1-term2
!
!-----Q2/Q0 = R2
!
          term1=a*rmx**3*exp(b*rmx)/b
     &          -(3*a*rmx**2*exp(b*rmx)/b**2
     &          -(6*a*rmx*exp(b*rmx)/b**3
     &          -6*a*exp(b*rmx)/b**4))
          term2=-6*a/b**4
          r2=term1-term2
        endif
      endif

      return
      end