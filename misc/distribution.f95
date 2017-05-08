!
!=====MAIN
!
      program main
      common/moments/q0(1,1,1),q1(1,1,1),q2(1,1,1),q3(1,1,1)&
        ,q4(1,1,1),q5(1,1,1),q6(1,1,1),q7(1,1,1)&
        ,iqmax,rexpm
      real qr(0:7)
!
      do iqq=0,7
        qr(iqq)=0.
      enddo
!
!       write(*,*)'iqmax,r32in,pkin'
!       read*,iqmax
!       read*,r32in
!       read*,pkin
!
      iqmax=6
      r32in=15.e-6
      pkin=2.
!
      i=1;j=1;k=1
      q3(i,j,k)=0.145
!
!-----MOMENTS AT INLET BASED ON GAMMA DISTRIBUTION
!
      q2(i,j,k)=q3(i,j,k)/r32in
      q1(i,j,k)=q2(i,j,k)*(pkin+2)/(r32in*(pkin+1))
      q0(i,j,k)=q1(i,j,k)*(pkin+2)/(r32in*pkin)
      q4(i,j,k)=((pkin+3)/(pkin+2))*r32in*q3(i,j,k)
      q5(i,j,k)=((pkin+4)/(pkin+2))*r32in*q4(i,j,k)
      q6(i,j,k)=((pkin+5)/(pkin+2))*r32in*q5(i,j,k)
      q7(i,j,k)=((pkin+6)/(pkin+2))*r32in*q6(i,j,k)

      rexpm=-1.
      qa=dmom(0.,0.,i,j,k)
      print*,q0(i,j,k),qa
!
      end
!
!=====MOMENTS OF THE SIZE DISTRIBUTION
!
      function dmom(rlb,rub,i,j,k)
      common/moments/q0(1,1,1),q1(1,1,1),q2(1,1,1),q3(1,1,1)&
        ,q4(1,1,1),q5(1,1,1),q6(1,1,1),q7(1,1,1)&
        ,iqmax,rexpm
      common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
      external flagu
!
!-----SIZE DISTRIBUTION NORMALISATION RADIUS
!
      if(iqmax.eq.4)rnorm=(q3(i,j,k)/q0(i,j,k))**(1/3.)
      if(iqmax.eq.5)rnorm=(q4(i,j,k)/q0(i,j,k))**(1/4.)
      if(iqmax.eq.6)rnorm=(q5(i,j,k)/q0(i,j,k))**(1/5.)
      if(iqmax.eq.7)rnorm=(q6(i,j,k)/q0(i,j,k))**(1/6.)
      if(iqmax.eq.8)rnorm=(q7(i,j,k)/q0(i,j,k))**(1/7.)
!
!-----NORMALISED MOMENTS
!
      a0=(q0(i,j,k)/q0(i,j,k))*(1./rnorm)**0
      a1=(q1(i,j,k)/q0(i,j,k))*(1./rnorm)**1
      a2=(q2(i,j,k)/q0(i,j,k))*(1./rnorm)**2
      if(iqmax.ge.4)a3=(q3(i,j,k)/q0(i,j,k))*(1./rnorm)**3
      if(iqmax.ge.5)a4=(q4(i,j,k)/q0(i,j,k))*(1./rnorm)**4
      if(iqmax.ge.6)a5=(q5(i,j,k)/q0(i,j,k))*(1./rnorm)**5
      if(iqmax.ge.7)a6=(q6(i,j,k)/q0(i,j,k))*(1./rnorm)**6
      if(iqmax.ge.8)a7=(q7(i,j,k)/q0(i,j,k))*(1./rnorm)**7
!
!-----LAGUERRE COEFFICIENTS
!
      b0=plagu(0,a0)
      b1=plagu(1,a1)
      b2=plagu(2,a2)
      if(iqmax.ge.4)b3=plagu(3,a3)
      if(iqmax.ge.5)b4=plagu(4,a4)
      if(iqmax.ge.6)b5=plagu(5,a5)
      if(iqmax.ge.7)b6=plagu(6,a6)
      if(iqmax.ge.8)b7=plagu(7,a7)
!
!-----LIMITS OF DISTRIBUTION FUNCTION
!
      delr=1.e-7;rrlb=1.e-6;rrub=2*q3(i,j,k)/q2(i,j,k)
!
      if(rlb.lt.1.e-7)then
        root1=0
        do r=rrlb,rrub,delr
          terml=flagu((r-delr)/rnorm)
          termc=flagu(r/rnorm)
          termu=flagu((r+delr)/rnorm)
          if(termc.gt.1.e-10.and.terml.lt.termu)then
            root1=r
            goto 100
          endif
        enddo
100     continue
      else
        if(rlb.gt.rrub)rlb=rrub
        root1=rlb
      endif
!
      if(rub.lt.1.e-7)then
        root2=0
        do r=rrub,rrlb,-delr
          terml=flagu((r-delr)/rnorm)
          termc=flagu(r/rnorm)
          termu=flagu((r+delr)/rnorm)
          if(termc.gt.1.e-10.and.terml.gt.termu)then
            root2=r
            goto 200
          endif
        enddo
200     continue
      else
        if(rub.gt.rrub)rlb=rrub
        root2=rub
      endif
!
      dmom=0.
      if(root1.gt.1.e-10.and.root2.gt.root1)then
        rlbn=root1/rnorm
        rubn=root2/rnorm
!
!-----NUMERICAL INTEGRATION
!
!         if(1+rexpm.lt.0.5)then
!         call qtrap(flagu,rlbn,rubn,qnum)
!
!-----ANALYTICAL INTEGRATION
!
!         else
          call qalt(rlbn,rubn,qnum)
!         endif
!
!-----UN-NORMALISE MOMENT
!
        qnorm=q0(i,j,k)/(1./rnorm)**(rexpm)
        dmom=qnum*qnorm
      endif
!
      return
      end
!
!=====LAGUERRE POLYNOMIALS
!
      function plagu(iqq,var)
!
      plagu=0.
      if(iqq.eq.0)then
        plagu=1.
      elseif(iqq.eq.1)then
        plagu=1.-var
      elseif(iqq.eq.2)then
        plagu=(1/2.)&
        *(2.-4.*var+var**2)
      elseif(iqq.eq.3)then
        plagu=(1/6.)&
        *(6.-18.*var+9.*var**2&
        -var**3)
      elseif(iqq.eq.4)then
        plagu=(1/24.)&
        *(24.-96.*var+72.*var**2&
        -16.*var**3+var**4)
      elseif(iqq.eq.5)then
        plagu=(1/120.)&
        *(120.-600.*var+600.*var**2&
        -200.*var**3+25.*var**4-var**5)
      elseif(iqq.eq.6)then
        plagu=(1/720.)&
        *(720.-4320.*var+5400.*var**2&
        -2400.*var**3+450.*var**4-36.*var**5&
        +var**6)
      elseif(iqq.eq.7)then
        plagu=(1/5040.)&
        *(5040.-35280.*var+52920.*var**2&
        -29400.*var**3+7350.*var**4-882.*var**5&
        +49.*var**6-var**7)
      endif
!
      return
      end
!
!=====SIZE DISTRIBUTION FUNCTION
!
      function flagu(xx)
      common/moments/q0(1,1,1),q1(1,1,1),q2(1,1,1),q3(1,1,1)&
        ,q4(1,1,1),q5(1,1,1),q6(1,1,1),q7(1,1,1)&
        ,iqmax,rexpm
      common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
!
      flagu=0.
      if(iqmax.eq.4)then
        flagu=exp(-xx)*xx**(rexpm)&
          *(b0*plagu(0,xx)&
          +b1*plagu(1,xx)&
          +b2*plagu(2,xx)&
          +b3*plagu(3,xx))
      elseif(iqmax.eq.5)then
        flagu=exp(-xx)*xx**(rexpm)&
          *(b0*plagu(0,xx)&
          +b1*plagu(1,xx)&
          +b2*plagu(2,xx)&
          +b3*plagu(3,xx)&
          +b4*plagu(4,xx))
      elseif(iqmax.eq.6)then
        flagu=exp(-xx)*xx**(rexpm)&
          *(b0*plagu(0,xx)&
          +b1*plagu(1,xx)&
          +b2*plagu(2,xx)&
          +b3*plagu(3,xx)&
          +b4*plagu(4,xx)&
          +b5*plagu(5,xx))
      elseif(iqmax.eq.7)then
        flagu=exp(-xx)*xx**(rexpm)&
          *(b0*plagu(0,xx)&
          +b1*plagu(1,xx)&
          +b2*plagu(2,xx)&
          +b3*plagu(3,xx)&
          +b4*plagu(4,xx)&
          +b5*plagu(5,xx)&
          +b6*plagu(6,xx))
      elseif(iqmax.eq.8)then
        flagu=exp(-xx)*xx**(rexpm)&
          *(b0*plagu(0,xx)&
          +b1*plagu(1,xx)&
          +b2*plagu(2,xx)&
          +b3*plagu(3,xx)&
          +b4*plagu(4,xx)&
          +b5*plagu(5,xx)&
          +b6*plagu(6,xx)&
          +b7*plagu(7,xx))
      endif
!       if(flagu.lt.0.)flagu=0.
!
      return
      end
!
!=====INTEGRATION REFINEMENT [Num. Res. in F77, Ch 4.2]
!
      subroutine qtrap(func,a,b,s)
!
      integer jmax
      real a,b,func,s,eps
      external func
      parameter (eps=1.e-6,jmax=5)
!
      integer j
      real olds
!
      olds=-1.e30
      do j=1,jmax
        call trapzd(func,a,b,s,j)
        if(j.gt.5)then
          if(abs(s-olds).lt.eps*abs(olds).or.&
            (s.eq.0. .and.olds.eq.0))return
        endif
        olds=s
      enddo
!       print*, 'too many steps in qtrap'
!
      return
      end
!
!=====NUMERICAL INTEGRATION
!     USING EXTENDED TRAPEZOIDAL RULE [Num. Res. in F77, Ch 4.2]
!
      subroutine trapzd(func,a,b,s,n)
!
      integer n
      real a,b,s,func
      external func
!
      integer it,j
      real del,sum,tnm,x
!
      if(n.eq.1)then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do j=1,it
          sum=sum+func(x)
          x=x+del
        enddo
        s=0.5*(s+(b-a)*sum/tnm)
      endif
!
      return
      end
!
!=====INTEGRAL FORM OF LAGUERRE FUNCTION
!
!       function filagu(x)
!       common/moments/q0(1,1,1),q1(1,1,1),q2(1,1,1),q3(1,1,1)&
!         ,q4(1,1,1),q5(1,1,1),q6(1,1,1),q7(1,1,1)&
!         ,iqmax,power
!       common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
! !
!       filagu=0.
!       filagu=0. &
!       
!         +(b0)*(-gammq(1+power,x)) &
!         
!         +(b1)*(-gammq(1+power,x)+gammq(2+power,x)) &
!       
!         +(b2/2)*(-2*gammq(1+power,x)+4*gammq(2+power,x)-gammq(3+power,x)) &
!         
!         +(b3/6)*(-6*gammq(1+power,x)+18*gammq(2+power,x)-9*gammq(3+power,x)+gammq(4+power,x))
! !
!       if(iqmax.ge.5)then
!         filagu=filagu+(b4/24)*(-24*gammq(1+power,x)+96*gammq(2+power,x)&
!           -72*gammq(3+power,x)+16*gammq(4+power,x)-gammq(5+power,x))
!       endif
! !
!       if(iqmax.ge.6)then
!         filagu=filagu+(b5/120)*(-120*gammq(1+power,x)+600*gammq(2+power,x)&
!           -600*gammq(3+power,x)+200*gammq(4+power,x)&
!           -25*gammq(5+power,x)+gammq(6+power,x))
!       endif
! !
!       return
!       end
!
!=====INTEGRAL OF LAGUERRE FUNCTION
!
      subroutine qalt(a,b,s)
!
!       s=(filagu(b)-filagu(a))*(1./(a-b))
!       s=(1./(a-b))*filagu(a,b)
!       s=xpower3(b)-xpower3(a)
!       s=xpowern(b)-xpowern(a)
!       s=(xpowern(b)-xpowern(a))/(b-a)
      s=xpowerm1(b)-xpowerm1(a)
!
      return
      end
!
!-ANALYTICAL INTERGRAL OF n(x)*x**3.dx
!
  function xpower3(x)
!
  common/moments/q0(1,1,1),q1(1,1,1),q2(1,1,1),q3(1,1,1)&
    ,q4(1,1,1),q5(1,1,1),q6(1,1,1),q7(1,1,1)&
    ,iqmax,power
  common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
!
  n=iqmax
  term1=0.;term2=0.;term3=0.;term4=0.;term5=0.;term6=0.;sum=0.
!
  if(n.ge.1)then
    a=1.;coeff=b0
    term1=a*(-6. -6.*x -3.*x**2 -x**3)
    term1=term1*coeff
  endif
!
  if(n.ge.2)then
    a=1.;b=-1.;coeff=b1
    term2=-a*(6. +6.*x +3.*x**2 +x**3) &
      -b*(24. +24.*x +12.*x**2 +4.*x**3 +x**4)
    term2=term2*coeff
  endif
!
  if(n.ge.3)then
    a=2.;b=-4.;c=1.;coeff=b2/2.
    term3=-6.*(a +4*b +20.*c) &
      -6.*(a +4.*b +20.*c)*x &
      -3.*(a +4.*b +20.*c)*x**2 &
      +(-a -4.*b -20.*c)*x**3 &
      +(-b -5.*c)*x**4 &
      -c*x**5
    term3=term3*coeff
  endif
!
  if(n.ge.4)then
    a=6.;b=-18.;c=9.;d=-1.;coeff=b3/6.
    term4=-6.*(a +4.*b +20.*c +120.*d) &
      -6.*(a +4.*b +20.*c +120.*d)*x &
      -3.*(a +4.*b +20.*c +120.*d)*x**2 &
      +(-a -4.*b -20.*c -120.*d)*x**3 &
      +(-b -5.*c -30.*d)*x**4 &
      +(-c -6.*d)*x**5 &
      -d*x**6
    term4=term4*coeff
  endif
!
  if(n.ge.5)then
    a=24.;b=-96.;c=72.;d=-16.;e=1.;coeff=b4/24.
    term5=-6.*(a +4.*b +20.*c +120.*d +840.*e) &
      -6.*(a +4.*b +20.*c +120.*d +840.*e)*x &
      -3.*(a +4.*b +20.*c +120.*d +840.*e)*x**2 &
      +(-a -4.*b -20.*c -120.*d -840.*e)*x**3 &
      +(-b -5.*c -30.*d -210.*e)*x**4 &
      +(-c -6.*d -42.*e)*x**5 &
      +(-d -7.*e)*x**6 &
      -e*x**7
    term5=term5*coeff
  endif
!
  if(n.ge.6)then
    a=120.;b=-600.;c=600.;d=-200.;e=25.;f=-1.;coeff=b5/120.
    term6=-6.*(a +4.*b +20.*c +120.*d +840.*e +6720.*f) &
      -6.*(a +4.*b +20.*c +120.*d +840.*e +6720.*f)*x &
      -3.*(a +4.*b +20.*c +120.*d +840.*e +6720.*f)*x**2 &
      +(-a -4.*b -20.*c -120.*d -840.*e -6720.*f)*x**3 &
      +(-b -5.*c -30.*d -210.*e -1680.*f)*x**4 &
      +(-c -6.*d -42.*e -336.*f)*x**5 &
      +(-d -7.*e -56.*f)*x**6 &
      +(-e -8.*f)*x**7 &
      -f*x**8
    term6=term6*coeff
  endif
!
  sum=exp(-x)*(term1+term2+term3+term4+term5+term6)
  xpower3=sum
!
  return
  end
!
!-ANALYTICAL INTERGRAL OF n(x)*x**power.dx
!
  function xpowern(x)
  common/moments/q0(1,1,1),q1(1,1,1),q2(1,1,1),q3(1,1,1)&
    ,q4(1,1,1),q5(1,1,1),q6(1,1,1),q7(1,1,1)&
    ,iqmax,power
  common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
  double precision g1,g2,g3,g4,g5,g6
  double precision a,b,c,d,e,f
!
  n=iqmax
  term1=0.;term2=0.;term3=0.;term4=0.;term5=0.;term6=0.
!
  g1=gammq(1.+power,x)
  g2=gammq(2.+power,x)
  g3=gammq(3.+power,x)
  g4=gammq(4.+power,x)
  g5=gammq(5.+power,x)
  g6=gammq(6.+power,x)
!
  if(n.ge.1)then
    a=1.;coeff=b0
    term1=-a*g1
    term1=term1*coeff
  endif
!
  if(n.ge.2)then
    a=1.;b=-1.;coeff=b1
    term2=-a*g1-b*g2
    term2=term2*coeff
  endif
!
  if(n.ge.3)then
    a=2.;b=-4.;c=1.;coeff=b2/2.
    term3=-a*g1-b*g2-c*g3
    term3=term3*coeff
  endif
!
  if(n.ge.4)then
    a=6.;b=-18.;c=9.;d=-1.;coeff=b3/6.
    term4=-a*g1-b*g2-c*g3-d*g4
    term4=term4*coeff
  endif
!
  if(n.ge.5)then
    a=24.;b=-96.;c=72.;d=-16.;e=1.;coeff=b4/24.
    term5=-a*g1-b*g2-c*g3-d*g4-e*g5
    term5=term5*coeff
  endif
!
  if(n.ge.6)then
    a=120.;b=-600.;c=600.;d=-200.;e=25.;f=-1.;coeff=b5/120.
    term6=-a*g1-b*g2-c*g3-d*g4-e*g5-f*g6
    term6=term6*coeff
  endif
!
  xpowern=exp(-x)*(term1+term2+term3+term4+term5+term6)
!
  return
  end
!
!=====INTEGRAL OF LAGUERRE FUNCTION
!
      function filagu(xlb,xub)
      common/moments/q0(1,1,1),q1(1,1,1),q2(1,1,1),q3(1,1,1)&
        ,q4(1,1,1),q5(1,1,1),q6(1,1,1),q7(1,1,1)&
        ,iqmax,power
      common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
      real a(6),b(6),c(6),d(6),e(6),f(6),coeff(6)
!
      sum=0.
      do n=1,iqmax
        a(n)=0;b(n)=0;c(n)=0;d(n)=0;e(n)=0;f(n)=0;coeff(n)=0
      enddo
!
      do n=1,iqmax
        if(n.eq.1)a(n)=1.;coeff(n)=b0
        if(n.eq.2)a(n)=1.;b(n)=1.;coeff(n)=b1
        if(n.eq.3)a(n)=2.;b(n)=-4.;c(n)=1.;coeff(n)=b2/2.
        if(n.eq.4)a(n)=6.;b(n)=-18.;c(n)=9.;d(n)=-1.;coeff(n)=b3/6.
        if(n.eq.5)a(n)=24.;b(n)=-96.;c(n)=72.;d(n)=-16.;e(n)=1.;coeff(n)=b4/24.
        if(n.eq.6)a(n)=120.;b(n)=-600.;c(n)=600.;d(n)=-200.;e(n)=25.;f(n)=-1.;coeff(n)=b5/120.
!
        sum=coeff(n)* &
          (a(n)*(gammq(1+power,xub)-gammq(1+power,xlb)) &
          +b(n)*(gammq(2+power,xub)-gammq(2+power,xlb)) &
          +c(n)*(gammq(3+power,xub)-gammq(3+power,xlb)) &
          +d(n)*(gammq(4+power,xub)-gammq(4+power,xlb)) &
          +e(n)*(gammq(5+power,xub)-gammq(5+power,xlb)) &
          +f(n)*(gammq(6+power,xub)-gammq(6+power,xlb)))+sum
      enddo
      filagu=sum
!
      return
      end
!
!=====INCOMPLETE GAMMA FUNCTION P(A,X) [Num. Res. in F77, Ch 6.2]
!
      function gammp(a,x)
!
!     GIVEN A AND X, RETURNS THE INCOMPLETE GAMMA FUNCTION P(A,X).
!     P(A,0.)=0., P(A,INFINITY)=1.
!
      real a,gammp,x
      real gammcf,gamser,gln
!
      if(a.lt.0. .or.x.le.0.)print*, 'bad arguments in gammp'
      if(x.lt.a+1)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.-gammcf
      endif
!
      if(gammp.lt.1.e-8)gammp=0.
!
      return
      end
!
!=====INCOMPLETE GAMMA FUNCTION Q(A,X) [Num. Res. in F77, Ch 6.2]
!
      function gammq(a,x)
!
!     GIVEN A AND X, RETURNS THE INCOMPLETE GAMMA FUNCTION Q(A,X).
!     Q(A,X)=1-P(A,X)
!     Q(A,0.)=1., Q(A,INFINITY)=0.
!
      real a,gammq,x
      real gammcf,gamser,gln
!
      if(a.lt.0. .or.x.le.0.)print*, 'bad arguments in gammp'
      if(x.lt.a+1)then
        call gser(gamser,a,x,gln)
        gammq=1.-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
!
      if(gammq.lt.1.e-30)gammq=0.
!
      return
      end
!
!=====SERIES REPRESENTATION
!     OF THE INCOMPLETE GAMMA FUNCTION [Num. Res. in F77, Ch 6.2]
!
      subroutine gser(gamser,a,x,gln)
!
!     RETURNS THE INCOMPLETE GAMMA FUNCTION P(A,X)
!
      integer itmax
      real a,gamser,gln,x,eps
!
      parameter(itmax=100,eps=3.e-7)
      integer n
      real ap,del,sum,gammln
!
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)print*, 'x < 0 in gser'
        gamser=0.
        return
      endif
!
      ap=a
      sum=1./a
      del=sum
      do n=1,itmax
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*eps)goto 100
      enddo
100   gamser=sum*exp(-x+a*log(x)-gln)
      if(gamser.lt.0.)gamser=0.
!
      return
      end
!
!=====CONTINUED FRACTION REPRESENTATION
!     OF THE INCOMPLETE GAMMA FUNCTION [Num. Res. in F77, Ch 6.2]
!
      subroutine gcf(gammcf,a,x,gln)
!
!     RETURNS THE INCOMPLETE GAMMA FUNCTION Q(A,X)
!
      integer itmax
      real a,gammcf,gln,x,eps,fpmin
!
      parameter(itmax=100,eps=3.e-7,fpmin=1.e-30)
      integer i
      real an,b,c,d,del,h,gammln
!
      gln=gammln(a)
      b=x+1.-a
      c=1./fpmin
      d=1./b
      h=d
      do i=1,itmax
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.fpmin)d=fpmin
        c=b+an/c
        if(abs(c).lt.fpmin)c=fpmin
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.eps)goto 100
      enddo
      print*, 'a is too large, itmax is too small in gcf'
100   gammcf=exp(-x+a*log(x)-gln)*h
!
      return
      end
!
!=====NATURAL LOGARITHM OF THE GAMMA FUNCTION [Num. Res. in F77, Ch 6.1]
!
      function gammln(xx)
!
!-----GIVEN XX, WHERE XX > 0, LN(G(XX)) IS RETURNED
!
      real gammln,xx
      integer j
      double precision ser,stp,tmp,x,y,cof(6)
      save cof,stp
!
      data cof,stp/&
        76.18009172947146d0&
        ,-86.50532032941677d0&
        ,24.01409824083091d0&
        ,-1.231739572450155d0&
        ,0.1208650973866179d-2&
        ,-0.5395239384953d-5&
        ,2.5066282746310005d0/
!
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.0000000001900d0
!
      do j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
      enddo
      gammln=tmp+log(stp*ser/x)
!
      return
      end
!
!-ANALYTICAL INTERGRAL OF n(x)*x**(-1).dx
!
  function xpowerm1(x)
  common/moments/q0(1,1,1),q1(1,1,1),q2(1,1,1),q3(1,1,1)&
    ,q4(1,1,1),q5(1,1,1),q6(1,1,1),q7(1,1,1)&
    ,iqmax,power
  common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
!
  n=iqmax
  term1=0.;term2=0.;term3=0.;term4=0.;term5=0.;term6=0.
  call eix(-x,expi)
  expx=exp(-x)
!
  if(n.ge.1)then
    a=1.;coeff=b0
    term1=0.
    term1=(term1+a*expi)*coeff
  endif
!
  if(n.ge.2)then
    a=1.;b=-1.;coeff=b1
    term2=-expx*(b)
    term2=(term2+a*expi)*coeff
  endif
!
  if(n.ge.3)then
    a=2.;b=-4.;c=1.;coeff=b2/2.
    term3=-expx*(b+c+c*x)
    term3=(term3+a*expi)*coeff
  endif
!
  if(n.ge.4)then
    a=6.;b=-18.;c=9.;d=-1.;coeff=b3/6.
    term4=-expx*(b+c*(1+x)+d*(2+2*x+x**2))
    term4=(term4+a*expi)*coeff
  endif
!
  if(n.ge.5)then
    a=24.;b=-96.;c=72.;d=-16.;e=1.;coeff=b4/24.
    term5=-expx*(b+c+2*d+6*e+c*x+2*d*x+6*e*x+d*x**2+3*e*x**2+e*x**3)
    term5=(term5+a*expi)*coeff
  endif
!
  if(n.ge.6)then
    a=120.;b=-600.;c=600.;d=-200.;e=25.;f=-1.;coeff=b5/120.
    term6=expx*(-b-c-2*d-6*e-24*f+(-c-2*d-6*e-24*f)*x+(-d-3*e-12*f)*x**2+(-e-4*f)*x**3-f*x**4)
    term6=(term6+a*expi)*coeff
  endif
!
  xpowerm1=(term1+term2+term3+term4+term5+term6)
!
  return
  end
!
!=EXPONENTIAL INTEGRAL, Ei(x)
!
  subroutine eix(x,ei)
!
!============================================
! input :  x  --- argument of ei(x)
! output:  ei --- ei(x) ( x > 0 )
!============================================
!
  implicit double precision (a-h,o-z)
  if(x.eq.0.0)then
    ei=-1.0d+30
  elseif(x.le.40.0)then
    ei=1.0d0
    r=1.0d0
    do 15 k=1,100
      r=r*k*x/(k+1.0d0)**2
      ei=ei+r
      if(dabs(r/ei).le.1.0d-15)goto 20
15  continue
20  ga=0.5772156649015328d0
    ei=ga+dlog(x)+x*ei
  else
    ei=1.0d0
    r=1.0d0
    do 25 k=1,20
      r=r*k/x
25  ei=ei+r
    ei=dexp(x)/x*ei
  endif
!
  return
  end