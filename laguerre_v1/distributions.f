!
!=====MOMENTS OF THE SIZE DISTRIBUTION
!
      function dmom(rlb,rub,i,j,k)
      include "common.inc"
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
      if(rlb.lt.1.e-7)rrlb=1.5e-6
      if(rub.lt.1.e-7)rrub=50.e-6
      delr=1.e-7
!
      root1=0;root2=0
      do r=rrlb,rrub,delr
        terml=flagu((r-delr)/rnorm)
        termc=flagu(r/rnorm)
        termu=flagu((r+delr)/rnorm)
        if(termc.gt.tiny.and.terml.lt.termu)then
          root1=r
          goto 100
        endif
      enddo
100   continue
      do r=rrub,rrlb,-delr
        terml=flagu((r-delr)/rnorm)
        termc=flagu(r/rnorm)
        termu=flagu((r+delr)/rnorm)
        if(termc.gt.tiny.and.terml.gt.termu)then
          root2=r
          goto 200
        endif
      enddo
200   continue
      if(root2.lt.root1)then
        root1=0;root2=0
      endif
      rlbn=root1/rnorm
      rubn=root2/rnorm
!
!-----NUMERICAL INTEGRATION
!
C       if(1+rexpm.lt.0.5)then
      call qtrap(flagu,rlbn,rubn,qnum)
!
!-----ANALYTICAL INTEGRATION
!
C       else
C         call qalt(rlbn,rubn,qnum)
C       endif
!
!-----UN-NORMALISE MOMENT
!
      qnorm=q0(i,j,k)/(1./rnorm)**(rexpm)
      dmom=qnum*qnorm
!
      return
      end
!
!=====SIZE DISTRIBUTION FUNCTION
!
      function flagu(xx)
      include "common.inc"
      common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
!
      flagu=0.
      if(iqmax.eq.4)then
        flagu=exp(-xx)*xx**(rexpm)
     &    *(b0*plagu(0,xx)
     &    +b1*plagu(1,xx)
     &    +b2*plagu(2,xx)
     &    +b3*plagu(3,xx))
      elseif(iqmax.eq.5)then
        flagu=exp(-xx)*xx**(rexpm)
     &    *(b0*plagu(0,xx)
     &    +b1*plagu(1,xx)
     &    +b2*plagu(2,xx)
     &    +b3*plagu(3,xx)
     &    +b4*plagu(4,xx))
      elseif(iqmax.eq.6)then
        flagu=exp(-xx)*xx**(rexpm)
     &    *(b0*plagu(0,xx)
     &    +b1*plagu(1,xx)
     &    +b2*plagu(2,xx)
     &    +b3*plagu(3,xx)
     &    +b4*plagu(4,xx)
     &    +b5*plagu(5,xx))
      elseif(iqmax.eq.7)then
        flagu=exp(-xx)*xx**(rexpm)
     &    *(b0*plagu(0,xx)
     &    +b1*plagu(1,xx)
     &    +b2*plagu(2,xx)
     &    +b3*plagu(3,xx)
     &    +b4*plagu(4,xx)
     &    +b5*plagu(5,xx)
     &    +b6*plagu(6,xx))
      elseif(iqmax.eq.8)then
        flagu=exp(-xx)*xx**(rexpm)
     &    *(b0*plagu(0,xx)
     &    +b1*plagu(1,xx)
     &    +b2*plagu(2,xx)
     &    +b3*plagu(3,xx)
     &    +b4*plagu(4,xx)
     &    +b5*plagu(5,xx)
     &    +b6*plagu(6,xx)
     &    +b7*plagu(7,xx))
      endif
!
      return
      end
!
!=====INTEGRAL FORM OF LAGUERRE FUNCTION
!
      function filagu(xx)
      include "common.inc"
      common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
!
      filagu=0.
     &  +b0*(-gammq(1+rexpm,xx))
     &  +b1*(-gammq(1+rexpm,xx)+gammq(2+rexpm,xx))
     &  +(b2/2)*(-2*gammq(1+rexpm,xx)+4*gammq(2+rexpm,xx)
     &  -gammq(3+rexpm,xx))
     &  +(b3/6)*(-6*gammq(1+rexpm,xx)+18*gammq(2+rexpm,xx)
     &  -9*gammq(3+rexpm,xx)+gammq(4+rexpm,xx))
     &  +(b4/24)*(-24*gammq(1+rexpm,xx)+96*gammq(2+rexpm,xx)
     &  -72*gammq(3+rexpm,xx)+16*gammq(4+rexpm,xx)-gammq(5+rexpm,xx))
     &  +(b5/120)*(-120*gammq(1+rexpm,xx)+600*gammq(2+rexpm,xx)
     &  -600*gammq(3+rexpm,xx)+200*gammq(4+rexpm,xx)
     &  -25*gammq(5+rexpm,xx)+gammq(6+rexpm,xx))
!
      return
      end
!
!=====INTEGRAL OF LAGUERRE FUNCTION
!
      subroutine qalt(a,b,s)
!
      s=filagu(b)-filagu(a)
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
        plagu=(1/2.)
     &  *(2.-4.*var+var**2)
      elseif(iqq.eq.3)then
        plagu=(1/6.)
     &  *(6.-18.*var+9.*var**2
     &  -var**3)
      elseif(iqq.eq.4)then
        plagu=(1/24.)
     &  *(24.-96.*var+72.*var**2
     &  -16.*var**3+var**4)
      elseif(iqq.eq.5)then
        plagu=(1/120.)
     &  *(120.-600.*var+600.*var**2
     &  -200.*var**3+25.*var**4-var**5)
      elseif(iqq.eq.6)then
        plagu=(1/720.)
     &  *(720.-4320.*var+5400.*var**2
     &  -2400.*var**3+450.*var**4-36.*var**5
     &  +var**6)
      elseif(iqq.eq.7)then
        plagu=(1/5040.)
     &  *(5040.-35280.*var+52920.*var**2
     &  -29400.*var**3+7350.*var**4-882.*var**5
     &  +49.*var**6-var**7)
      endif
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
      parameter (eps=1.e-6,jmax=15)
!
      integer j
      real olds
!
      olds=-1.e30
      do j=1,jmax
        call trapzd(func,a,b,s,j)
        if(j.gt.5)then
          if(abs(s-olds).lt.eps*abs(olds).or.
     &      (s.eq.0. .and.olds.eq.0))return
        endif
        olds=s
      enddo
C       print*, 'too many steps in qtrap'
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
      if(gammq.lt.1.e-8)gammq=0.
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
      data cof,stp/
     &  76.18009172947146d0
     &  ,-86.50532032941677d0
     &  ,24.01409824083091d0
     &  ,-1.231739572450155d0
     &  ,0.1208650973866179d-2
     &  ,-0.5395239384953d-5
     &  ,2.5066282746310005d0/
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
!=====LAGRANGE'S POLYNOMIAL INTERPOLATION [Num. Res. in F77, Ch 3.1]
!
      subroutine polint(xa,ya,n,x,y,dy)
!
!     GIVEN ARRAYS XA AND YA, EACH OF LENGTH N, AND GIVEN A VALUE X,
!     THIS ROUTINE RETURNS A VALUE Y, AND AN ERROR ESTIMATE DY.
!
!     e.g.
!       x(1)=0.;y(1)=0.
!       x(2)=2.;y(2)=2.
!       xp=1.
!       call polint(x(1),y(1),2,xp,yp,dy)
!       print*,xp,yp  [-> yp=1.0]
!
      integer n,nmax
      real dy,x,y,xa(n),ya(n)
      parameter (nmax=10)
!
      integer i,m,ns
      real den,dif,dift,ho,hp,w,c(nmax),d(nmax)
!
      ns=1
      dif=abs(x-xa(1))
!
      do i=1,n
        dift=abs(x-xa(i))
        if(dift.lt.dif)then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
      enddo
!
      y=ya(ns)
      ns=ns-1
!
      do m=1,n-1
        do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)print*, 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
        enddo
        if(2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
      enddo
!
      return
      end