!
!=====MOMENTS OF THE SIZE DISTRIBUTION
!
      subroutine moment(qq,power,rlb,rub,i,j,k)
      include "common.inc"
      common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
      external flagu
!
      qq=0.
!
!-----SIZE DISTRIBUTION NORMALISATION RADIUS
!
      if(iqmax.eq.4)rnorm=(q3(i,j,k)/q0(i,j,k))**(1/3.)
      if(iqmax.eq.5)rnorm=(q4(i,j,k)/q0(i,j,k))**(1/4.)
      if(iqmax.eq.6)rnorm=(q5(i,j,k)/q0(i,j,k))**(1/5.)
      if(iqmax.eq.7)rnorm=(q6(i,j,k)/q0(i,j,k))**(1/6.)
      if(iqmax.eq.8)rnorm=(q7(i,j,k)/q0(i,j,k))**(1/7.)
      if(rnorm.lt.1.e-6.or.rnorm.gt.0.001)return
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
!       print*,'A'
!
!-----LIMITS OF DISTRIBUTION FUNCTION
!
      root1=0; root2=0; delr=0.5e-6
!
      if(rlb.lt.1.e-6)then
        rrlb=1.e-6
      else
        root1=rlb
        rrlb=rlb
      endif
!
      if(rub.lt.1.e-6)then
        rrub=50.e-6
      else
        root2=rub
        rrub=rub
      endif
!
      if(abs(root1).lt.tiny)then
      pfunc=3.
      do r=rrlb,rrub,delr
        terml=flagu((r-delr)/rnorm,pfunc)
        termc=flagu(r/rnorm,pfunc)
        termu=flagu((r+delr)/rnorm,pfunc)
        if(termc.gt.tiny.and.terml.lt.termu)then
          root1=r
          goto 100
        endif
      enddo
!
      do r=rrlb,rrub,delr
        terml=flagu((r-delr)/rnorm,pfunc)
        termc=flagu(r/rnorm,pfunc)
        termu=flagu((r+delr)/rnorm,pfunc)
        if(terml.gt.tiny.and.termc.gt.tiny.and.termu.gt.tiny)then
          root1=r
          goto 100
        endif
      enddo
100   continue
      endif
!
      if(abs(root2).lt.tiny)then
      do r=rrub,rrlb,-delr
        terml=flagu((r-delr)/rnorm,pfunc)
        termc=flagu(r/rnorm,pfunc)
        termu=flagu((r+delr)/rnorm,pfunc)
        if(termc.gt.tiny.and.terml.gt.termu.and.r.gt.root1)then
          root2=r
          goto 200
        endif
      enddo
200   continue
      endif
!       print*,'B'
!
      if(root1.gt.tiny.and.root2.gt.root1)then
        rlbn=root1/rnorm
        rubn=root2/rnorm
!
        icount=0;icountmx=10
300     continue
        icount=icount+1
!
        pfunc=3.
        call qalt(flagu,pfunc,rlbn,rubn,qan)
        qnorm=q0(i,j,k)/(1./rnorm)**(pfunc)
        qa=qan*qnorm
        ratio=qa/q3(i,j,k)
!
        if(ratio.lt.0.98 .and.icount.lt.icountmx)then
          if(ratio.gt.0.8)rubn=rubn+0.02
          if(ratio.lt.0.8)rubn=rubn+0.2
          goto 300
        elseif(ratio.gt.1.02 .and.icount.lt.icountmx)then
          if(ratio.lt.1.2)rubn=rubn-0.02
          if(ratio.gt.1.2)rubn=rubn-0.2
          goto 300
        endif
!
        if(rlbn.gt.tiny.and.rubn.gt.rlbn)then
!           print*,'C'
!
!-----INTEGRATION
!
          call qalt(flagu,power,rlbn,rubn,qqn)
!
!-----UN-NORMALISE MOMENT
!
          qnorm=q0(i,j,k)/(1./rnorm)**(power)
          qq=qqn*qnorm
          rlb=root1
          rub=root2
          if(qq.lt.tiny)qq=0.
!           print*,'D'
        endif
      endif
!
      return
      end
!
!=====LAGUERRE SIZE DISTRIBUTION FUNCTION
!
      function flagu(x,power)
      common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
      common/qmax/iqmax
!
      flagu=0.
      if(iqmax.ge.4)then
        flagu=exp(-x)*x**(power)&
          *(b0*plagu(0,x)+b1*plagu(1,x)&
          +b2*plagu(2,x)+b3*plagu(3,x))
      endif
      if(iqmax.ge.5)flagu=flagu+exp(-x)*x**(power)*b4*plagu(4,x)
      if(iqmax.ge.6)flagu=flagu+exp(-x)*x**(power)*b5*plagu(5,x)
      if(iqmax.ge.7)flagu=flagu+exp(-x)*x**(power)*b6*plagu(6,x)
      if(iqmax.ge.8)flagu=flagu+exp(-x)*x**(power)*b7*plagu(7,x)
!       if(flagu.lt.0.)flagu=0.
!
      return
      end
!
!=====GAMMA SIZE DISTRIBUTION FUNCTION
!
!       function fgam(x,power)
!       common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
!       include "common.inc"
! !
!       fgam=0.
!       r=x*rnorm
!       rr21=q2(i,j,k)/q1(i,j,k)
!       rr32=q3(i,j,k)/q2(i,j,k)
!       ratio=rr21/rr32
!       if(ratio.gt.0.5 .and.ratio.lt.0.95)then
!         pk=(1-2*ratio)/(1-ratio)
!         fgam=(exp(pk*log(pk+2)-gammln(pk))*(r**(pk-1)/rr32**pk) &
!           exp(-(pk+2)*(r/rr32)))/rnorm
!       endif
! !
!       return
!       end
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
!=INTEGRAL OF LAGUERRE FUNCTION
!
  subroutine qalt(func,power,a,b,s)
  external func
!
  if(power.lt.-1.999 .and.power.gt.-2.001)then
    s=xm2(b)-xm2(a)
  elseif(power.lt.-0.999 .and.power.gt.-1.001)then
    s=xm1(b)-xm1(a)
  elseif(power.gt.-0.001 .and.power.lt.0.001)then
    s=x0(b)-x0(a)
  elseif(power.gt.0.999 .and.power.lt.1.001)then
    s=x1(b)-x1(a)
  elseif(power.gt.1.999 .and.power.lt.2.001)then
    s=x2(b)-x2(a)
  elseif(power.gt.2.999 .and.power.lt.3.001)then
    s=x3(b)-x3(a)
  elseif(power.gt.3.999 .and.power.lt.4.001)then
    s=x4(b)-x4(a)
  elseif(power.gt.4.999 .and.power.lt.5.001)then
    s=x5(b)-x5(a)
  elseif(power.gt.-0.5)then
    s=(xn(b,power)-xn(a,power))/(b-a)
  else
    call qtrap(func,power,a,b,s)
  endif
!
  return
  end
!
!=====INTEGRATION REFINEMENT [Num. Res. in F77, Ch 4.2]
!
      subroutine qtrap(func,power,a,b,s)
!
      integer jmax
      real func,power,a,b,s,eps
      external func
      parameter (eps=1.e-6,jmax=5)
!
      integer j
      real olds
!
      olds=-1.e30
      do j=1,jmax
        call trapzd(func,power,a,b,s,j)
        if(j.gt.5)then
          if(abs(s-olds).lt.eps*abs(olds).or.&
            (s.eq.0. .and.olds.eq.0))return
        endif
        olds=s
      enddo
!
      return
      end
!
!=====NUMERICAL INTEGRATION
!     USING EXTENDED TRAPEZOIDAL RULE [Num. Res. in F77, Ch 4.2]
!
      subroutine trapzd(func,power,a,b,s,n)
!
      integer n
      real func,power,a,b,s
      external func
!
      integer it,j
      real del,sum,tnm,x
!
      if(n.eq.1)then
        s=0.5*(b-a)*(func(a,power)+func(b,power))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do j=1,it
          sum=sum+func(x,power)
          x=x+del
        enddo
        s=0.5*(s+(b-a)*sum/tnm)
      endif
!
      return
      end
!
!-ANALYTICAL INTERGRAL OF n(x)*x**(-2).dx
!
  function xm2(x)
  common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
  common/qmax/iqmax
!
  n=iqmax
  expi=(-0.0268/5)*x**5+(0.2343/4)*x**4+(-0.7457/3)*x**3+(0.8977/2)*x**2+(0.0108)*x
!
  if(n.ge.1)then
    a=1.;coeff=b0
    term=a*(-exp(-x)/x-expi)
    term1=term*coeff
  endif
!
  if(n.ge.2)then
    a=1.;b=-1.;coeff=b1
    term=-a*exp(-x)/x+(-a+b)*expi
    term2=term*coeff
  endif
!
  if(n.ge.3)then
    a=2.;b=-4.;c=1.;coeff=b2/2.
    term=-exp(-x)*(a+c*x)/x+(-a+b)*expi
    term3=term*coeff
  endif
!
  if(n.ge.4)then
    a=6.;b=-18.;c=9.;d=-1.;coeff=b3/6.
    term=exp(-x)*(-c-d-a/x-d*x)-a*expi+b*expi
    term4=term*coeff
  endif
!
  if(n.ge.5)then
    a=24.;b=-96.;c=72.;d=-16.;e=1.;coeff=b4/24.
    term=exp(-x)*(-a-x*(c+d*(1+x)+e*(2+2*x+x**2))-(a-b)*exp(-x)*x*expi)/x
    term5=term*coeff
  endif
!
  if(n.ge.6)then
    a=120.;b=-600.;c=600.;d=-200.;e=25.;f=-1.;coeff=b5/120.
    term=exp(-x)*(-c-d-2*e-6*f-a/x+(-d-2*e-6*f)*x+(-e-3*f)*x**2-f*x**3)-a*expi+b*expi
    term6=term*coeff
  endif
!
  xm2=(term1+term2+term3+term4+term5+term6)
!
  return
  end
!
!-ANALYTICAL INTERGRAL OF n(x)*x**(-1).dx
!
  function xm1(x)
  common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
  common/qmax/iqmax
!
  n=iqmax
  expi=(-0.0268/5)*x**5+(0.2343/4)*x**4+(-0.7457/3)*x**3+(0.8977/2)*x**2+(0.0108)*x
!
  if(n.ge.1)then
    a=1.;coeff=b0
    coeff0=0.
    term=coeff0*x**0
    term=exp(-x)*term
    term1=(term+a*expi)*coeff
  endif
!
  if(n.ge.2)then
    a=1.;b=-1.;coeff=b1
    coeff0=(-b)
    term=coeff0*x**0
    term=exp(-x)*term
    term2=(term+a*expi)*coeff
  endif
!
  if(n.ge.3)then
    a=2.;b=-4.;c=1.;coeff=b2/2.
    coeff0=(-b-c)
    coeff1=(-c)
    term=coeff0*x**0+coeff1*x**1
    term=exp(-x)*term
    term3=(term+a*expi)*coeff
  endif
!
  if(n.ge.4)then
    a=6.;b=-18.;c=9.;d=-1.;coeff=b3/6.
    coeff0=(-b-c-2*d)
    coeff1=(-c-2*d)
    coeff2=(-d)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2
    term=exp(-x)*term
    term4=(term+a*expi)*coeff
  endif
!
  if(n.ge.5)then
    a=24.;b=-96.;c=72.;d=-16.;e=1.;coeff=b4/24.
    coeff0=(-b-c-2*d-6*e)
    coeff1=(-c-2*d-6*e)
    coeff2=(-d-3*e)
    coeff3=(-e)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3
    term=exp(-x)*term
    term5=(term+a*expi)*coeff
  endif
!
  if(n.ge.6)then
    a=120.;b=-600.;c=600.;d=-200.;e=25.;f=-1.;coeff=b5/120.
    coeff0=(-b-c-2*d-6*e-24*f)
    coeff1=(-c-2*d-6*e-24*f)
    coeff2=(-d-3*e-12*f)
    coeff3=(-e-4*f)
    coeff4=(-f)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4
    term=exp(-x)*term
    term6=(term+a*expi)*coeff
  endif
!
  xm1=(term1+term2+term3+term4+term5+term6)
!
  return
  end
!
!-ANALYTICAL INTERGRAL OF n(x)*x**0.dx
!
  function x0(x)
  common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
  common/qmax/iqmax
!
  n=iqmax
!
  if(n.ge.1)then
    a=1.;coeff=b0
    coeff0=(-a)
    term=coeff0*x**0
    term1=coeff*term
  endif
!
  if(n.ge.2)then
    a=1.;b=-1.;coeff=b1
    coeff0=(-a-b)
    coeff1=(-b)
    term=coeff0*x**0+coeff1*x**1
    term2=coeff*term
  endif
!
  if(n.ge.3)then
    a=2.;b=-4.;c=1.;coeff=b2/2.
    coeff0=(-a-b-2*c)
    coeff1=(-b-2*c)
    coeff2=(-c)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2
    term3=coeff*term
  endif
!
  if(n.ge.4)then
    a=6.;b=-18.;c=9.;d=-1.;coeff=b3/6.
    coeff0=(-a-b-2*c-6*d)
    coeff1=(-b-2*c-6*d)
    coeff2=(-c-3*d)
    coeff3=(-d)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3
    term4=coeff*term
  endif
!
  if(n.ge.5)then
    a=24.;b=-96.;c=72.;d=-16.;e=1.;coeff=b4/24.
    coeff0=(-a-b-2*c-6*d-24*e)
    coeff1=(-b-2*c-6*d-24*e)
    coeff2=(-c-3*d-12*e)
    coeff3=(-d-4*e)
    coeff4=(-e)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4
    term5=coeff*term
  endif
!
  if(n.ge.6)then
    a=120.;b=-600.;c=600.;d=-200.;e=25.;f=-1.;coeff=b5/120.
    coeff0=(-a-b-2*c-6*d-24*e-120*f)
    coeff1=(-b-2*c-6*d-24*e-120*f)
    coeff2=(-c-3*d-12*e-60*f)
    coeff3=(-d-4*e-20*f)
    coeff4=(-e-5*f)
    coeff5=(-f)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5
    term6=coeff*term
  endif
!
  x0=exp(-x)*(term1+term2+term3+term4+term5+term6)
!
  return
  end
!
!-ANALYTICAL INTERGRAL OF n(x)*x**1.dx
!
  function x1(x)
  common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
  common/qmax/iqmax
!
  n=iqmax
!
  if(n.ge.1)then
    a=1.;coeff=b0
    coeff0=(-a)
    coeff1=(-a)
    term=coeff0*x**0+coeff1*x**1
    term1=coeff*term
  endif
!
  if(n.ge.2)then
    a=1.;b=-1.;coeff=b1
    coeff0=(-a-2*b)
    coeff1=(-a-2*b)
    coeff2=(-b)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2
    term2=coeff*term
  endif
!
  if(n.ge.3)then
    a=2.;b=-4.;c=1.;coeff=b2/2.
    coeff0=(-a-2*b-6*c)
    coeff1=(-a-2*b-6*c)
    coeff2=(-b-3*c)
    coeff3=(-c)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3
    term3=coeff*term
  endif
!
  if(n.ge.4)then
    a=6.;b=-18.;c=9.;d=-1.;coeff=b3/6.
    coeff0=(-a-2*b-6*c-24*d)
    coeff1=(-a-2*b-6*c-24*d)
    coeff2=(-b-3*c-12*d)
    coeff3=(-c-4*d)
    coeff4=(-d)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4
    term4=coeff*term
  endif
!
  if(n.ge.5)then
    a=24.;b=-96.;c=72.;d=-16.;e=1.;coeff=b4/24.
    coeff0=(-a-2*b-6*c-24*d-120*e)
    coeff1=(-a-2*b-6*c-24*d-120*e)
    coeff2=(-b-3*c-12*d-60*e)
    coeff3=(-c-4*d-20*e)
    coeff4=(-d-5*e)
    coeff5=(-e)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5
    term5=coeff*term
  endif
!
  if(n.ge.6)then
    a=120.;b=-600.;c=600.;d=-200.;e=25.;f=-1.;coeff=b5/120.
    coeff0=(-a-2*b-6*c-24*d-120*e-720*f)
    coeff1=(-a-2*b-6*c-24*d-120*e-720*f)
    coeff2=(-b-3*c-12*d-60*e-360*f)
    coeff3=(-c-4*d-20*e-120*f)
    coeff4=(-d-5*e-30*f)
    coeff5=(-e-6*f)
    coeff6=(-f)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5+coeff6*x**6
    term6=coeff*term
  endif
!
  x1=exp(-x)*(term1+term2+term3+term4+term5+term6)
!
  return
  end
!
!-ANALYTICAL INTERGRAL OF n(x)*x**2.dx
!
  function x2(x)
  common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
  common/qmax/iqmax
!
  n=iqmax
!
  if(n.ge.1)then
    a=1.;coeff=b0
    coeff0=(-2*a)
    coeff1=(-2*a)
    coeff2=(-a)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2
    term1=coeff*term
  endif
!
  if(n.ge.2)then
    a=1.;b=-1.;coeff=b1
    coeff0=(-2*a-6*b)
    coeff1=(-2*a-6*b)
    coeff2=(-a-3*b)
    coeff3=(-b)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3
    term2=coeff*term
  endif
!
  if(n.ge.3)then
    a=2.;b=-4.;c=1.;coeff=b2/2.
    coeff0=(-2*a-6*b-24*c)
    coeff1=(-2*a-6*b-24*c)
    coeff2=(-a-3*b-12*c)
    coeff3=(-b-4*c)
    coeff4=(-c)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4
    term3=coeff*term
  endif
!
  if(n.ge.4)then
    a=6.;b=-18.;c=9.;d=-1.;coeff=b3/6.
    coeff0=(-2*a-6*b-24*c-120*d)
    coeff1=(-2*a-6*b-24*c-120*d)
    coeff2=(-a-3*b-12*c-60*d)
    coeff3=(-b-4*c-20*d)
    coeff4=(-c-5*d)
    coeff5=(-d)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5
    term4=coeff*term
  endif
!
  if(n.ge.5)then
    a=24.;b=-96.;c=72.;d=-16.;e=1.;coeff=b4/24.
    coeff0=(-2*a-6*b-24*c-120*d-720*e)
    coeff1=(-2*a-6*b-24*c-120*d-720*e)
    coeff2=(-a-3*b-12*c-60*d-360*e)
    coeff3=(-b-4*c-20*d-120*e)
    coeff4=(-c-5*d-30*e)
    coeff5=(-d-6*e)
    coeff6=(-e)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5+coeff6*x**6
    term5=coeff*term
  endif
!
  if(n.ge.6)then
    a=120.;b=-600.;c=600.;d=-200.;e=25.;f=-1.;coeff=b5/120.
    coeff0=(-2*a-6*b-24*c-120*d-720*e-5040*f)
    coeff1=(-2*a-6*b-24*c-120*d-720*e-5040*f)
    coeff2=(-a-3*b-12*c-60*d-360*e-2520*f)
    coeff3=(-b-4*c-20*d-120*e-840*f)
    coeff4=(-c-5*d-30*e-210*f)
    coeff5=(-d-6*e-42*f)
    coeff6=(-e-7*f)
    coeff7=(-f)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5+coeff6*x**6+coeff7*x**7
    term6=coeff*term
  endif
!
  x2=exp(-x)*(term1+term2+term3+term4+term5+term6)
!
  return
  end
!
!-ANALYTICAL INTERGRAL OF n(x)*x**3.dx
!
  function x3(x)
  common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
  common/qmax/iqmax
!
  n=iqmax
!
  if(n.ge.1)then
    a=1.;coeff=b0
    coeff0=(-6*a)
    coeff1=(-6*a)
    coeff2=(-3*a)
    coeff3=(-a)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3
    term1=coeff*term
  endif
!
  if(n.ge.2)then
    a=1.;b=-1.;coeff=b1
    coeff0=(-6*a-24*b)
    coeff1=(-6*a-24*b)
    coeff2=(-3*a-12*b)
    coeff3=(-a-4*b)
    coeff4=(-b)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4
    term2=coeff*term
  endif
!
  if(n.ge.3)then
    a=2.;b=-4.;c=1.;coeff=b2/2.
    coeff0=(-6*a-24*b-120*c)
    coeff1=(-6*a-24*b-120*c)
    coeff2=(-3*a-12*b-60*c)
    coeff3=(-a-4*b-20*c)
    coeff4=(-b-5*c)
    coeff5=(-c)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5
    term3=coeff*term
  endif
!
  if(n.ge.4)then
    a=6.;b=-18.;c=9.;d=-1.;coeff=b3/6.
    coeff0=(-6*a-24*b-120*c-720*d)
    coeff1=(-6*a-24*b-120*c-720*d)
    coeff2=(-3*a-12*b-60*c-360*d)
    coeff3=(-a-4*b-20*c-120*d)
    coeff4=(-b-5*c-30*d)
    coeff5=(-c-6*d)
    coeff6=(-d)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5+coeff6*x**6
    term4=coeff*term
  endif
!
  if(n.ge.5)then
    a=24.;b=-96.;c=72.;d=-16.;e=1.;coeff=b4/24.
    coeff0=(-6*a-24*b-120*c-720*d-5040*e)
    coeff1=(-6*a-24*b-120*c-720*d-5040*e)
    coeff2=(-3*a-12*b-60*c-360*d-2520*e)
    coeff3=(-a-4*b-20*c-120*d-840*e)
    coeff4=(-b-5*c-30*d-210*e)
    coeff5=(-c-6*d-42*e)
    coeff6=(-d-7*e)
    coeff7=(-e)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5+coeff6*x**6+coeff7*x**7
    term5=coeff*term
  endif
!
  if(n.ge.6)then
    a=120.;b=-600.;c=600.;d=-200.;e=25.;f=-1.;coeff=b5/120.
    coeff0=(-6*a-24*b-120*c-720*d-5040*e-40320*f)
    coeff1=(-6*a-24*b-120*c-720*d-5040*e-40320*f)
    coeff2=(-3*a-12*b-60*c-360*d-2520*e-20160*f)
    coeff3=(-a-4*b-20*c-120*d-840*e-6720*f)
    coeff4=(-b-5*c-30*d-210*e-1680*f)
    coeff5=(-c-6*d-42*e-336*f)
    coeff6=(-d-7*e-56*f)
    coeff7=(-e-8*f)
    coeff8=(-f)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5+coeff6*x**6+coeff7*x**7 &
      +coeff8*x**8
    term6=coeff*term
  endif
!
  x3=exp(-x)*(term1+term2+term3+term4+term5+term6)
!
  return
  end
!
!-ANALYTICAL INTERGRAL OF n(x)*x**4.dx
!
  function x4(x)
  common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
  common/qmax/iqmax
!
  n=iqmax
!
  if(n.ge.1)then
    a=1.;coeff=b0
    coeff0=(-24*a)
    coeff1=(-24*a)
    coeff2=(-12*a)
    coeff3=(-4*a)
    coeff4=(-a)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4
    term1=coeff*term
  endif
!
  if(n.ge.2)then
    a=1.;b=-1.;coeff=b1
    coeff0=(-24*a-120*b)
    coeff1=(-24*a-120*b)
    coeff2=(-12*a-60*b)
    coeff3=(-4*a-20*b)
    coeff4=(-a-5*b)
    coeff5=(-b)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5
    term2=coeff*term
  endif
!
  if(n.ge.3)then
    a=2.;b=-4.;c=1.;coeff=b2/2.
    coeff0=(-24*a-120*b-720*c)
    coeff1=(-24*a-120*b-720*c)
    coeff2=(-12*a-60*b-360*c)
    coeff3=(-4*a-20*b-120*c)
    coeff4=(-a-5*b-30*c)
    coeff5=(-b-6*c)
    coeff6=(-c)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5+coeff6*x**6
    term3=coeff*term
  endif
!
  if(n.ge.4)then
    a=6.;b=-18.;c=9.;d=-1.;coeff=b3/6.
    coeff0=(-24*a-120*b-720*c-5040*d)
    coeff1=(-24*a-120*b-720*c-5040*d)
    coeff2=(-12*a-60*b-360*c-2520*d)
    coeff3=(-4*a-20*b-120*c-840*d)
    coeff4=(-a-5*b-30*c-210*d)
    coeff5=(-b-6*c-42*d)
    coeff6=(-c-7*d)
    coeff7=(-d)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5+coeff6*x**6+coeff7*x**7
    term4=coeff*term
  endif
!
  if(n.ge.5)then
    a=24.;b=-96.;c=72.;d=-16.;e=1.;coeff=b4/24.
    coeff0=(-24*a-120*b-720*c-5040*d-40320*e)
    coeff1=(-24*a-120*b-720*c-5040*d-40320*e)
    coeff2=(-12*a-60*b-360*c-2520*d-20160*e)
    coeff3=(-4*a-20*b-120*c-840*d-6720*e)
    coeff4=(-a-5*b-30*c-210*d-1680*e)
    coeff5=(-b-6*c-42*d-336*e)
    coeff6=(-c-7*d-56*e)
    coeff7=(-d-8*e)
    coeff8=(-e)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5+coeff6*x**6+coeff7*x**7 &
      +coeff8*x**8
    term5=coeff*term
  endif
!
  if(n.ge.6)then
    a=120.;b=-600.;c=600.;d=-200.;e=25.;f=-1.;coeff=b5/120.
    coeff0=(-24*a-120*b-720*c-5040*d-40320*e-362880*f)
    coeff1=(-24*a-120*b-720*c-5040*d-40320*e-362880*f)
    coeff2=(-12*a-60*b-360*c-2520*d-20160*e-181440*f)
    coeff3=(-4*a-20*b-120*c-840*d-6720*e-60480*f)
    coeff4=(-a-5*b-30*c-210*d-1680*e-15120*f)
    coeff5=(-b-6*c-42*d-336*e-3024*f)
    coeff6=(-c-7*d-56*e-504*f)
    coeff7=(-d-8*e-72*f)
    coeff8=(-e-9*f)
    coeff9=(-f)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5+coeff6*x**6+coeff7*x**7 &
      +coeff8*x**8+coeff9*x**9
    term6=coeff*term
  endif
!
  x4=exp(-x)*(term1+term2+term3+term4+term5+term6)
!
  return
  end
!
!-ANALYTICAL INTERGRAL OF n(x)*x**5.dx
!
  function x5(x)
  common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
  common/qmax/iqmax
!
  n=iqmax
!
  if(n.ge.1)then
    a=1.;coeff=b0
    coeff0=(-120*a)
    coeff1=(-120*a)
    coeff2=(-60*a)
    coeff3=(-20*a)
    coeff4=(-5*a)
    coeff5=(-a)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5
    term1=coeff*term
  endif
!
  if(n.ge.2)then
    a=1.;b=-1.;coeff=b1
    coeff0=(-120*a-720*b)
    coeff1=(-120*a-720*b)
    coeff2=(-60*a-360*b)
    coeff3=(-20*a-120*b)
    coeff4=(-5*a-30*b)
    coeff5=(-a-6*b)
    coeff6=(-b)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5+coeff6*x**6
    term2=coeff*term
  endif
!
  if(n.ge.3)then
    a=2.;b=-4.;c=1.;coeff=b2/2.
    coeff0=(-120*a-720*b-5040*c)
    coeff1=(-120*a-720*b-5040*c)
    coeff2=(-60*a-360*b-2520*c)
    coeff3=(-20*a-120*b-840*c)
    coeff4=(-5*a-30*b-210*c)
    coeff5=(-a-6*b-42*c)
    coeff6=(-b-7*c)
    coeff7=(-c)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5+coeff6*x**6+coeff7*x**7
    term3=coeff*term
  endif
!
  if(n.ge.4)then
    a=6.;b=-18.;c=9.;d=-1.;coeff=b3/6.
    coeff0=(-120*a-720*b-5040*c-40320*d)
    coeff1=(-120*a-720*b-5040*c-40320*d)
    coeff2=(-60*a-360*b-2520*c-20160*d)
    coeff3=(-20*a-120*b-840*c-6720*d)
    coeff4=(-5*a-30*b-210*c-1680*d)
    coeff5=(-a-6*b-42*c-336*d)
    coeff6=(-b-7*c-56*d)
    coeff7=(-c-8*d)
    coeff8=(-d)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5+coeff6*x**6+coeff7*x**7 &
      +coeff8*x**8
    term4=coeff*term
  endif
!
  if(n.ge.5)then
    a=24.;b=-96.;c=72.;d=-16.;e=1.;coeff=b4/24.
    coeff0=(-120*a-720*b-5040*c-40320*d-362880*e)
    coeff1=(-120*a-720*b-5040*c-40320*d-362880*e)
    coeff2=(-60*a-360*b-2520*c-20160*d-181440*e)
    coeff3=(-20*a-120*b-840*c-6720*d-60480*e)
    coeff4=(-5*a-30*b-210*c-1680*d-15120*e)
    coeff5=(-a-6*b-42*c-336*d-3024*e)
    coeff6=(-b-7*c-56*d-504*e)
    coeff7=(-c-8*d-72*e)
    coeff8=(-d-9*e)
    coeff9=(-e)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5+coeff6*x**6+coeff7*x**7 &
      +coeff8*x**8+coeff9*x**9
    term5=coeff*term
  endif
!
  if(n.ge.6)then
    a=120.;b=-600.;c=600.;d=-200.;e=25.;f=-1.;coeff=b5/120.
    coeff0=(-120*a-720*b-5040*c-40320*d-362880*e-3628800*f)
    coeff1=(-120*a-720*b-5040*c-40320*d-362880*e-3628800*f)
    coeff2=(-60*a-360*b-2520*c-20160*d-181440*e-1814400*f)
    coeff3=(-20*a-120*b-840*c-6720*d-60480*e-604800*f)
    coeff4=(-5*a-30*b-210*c-1680*d-15120*e-151200*f)
    coeff5=(-a-6*b-42*c-336*d-3024*e-30240*f)
    coeff6=(-b-7*c-56*d-504*e-5040*f)
    coeff7=(-c-8*d-72*e-720*f)
    coeff8=(-d-9*e-90*f)
    coeff9=(-e-10*f)
    coeff10=(-f)
    term=coeff0*x**0+coeff1*x**1+coeff2*x**2+coeff3*x**3 &
      +coeff4*x**4+coeff5*x**5+coeff6*x**6+coeff7*x**7 &
      +coeff8*x**8+coeff9*x**9+coeff10*x**10
    term6=coeff*term
  endif
!
  x5=exp(-x)*(term1+term2+term3+term4+term5+term6)
!
  return
  end
!
!-ANALYTICAL INTERGRAL OF f(x)*x**n.dx
!
  function xn(x,power)
  common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
  common/qmax/iqmax
!
  n=iqmax
!
  if(n.ge.1)then
    a=1.;coeff=b0
    g1=gammq(1.+power,x)
    term=-a*g1
    term1=term*coeff
  endif
!
  if(n.ge.2)then
    a=1.;b=-1.;coeff=b1
    g2=gammq(2.+power,x)
    term=-a*g1-b*g2
    term2=term*coeff
  endif
!
  if(n.ge.3)then
    a=2.;b=-4.;c=1.;coeff=b2/2.
    g3=gammq(3.+power,x)
    term=-a*g1-b*g2-c*g3
    term3=term*coeff
  endif
!
  if(n.ge.4)then
    a=6.;b=-18.;c=9.;d=-1.;coeff=b3/6.
    g4=gammq(4.+power,x)
    term=-a*g1-b*g2-c*g3-d*g4
    term4=term*coeff
  endif
!
  if(n.ge.5)then
    a=24.;b=-96.;c=72.;d=-16.;e=1.;coeff=b4/24.
    g5=gammq(5.+power,x)
    term=-a*g1-b*g2-c*g3-d*g4-e*g5
    term5=term*coeff
  endif
!
  if(n.ge.6)then
    a=120.;b=-600.;c=600.;d=-200.;e=25.;f=-1.;coeff=b5/120.
    g6=gammq(6.+power,x)
    term=-a*g1-b*g2-c*g3-d*g4-e*g5-f*g6
    term6=term*coeff
  endif
!
  xn=exp(-x)*(term1+term2+term3+term4+term5+term6)
!
  return
  end
!
!=================
! GAMMA FUNCTIONS
!=================
!
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