!
!-ANALYTICAL INTERGRAL OF n(x)*x**4.dx
!
  function x4(x)
  common/moments/q0(1,1,1),q1(1,1,1),q2(1,1,1),q3(1,1,1)&
    ,q4(1,1,1),q5(1,1,1),q6(1,1,1),q7(1,1,1)&
    ,iqmax,rexpm
  common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
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