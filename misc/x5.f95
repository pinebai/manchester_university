!
!-ANALYTICAL INTERGRAL OF n(x)*x**5.dx
!
  function x5(x)
  common/moments/q0(1,1,1),q1(1,1,1),q2(1,1,1),q3(1,1,1)&
    ,q4(1,1,1),q5(1,1,1),q6(1,1,1),q7(1,1,1)&
    ,iqmax,rexpm
  common/lcoef/b0,b1,b2,b3,b4,b5,b6,b7
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