!
!-ANALYTICAL INTERGRAL OF n(x)*x**1.dx
!
  function xpower0(x)
!
  term1=0.;term2=0.;term3=0.;term4=0.;term5=0.;term6=0.
!
  if(n.ge.1)then
    a=1.;coeff=b0
    term1=-a*(1. +x)
    term1=term1*coeff
  endif
!
  if(n.ge.2)then
    a=1.;b=-1.;coeff=b1
    term2=-a*(1. +x) &
      -b*x
    term2=term2*coeff
  endif
!
  if(n.ge.3)then
    a=2.;b=-4.;c=1.;coeff=b2/2.
    term3=-a -b-2.*c
      +(
    term3=term3*coeff
  endif
!
  if(n.ge.4)then
    a=6.;b=-18.;c=9.;d=-1.;coeff=b3/6.
    term4=(-a -b -2.*c -6.*d) &
      +(-b -2.*c -6.*d)*x &
      +(-c -3.*d)*x**2 &
      +(-d)*x**3
    term4=term4*coeff
  endif
!
  if(n.ge.5)then
    a=24.;b=-96.;c=72.;d=-16.;e=1.;coeff=b4/24.
    term5=(-a -b -2.*c -6.*d -24.*e) &
      +(-b -2.*c -6.*d -24.*e)*x &
      +(-c -3.*d -12.*e)*x**2 &
      +(-d -4.*e)*x**3
      +(-e)*x**4
    term5=term5*coeff
  endif
!
  if(n.ge.6)then
    a=120.;b=-600.;c=600.;d=-200.;e=25.;f=-1.;coeff=b5/120.
    term6=(-a -b -2.*c -6.*d -24.*e -120.*f) &
      +(-b -2.*c -6.*d -24.*e -120.*f)*x &
      +(-c -3.*d -12.*e -60.*f)*x**2 &
      +(-d -4.*e -20.*f)*x**3
      +(-e -5.*f)*x**4
      +(-f)*x**5
    term6=term6*coeff
  endif
!
  xpower0=exp(-x)*(term1+term2+term3+term4+term5+term6)
!
  return
  end