!
!-ANALYTICAL INTERGRAL OF n(x)*x**3.dx
!
  function xpower3(x)
!
  term1=0.;term2=0.;term3=0.;term4=0.;term5=0.;term6=0.
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
  xpower3=exp(-x)*(term1+term2+term3+term4+term5+term6)
!
  return
  end
!
!-ANALYTICAL INTERGRAL OF n(x)*x**power.dx
!
  function xpowern(x,power)
!
  term1=0.;term2=0.;term3=0.;term4=0.;term5=0.;term6=0.
!
  if(n.ge.1)then
    g1=gammq(1.+power,x)
    a=1.;coeff=b0
    term1=-a*g1
    term1=term1*coeff
  endif
!
  if(n.ge.2)then
    g2=gammq(2.+power,x)
    a=1.;b=-1.;coeff=b1
    term2=-a*g1-b*g2
    term2=term2*coeff
  endif
!
  if(n.ge.3)then
    g3=gammq(3.+power,x)
    a=2.;b=-4.;c=1.;coeff=b2/2.
    term3=-a*g1-b*g2-c*g3
    term3=term3*coeff
  endif
!
  if(n.ge.4)then
    g4=gammq(4.+power,x)
    a=6.;b=-18.;c=9.;d=-1.;coeff=b3/6.
    term4=-a*g1-b*g2-c*g3-d*g4
    term4=term4*coeff
  endif
!
  if(n.ge.5)then
    g5=gammq(5.+power,x)
    a=24.;b=-96.;c=72.;d=-16.;e=1.;coeff=b4/24.
    term5=-a*g1-b*g2-c*g3-d*g4-e*g5
    term5=term5*coeff
  endif
!
  if(n.ge.6)then
    g6=gammq(6.+power,x)
    a=120.;b=-600.;c=600.;d=-200.;e=25.;f=-1.;coeff=b5/120.
    term6=-a*g1-b*g2-c*g3-d*g4-e*g5-f*g6
    term6=term6*coeff
  endif
!
  xpowern=exp(-x)*(term1+term2+term3+term4+term5+term6)
!
  return
  end
