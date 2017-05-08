!
!-ANALYTICAL INTERGRAL OF n(x)*x**(-1).dx
!
  function xpowerm1(x)
!
  term1=0.;term2=0.;term3=0.;term4=0.;term5=0.;term6=0.
  expi=eix(-x)
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
    term3=-expx(b+c+c*x)
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
  xpowerm1=exp(-x)*(term1+term2+term3+term4+term5+term6)
!
  return
  end