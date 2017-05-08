!
!=====GAMMA FUNCTION
!
      function gmf(pk)
      
      pi=4.*atan(1.)
      gmf=((pk/exp(1.)*(pk*sinh(1/pk
     &    +1./(810.*pk**6)))**0.5)**pk)
     &    *(2.*pi/pk)**0.5

      return
      end