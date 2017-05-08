!
!=====INTEGRATE FUNCTION
!
      program area
      external fx
!
      a=0.05
      b=2.
      do b=0.05,2.,0.05
        call qtrap(fx,a,b,s)
        print*,'s =',s
      enddo
!
      end
!
!=====FUNCTION
!
      function fx(x)
!
      fx=exp(-x)/x
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