!
!=====TRUNCATED MOMENT
!
      function qtrunc(n,xlb,xub,i,j,k)
      include "common.inc"
      
      term1=pk(i,j,k)+n/2.
      term2=q0(i,j,k)
      term3=(r32(i,j,k)/(pk(i,j,k)+2.))**(n/2.)
      term4=exp(gammln(term1)-gammln(pk(i,j,k)))
      term5=gammp(term1,xub)-gammp(term1,xlb)
      if(term5.lt.0.)term5=0.
      
      qtrunc=term2*term3*term4*term5
      
      return
      end
!
!-----FUNCTION GAMMP(A,X)
!
      function gammp(a,x)
!
!.....uses GCF, GSER
!     returns the incomplete gamma function P(a,x)
!
      real a,gammp,x
      real gamser,gln
      if(a.lt.0.)a=0.0001
      if(x.lt.0.)x=0.0001
      if(x.gt.50.)x=50.
      
      call gser(gamser,a,x,gln)
      gammp=gamser

      return
      end
!
!-----SUBROUTINE GSER(GAMSER,A,X,GLN)
!
      subroutine gser(gamser,a,x,gln)
!
!.....uses GAMMLN
!     returns the incomplete gamma function P(a,x)
!
      integer itmax
      real a,gamser,gln,x,eps
      
      parameter(itmax=100,eps=3.e-7)
      integer n
      real ap,del,sum1,sum2,gammln
!
      anew=a
      sum1=0.
      if(int(a).gt.1)then
        mmax=int(a)-1
        do m=1,mmax
          sum1=sum1+exp((a-m)*log(x)-gammln(a-m+1))
        enddo
        anew=a-mmax
      endif
      sum1=-exp(-x)*sum1
!
      gln=gammln(anew)
      if(x.lt.0.)x=0.
      ap=anew
      sum2=1./anew
      del=sum2
      do n=1,itmax
        ap=ap+1
        del=del*x/ap
        sum2=sum2+del
        if(abs(del).lt.abs(sum2)*eps)goto 1
      enddo
1     gamser=sum1+sum2*exp(-x+anew*log(x)-gln)
      if(gamser.lt.0.)gamser=0.
      
      return
      end
!
!-----GAMMLN
!
      function gammln(xx)
      real gammln,xx
      integer j
      double precision ser,stp,tmp,x,y,cof(6)
      save cof,stp
      data cof,stp/
     &  76.18009172947146d0
     &  ,-86.50532032941677d0
     &  ,24.01409824083091d0
     &  ,-1.231739572450155d0
     &  ,0.1208650973866179d-2
     &  ,-0.5395239384953d-5
     &  ,2.5066282746310005d0/
      
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.0000000001900d0
      do j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
      enddo
      gammln=tmp+log(stp*ser/x)
      
      return
      end
