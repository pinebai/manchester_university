!
!=====MAIN COEFFICIENTS
!
      subroutine mcf(phi0)
      include "common.inc"

      real phi0(ia,ja)
      
      do i=imn,imx
        do j=jmn,jmx
!
!-----TIME STEP COEFFICIENT
!
          ap0(i,j)=den0(i,j)*qq0(i,j)*vol(i,j)/delt
!
!-----MAIN COEFFICIENT
!
          if(iap.eq.1)then
            ap(i,j)=aw(i,j)+ae(i,j)+as(i,j)+an(i,j)
     &              -sp(i,j)+ap0(i,j)
          elseif(iap.eq.2)then
            ap(i,j)=aw(i,j)+ae(i,j)+as(i,j)+an(i,j)
     &              +fe(i,j)-fw(i,j)+fn(i,j)-fs(i,j)
     &              -sp(i,j)+ap0(i,j)
          endif
!
!-----MOMENT SOURCE TERM
!
          if(ivel.eq.0.and.iph.eq.2)then
            su(i,j)=su(i,j)+ap0(i,j)*phi0(i,j)
!
!-----GAS PHASE MOMENTUM SOURCE TERM
!
          elseif(ivel.eq.1.and.iph.eq.1)then
            su(i,j)=su(i,j)+ap0(i,j)*phi0(i,j)+drg(i,j)
     &              -(p(i+1,j)-p(i-1,j))*0.5*(areae(i,j)+areaw(i,j))
     &              *(1-cnt*q3(i,j))
            awu(i,j)=aw(i,j)
            aeu(i,j)=ae(i,j)
            asu(i,j)=as(i,j)
            anu(i,j)=an(i,j)
            suu(i,j)=su(i,j)
            apu(i,j)=ap(i,j)
          elseif(ivel.eq.2.and.iph.eq.1)then
            su(i,j)=su(i,j)+ap0(i,j)*phi0(i,j)+drg(i,j)
     &              -(p(i,j+1)-p(i,j-1))*0.5*(arean(i,j)+areas(i,j))
     &              *(1-cnt*q3(i,j))
            awv(i,j)=aw(i,j)
            aev(i,j)=ae(i,j)
            asv(i,j)=as(i,j)
            anv(i,j)=an(i,j)
            suv(i,j)=su(i,j)
            apv(i,j)=ap(i,j)
!
!-----LIQUID PHASE MOMENTUM SOURCE TERM
!
          elseif(ivel.ge.1.and.iph.eq.2.and.iq.ge.1)then
            su(i,j)=su(i,j)+ap0(i,j)*phi0(i,j)-drg(i,j)
          endif
        enddo
      enddo

      return
      end