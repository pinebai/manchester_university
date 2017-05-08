!
!=====SPATIAL DISCRETISATION SCHEME: HYBRID
!
      subroutine hybrid
      include "common.inc"
      
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            aw(i,j,k)=max(fw(i,j,k),(dw(i,j,k)+fw(i,j,k)/2),0)
            ae(i,j,k)=max(-fe(i,j,k),(de(i,j,k)-fe(i,j,k)/2),0)
            as(i,j,k)=max(fs(i,j,k),(ds(i,j,k)+fs(i,j,k)/2),0)
            an(i,j,k)=max(-fn(i,j,k),(dn(i,j,k)-fn(i,j,k)/2),0)
            ab(i,j,k)=max(fb(i,j,k),(db(i,j,k)+fb(i,j,k)/2),0)
            at(i,j,k)=max(-ft(i,j,k),(dt(i,j,k)-ft(i,j,k)/2),0)
          enddo
        enddo
      enddo

      return
      end
!
!=====SPATIAL DISCRETISATION SCHEME: QUICK [HAYASE]
!
      subroutine quick(phi)
      include "common.inc"
      real phi(ia,ja,ka)
      
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            if(fw(i,j,k).lt.tiny)qw=0
            if(fw(i,j,k).gt.-tiny)qw=1
            if(fe(i,j,k).lt.tiny)qe=0
            if(fe(i,j,k).gt.-tiny)qe=1
            if(fs(i,j,k).lt.tiny)qs=0
            if(fs(i,j,k).gt.-tiny)qs=1
            if(fn(i,j,k).lt.tiny)qn=0
            if(fn(i,j,k).gt.-tiny)qn=1
            if(fb(i,j,k).lt.tiny)qb=0
            if(fb(i,j,k).gt.-tiny)qb=1
            if(ft(i,j,k).lt.tiny)qt=0
            if(ft(i,j,k).gt.-tiny)qt=1
!
!-----MAIN COEFFICIENTS
!
            aw(i,j,k)=dw(i,j,k)+qw*fw(i,j,k)
            ae(i,j,k)=de(i,j,k)-(1-qe)*fe(i,j,k)
            as(i,j,k)=ds(i,j,k)+qs*fs(i,j,k)
            an(i,j,k)=dn(i,j,k)-(1-qn)*fn(i,j,k)
            ab(i,j,k)=db(i,j,k)+qb*fb(i,j,k)
            at(i,j,k)=dt(i,j,k)-(1-qt)*ft(i,j,k)
!
!-----phi(ww) at the west boundary is linearly interpolated
!     from phi(w) & phi(p), etc.
!
            term1a=1/8.*(3*phi(i,j,k)-2*phi(i-1,j,k)-phi(i-2,j,k))&
                   *qw*fw(i,j,k)
            if(i.eq.imn)term1a=1/8.*(4*phi(i,j,k)-4*phi(i-1,j,k))&
                               *qw*fw(i,j,k)
            term1b=1/8.*(phi(i-1,j,k)+2*phi(i,j,k)-3*phi(i+1,j,k))&
                   *qe*fe(i,j,k)
            term2a=1/8.*(3*phi(i-1,j,k)-2*phi(i,j,k)-phi(i+1,j,k))&
                   *(1-qw)*fw(i,j,k)
            term2b=1/8.*(2*phi(i+1,j,k)+phi(i+2,j,k)-3*phi(i,j,k))&
                   *(1-qe)*fe(i,j,k)
            if(i.eq.imx)term2b=1/8.*(4*phi(i+1,j,k)-4*phi(i,j,k))&
                               *(1-qe)*fe(i,j,k)
!
            term3a=1/8.*(3*phi(i,j,k)-2*phi(i,j-1,k)-phi(i,j-2,k))&
                   *qs*fs(i,j,k)
            if(j.eq.jmn)term3a=1/8.*(4*phi(i,j,k)-4*phi(i,j-1,k))&
                               *qs*fs(i,j,k)
            term3b=1/8.*(phi(i,j-1,k)+2*phi(i,j,k)-3*phi(i,j+1,k))&
                   *qn*fn(i,j,k)
            term4a=1/8.*(3*phi(i,j-1,k)-2*phi(i,j,k)-phi(i,j+1,k))&
                   *(1-qs)*fs(i,j,k)
            term4b=1/8.*(2*phi(i,j+1,k)+phi(i,j+2,k)-3*phi(i,j,k))&
                   *(1-qn)*fn(i,j,k)
            if(j.eq.jmx)term4b=1/8.*(4*phi(i,j+1,k)-4*phi(i,j,k))&
                               *(1-qn)*fn(i,j,k)
!
            term5a=1/8.*(3*phi(i,j,k)-2*phi(i,j,k-1)-phi(i,j,k-2))&
                   *qb*fb(i,j,k)
            if(k.eq.kmn)term5a=1/8.*(4*phi(i,j,k)-4*phi(i,j,k-1))&
                               *qb*fb(i,j,k)
            term5b=1/8.*(phi(i,j,k-1)+2*phi(i,j,k)-3*phi(i,j,k+1))&
                   *qt*ft(i,j,k)
            term6a=1/8.*(3*phi(i,j,k-1)-2*phi(i,j,k)-phi(i,j,k+1))&
                   *(1-qb)*fb(i,j,k)
            term6b=1/8.*(2*phi(i,j,k+1)+phi(i,j,k+2)-3*phi(i,j,k))&
                   *(1-qt)*ft(i,j,k)
            if(k.eq.kmx)term6b=1/8.*(4*phi(i,j,k+1)-4*phi(i,j,k))&
                               *(1-qt)*ft(i,j,k)
!
!-----SOURCE TERM
!
            su(i,j,k)=su(i,j,k)&
                      +term1a+term1b+term2a+term2b+term3a+term3b&
                      +term4a+term4b+term5a+term5b+term6a+term6b
          enddo
        enddo
      enddo

      return
      end
!
!=====TIME DISCRETISATION SCHEME
!
      subroutine temporal(phi0)
      include "common.inc"
      real phi0(ia,ja,ka)
      
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
!
!-----MASS FLUX IMBALACE TERM
!
            delf=fe(i,j,k)-fw(i,j,k)&
                 +fn(i,j,k)-fs(i,j,k)&
                 +ft(i,j,k)-fb(i,j,k)
!
!-----TIME STEP COEFFICIENT
!
            rdt0=ff0(i,j,k)*vol(i,j,k)/delt
            rdt=ff(i,j,k)*vol(i,j,k)/delt
!
!-----MAIN COEFFICIENT
!
            if(ivel.eq.0.and.iph.eq.2)then 
              ap(i,j,k)=rdt&
                        +aw(i,j,k)+ae(i,j,k)&
                        +as(i,j,k)+an(i,j,k)&
                        +ab(i,j,k)+at(i,j,k)&
                        -sp(i,j,k)+delf
            else
              ap(i,j,k)=rdt0&
                        +aw(i,j,k)+ae(i,j,k)&
                        +as(i,j,k)+an(i,j,k)&
                        +ab(i,j,k)+at(i,j,k)&
                        -sp(i,j,k)
            endif
!
!-----SOURCE TERM
!
            su(i,j,k)=su(i,j,k)+rdt0*phi0(i,j,k)
          enddo
        enddo
      enddo

      return
      end

      