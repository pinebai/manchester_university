!
!=====GAS PRESSURE
!
      subroutine calcp
      include "common.inc"
      real vel(6,ia,ja,ka)
!
      ivel=0;iph=1;iq=3
!
      call zerocoef
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff0(i,j,k)=dng0(i,j,k)*(1.-cnt*q30(i,j,k))
            ff(i,j,k)=dng(i,j,k)*(1.-cnt*q3(i,j,k))
            do nf=1,6
              vel(nf,i,j,k)=0.
            enddo
            pc(i,j,k)=0.
          enddo
        enddo
      enddo
!
      call zerocoef
!
!-----SOLVE COEFFICIENTS
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            aw(i,j,k)=area(1,i,j,k)**2*(1.-cnt*f(1,i,j,k,q3))**2&
              *f(1,i,j,k,dng)/f(1,i,j,k,apu)
            ae(i,j,k)=area(2,i,j,k)**2*(1.-cnt*f(2,i,j,k,q3))**2&
              *f(2,i,j,k,dng)/f(2,i,j,k,apu)
            as(i,j,k)=area(3,i,j,k)**2*(1.-cnt*f(3,i,j,k,q3))**2&
              *f(3,i,j,k,dng)/f(3,i,j,k,apv)
            an(i,j,k)=area(4,i,j,k)**2*(1.-cnt*f(4,i,j,k,q3))**2&
              *f(4,i,j,k,dng)/f(4,i,j,k,apv)
            ab(i,j,k)=area(5,i,j,k)**2*(1.-cnt*f(5,i,j,k,q3))**2&
              *f(5,i,j,k,dng)/f(5,i,j,k,apw)
            at(i,j,k)=area(6,i,j,k)**2*(1.-cnt*f(6,i,j,k,q3))**2&
              *f(6,i,j,k,dng)/f(6,i,j,k,apw)
          enddo
        enddo
      enddo
!
!-----ZERO BOUNDARY COEFFICIENTS
!
      call zerobcoef
!
!-----MAIN COEFFICIENT
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            ap(i,j,k)=aw(i,j,k)+ae(i,j,k)+as(i,j,k)&
              +an(i,j,k)+ab(i,j,k)+at(i,j,k)-sp(i,j,k)
          enddo
        enddo
      enddo
!
!-----CELL FACE VELOCITIES
!
!
!-----WEST
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            if(i.ne.imn)then
              term1=f(1,i,j,k,ug)
              coef1=(1.-cnt*f(1,i,j,k,q3))*f(1,i,j,k,vol)/f(1,i,j,k,apu)
              grad1=(p(i,j,k)-p(i-1,j,k))/(xc(i)-xc(i-1))
              coef2=(1.-cnt*q3(i,j,k))*vol(i,j,k)/apu(i,j,k)
              grad2=dpx(i,j,k)
              coef3=(1.-cnt*q3(i-1,j,k))*vol(i-1,j,k)/apu(i-1,j,k)
              grad3=dpx(i-1,j,k)
              vel(1,i,j,k)=term1-coef1*grad1+(coef2*grad2+coef3*grad3)/2
            endif
            if(i.eq.imn)vel(1,i,j,k)=ug(i-1,j,k)
          enddo
        enddo
      enddo
!
!-----EAST
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            if(i.ne.imx)then
              term1=f(2,i,j,k,ug)
              coef1=(1.-cnt*f(2,i,j,k,q3))*f(2,i,j,k,vol)/f(2,i,j,k,apu)
              grad1=(p(i+1,j,k)-p(i,j,k))/(xc(i+1)-xc(i))
              coef2=(1.-cnt*q3(i,j,k))*vol(i,j,k)/apu(i,j,k)
              grad2=dpx(i,j,k)
              coef3=(1.-cnt*q3(i+1,j,k))*vol(i+1,j,k)/apu(i+1,j,k)
              grad3=dpx(i+1,j,k)
              vel(2,i,j,k)=term1-coef1*grad1+(coef2*grad2+coef3*grad3)/2
            endif
            if(i.eq.imx)vel(2,i,j,k)=ug(i+1,j,k)
          enddo
        enddo
      enddo
!
!-----SOUTH
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            if(j.ne.jmn)then
              term1=f(3,i,j,k,vg)
              coef1=(1.-cnt*f(3,i,j,k,q3))*f(3,i,j,k,vol)/f(3,i,j,k,apv)
              grad1=(p(i,j,k)-p(i,j-1,k))/(yc(j)-yc(j-1))
              coef2=(1.-cnt*q3(i,j,k))*vol(i,j,k)/apv(i,j,k)
              grad2=dpy(i,j,k)
              coef3=(1.-cnt*q3(i,j-1,k))*vol(i,j-1,k)/apv(i,j-1,k)
              grad3=dpy(i,j-1,k)
              vel(3,i,j,k)=term1-coef1*grad1+(coef2*grad2+coef3*grad3)/2
            endif
            if(j.eq.jmn)vel(3,i,j,k)=vg(i,j-1,k)
          enddo
        enddo
      enddo
!
!-----NORTH
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            if(j.ne.jmx)then
              term1=f(4,i,j,k,vg)
              coef1=(1.-cnt*f(4,i,j,k,q3))*f(4,i,j,k,vol)/f(4,i,j,k,apv)
              grad1=(p(i,j+1,k)-p(i,j,k))/(yc(j+1)-yc(j))
              coef2=(1.-cnt*q3(i,j,k))*vol(i,j,k)/apv(i,j,k)
              grad2=dpy(i,j,k)
              coef3=(1.-cnt*q3(i,j+1,k))*vol(i,j+1,k)/apv(i,j+1,k)
              grad3=dpy(i,j+1,k)
              vel(4,i,j,k)=term1-coef1*grad1+(coef2*grad2+coef3*grad3)/2
            endif
            if(j.eq.jmx)vel(4,i,j,k)=vg(i,j+1,k)
          enddo
        enddo
      enddo
!
!-----BOTTOM
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            if(k.ne.kmn)then
              term1=f(5,i,j,k,wg)
              coef1=(1.-cnt*f(5,i,j,k,q3))*f(5,i,j,k,vol)/f(5,i,j,k,apw)
              grad1=(p(i,j,k)-p(i,j,k-1))/(zc(k)-zc(k-1))
              coef2=(1.-cnt*q3(i,j,k))*vol(i,j,k)/apw(i,j,k)
              grad2=dpz(i,j,k)
              coef3=(1.-cnt*q3(i,j,k-1))*vol(i,j,k-1)/apw(i,j,k-1)
              grad3=dpz(i,j,k-1)
              vel(5,i,j,k)=term1-coef1*grad1+(coef2*grad2+coef3*grad3)/2
            endif
            if(k.eq.kmn)vel(5,i,j,k)=wg(i,j,k-1)
          enddo
        enddo
      enddo
!
!-----TOP
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            if(k.ne.kmx)then
              term1=f(6,i,j,k,wg)
              coef1=(1.-cnt*f(6,i,j,k,q3))*f(6,i,j,k,vol)/f(6,i,j,k,apw)
              grad1=(p(i,j,k+1)-p(i,j,k))/(zc(k+1)-zc(k))
              coef2=(1.-cnt*q3(i,j,k))*vol(i,j,k)/apw(i,j,k)
              grad2=dpz(i,j,k)
              coef3=(1.-cnt*q3(i,j,k+1))*vol(i,j,k+1)/apw(i,j,k+1)
              grad3=dpz(i,j,k+1)
              vel(6,i,j,k)=term1-coef1*grad1+(coef2*grad2+coef3*grad3)/2
            endif
            if(k.eq.kmx)vel(6,i,j,k)=wg(i,j,k+1)
          enddo
        enddo
      enddo
!
!-----CONTINUITY SOURCE TERM
!
!       res(iequ)=0
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            fw(i,j,k)=area(1,i,j,k)*f(1,i,j,k,ff)*vel(1,i,j,k)
            fe(i,j,k)=area(2,i,j,k)*f(2,i,j,k,ff)*vel(2,i,j,k)
            fs(i,j,k)=area(3,i,j,k)*f(3,i,j,k,ff)*vel(3,i,j,k)
            fn(i,j,k)=area(4,i,j,k)*f(4,i,j,k,ff)*vel(4,i,j,k)
            fb(i,j,k)=area(5,i,j,k)*f(5,i,j,k,ff)*vel(5,i,j,k)
            ft(i,j,k)=area(6,i,j,k)*f(6,i,j,k,ff)*vel(6,i,j,k)
            rdt=(ff(i,j,k)-ff0(i,j,k))*vol(i,j,k)/delt
            su(i,j,k)=fw(i,j,k)-fe(i,j,k)+fs(i,j,k)&
              -fn(i,j,k)+fb(i,j,k)-ft(i,j,k)-rdt
!             res(iequ)=res(iequ)+abs(su(i,j,k))+tiny
          enddo
        enddo
      enddo
!
!-----SOLVE PRESSURE CORRECTION EQUATION
!
      call solve(pc,1.0)
!
!-----UPDATE PRESSURE
!
      iequ=2
      ppref=pc(imx,jmx,kmx)
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            p(i,j,k)=p(i,j,k)+urf(iequ)*(pc(i,j,k)-ppref)
          enddo
        enddo
      enddo
!
!-----EXTRAPOLATE BOUNDARY VALUES
!
      do nf=1,6
        call nbc(nf,pc)
        call nbc(nf,p)
      enddo
!
!-----UPDATE VELOCITIES
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            dpx(i,j,k)=(f(2,i,j,k,p)-f(1,i,j,k,p))/dx(i)
            dppx=(f(2,i,j,k,pc)-f(1,i,j,k,pc))/dx(i)
            ug(i,j,k)=ug(i,j,k)-(1.-cnt*q3(i,j,k))*vol(i,j,k)&
              *dppx/apu(i,j,k)
!
            dpy(i,j,k)=(f(4,i,j,k,p)-f(3,i,j,k,p))/dy(j)
            dppy=(f(4,i,j,k,pc)-f(3,i,j,k,pc))/dy(j)
            vg(i,j,k)=vg(i,j,k)-(1.-cnt*q3(i,j,k))*vol(i,j,k)&
              *dppy/apv(i,j,k)
!
            dpz(i,j,k)=(f(6,i,j,k,p)-f(5,i,j,k,p))/dz(k)
            dppz=(f(6,i,j,k,pc)-f(5,i,j,k,pc))/dz(k)
            wg(i,j,k)=wg(i,j,k)-(1.-cnt*q3(i,j,k))*vol(i,j,k)&
              *dppz/apw(i,j,k)
          enddo
        enddo
      enddo
!
      return
      end