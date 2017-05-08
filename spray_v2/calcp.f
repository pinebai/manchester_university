!
!=====GAS PRESSURE
!
      subroutine calcp
      include "common.inc"
      real vel(6,ia,ja,ka)
      
      do i=imn-1,imx+1
        do j=jmn-1,jmx+1
          do k=kmn-1,kmx+1
            ff00(i,j,k)=dng00(i,j,k)*(1.-cnt*q300(i,j,k))
            ff0(i,j,k)=dng0(i,j,k)*(1.-cnt*q30(i,j,k))
            ff(i,j,k)=dng(i,j,k)*(1.-cnt*q3(i,j,k))
            do nf=1,6
              vel(nf,i,j,k)=0.
            enddo
            pc(i,j,k)=0.
          enddo
        enddo
      enddo
      ivel=0;iph=1;iq=3;iequ=16
      call zerocoef
!
!-----SOLVE COEFFICIENTS
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            aw(i,j,k)=area(1,i,j,k)**2*(1.-cnt*f(1,i,j,k,q3))**2
     &                *f(1,i,j,k,dng)/f(1,i,j,k,apu)
            ae(i,j,k)=area(2,i,j,k)**2*(1.-cnt*f(2,i,j,k,q3))**2
     &                *f(2,i,j,k,dng)/f(2,i,j,k,apu)
            as(i,j,k)=area(3,i,j,k)**2*(1.-cnt*f(3,i,j,k,q3))**2
     &                *f(3,i,j,k,dng)/f(3,i,j,k,apv)
            an(i,j,k)=area(4,i,j,k)**2*(1.-cnt*f(4,i,j,k,q3))**2
     &                *f(4,i,j,k,dng)/f(4,i,j,k,apv)
            ab(i,j,k)=area(5,i,j,k)**2*(1.-cnt*f(5,i,j,k,q3))**2
     &                *f(5,i,j,k,dng)/f(5,i,j,k,apw)
            at(i,j,k)=area(6,i,j,k)**2*(1.-cnt*f(6,i,j,k,q3))**2
     &                *f(6,i,j,k,dng)/f(6,i,j,k,apw)
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
            ap(i,j,k)=aw(i,j,k)+ae(i,j,k)+as(i,j,k)
     &                +an(i,j,k)+ab(i,j,k)+at(i,j,k)-sp(i,j,k)
          enddo
        enddo
      enddo
!
!=====CELL FACE VELOCITIES
!
!
!-----WEST
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            dxwp=xc(i)-xc(i-1)
            volw=(1.-cnt*f(1,i,j,k,q3))*area(1,i,j,k)*dxwp
            dpxw=(p(i,j,k)-p(i-1,j,k))/dxwp
            dpxwa=(dpx(i-1,j,k)+dpx(i,j,k))/2
            vel(1,i,j,k)=f(1,i,j,k,ug)
     &                   -volw/f(1,i,j,k,apu)*(dpxw-dpxwa)
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
            dxpe=xc(i+1)-xc(i)
            vole=(1.-cnt*f(2,i,j,k,q3))*area(2,i,j,k)*dxpe
            dpxe=(p(i+1,j,k)-p(i,j,k))/dxpe
            dpxea=(dpx(i,j,k)+dpx(i+1,j,k))/2
            vel(2,i,j,k)=f(2,i,j,k,ug)
     &                   -vole/f(2,i,j,k,apu)*(dpxe-dpxea)
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
            dysp=yc(j)-yc(j-1)
            vols=(1.-cnt*f(3,i,j,k,q3))*area(3,i,j,k)*dysp
            dpys=(p(i,j,k)-p(i,j-1,k))/dysp
            dpysa=(dpy(i,j-1,k)+dpy(i,j,k))/2
            vel(3,i,j,k)=f(3,i,j,k,vg)
     &                   -vols/f(3,i,j,k,apv)*(dpys-dpysa)
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
            dypn=yc(j+1)-yc(j)
            voln=(1.-cnt*f(4,i,j,k,q3))*area(4,i,j,k)*dypn
            dpyn=(p(i,j+1,k)-p(i,j,k))/dypn
            dpyna=(dpy(i,j,k)+dpy(i,j+1,k))/2
            vel(4,i,j,k)=f(4,i,j,k,vg)
     &                   -voln/f(4,i,j,k,apv)*(dpyn-dpyna)
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
            dzbp=zc(k)-zc(k-1)
            volb=(1.-cnt*f(5,i,j,k,q3))*area(5,i,j,k)*dzbp
            dpzb=(p(i,j,k)-p(i,j,k-1))/dzbp
            dpzba=(dpz(i,j,k-1)+dpz(i,j,k))/2
            vel(5,i,j,k)=f(5,i,j,k,wg)
     &                   -volb/f(5,i,j,k,apw)*(dpzb-dpzba)
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
            dzpt=zc(k+1)-zc(k)
            volt=(1.-cnt*f(6,i,j,k,q3))*area(6,i,j,k)*dzpt
            dpzt=(p(i,j,k+1)-p(i,j,k))/dzpt
            dpzta=(dpz(i,j,k)+dpz(i,j,k+1))/2
            vel(6,i,j,k)=f(6,i,j,k,wg)
     &                   -volt/f(6,i,j,k,apw)*(dpzt-dpzta)
            if(k.eq.kmx)vel(6,i,j,k)=wg(i,j,k+1)
          enddo
        enddo
      enddo
!
!=====CONTINUITY SOURCE TERM
!
      res(iequ)=0
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
            if(tdcoef.lt.0)rdt=(1.5*ff(i,j,k)-2*ff0(i,j,k)
     &                         +0.5*ff00(i,j,k))
     &                         *vol(i,j,k)/delt
            su(i,j,k)=fw(i,j,k)-fe(i,j,k)+fs(i,j,k)
     &                -fn(i,j,k)+fb(i,j,k)-ft(i,j,k)-rdt
            res(iequ)=res(iequ)+abs(su(i,j,k))
          enddo
        enddo
      enddo
!
!-----SOLVE PRESSURE CORRECTION EQUATION
!
      call solve(pc)
!
!=====UPDATE PRESSURE
!
      ppref=pc(imx,jmx,kmx)
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            p(i,j,k)=p(i,j,k)+urf(12)*(pc(i,j,k)-ppref)
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
!=====UPDATE VELOCITIES
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            dpx(i,j,k)=(f(2,i,j,k,p)-f(1,i,j,k,p))/dx(i)
            dppx=(f(2,i,j,k,pc)-f(1,i,j,k,pc))/dx(i)
            ug(i,j,k)=ug(i,j,k)-(1.-cnt*q3(i,j,k))*vol(i,j,k)
     &                *dppx/apu(i,j,k)
!
            dpy(i,j,k)=(f(4,i,j,k,p)-f(3,i,j,k,p))/dy(j)
            dppy=(f(4,i,j,k,pc)-f(3,i,j,k,pc))/dy(j)
            vg(i,j,k)=vg(i,j,k)-(1.-cnt*q3(i,j,k))*vol(i,j,k)
     &                *dppy/apv(i,j,k)
!
            dpz(i,j,k)=(f(6,i,j,k,p)-f(5,i,j,k,p))/dz(k)
            dppz=(f(6,i,j,k,pc)-f(5,i,j,k,pc))/dz(k)
            wg(i,j,k)=wg(i,j,k)-(1.-cnt*q3(i,j,k))*vol(i,j,k)
     &                *dppz/apw(i,j,k)
          enddo
        enddo
      enddo

      return
      end