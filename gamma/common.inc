!
!=====COMMON
!
      parameter(ia=40,ja=20,ka=20)

      common/grd/
     &  xc(ia),x(ia),dx(ia)
     & ,yc(ja),y(ja),dy(ja)
     & ,zc(ka),z(ka),dz(ka)
     & ,fxw(ia),fxe(ia)
     & ,fys(ja),fyn(ja)
     & ,fzb(ka),fzt(ka)
     & ,area(6,ia,ja,ka)
     & ,vol(ia,ja,ka)
     & ,imn,imx,jmn,jmx,kmn,kmx
     & ,iinj,jinj,kinj
     & ,imon,jmon,kmon

      common/coeff/
     &  fw(ia,ja,ka),fe(ia,ja,ka)
     & ,fs(ia,ja,ka),fn(ia,ja,ka)
     & ,fb(ia,ja,ka),ft(ia,ja,ka)
     & ,dw(ia,ja,ka),de(ia,ja,ka)
     & ,ds(ia,ja,ka),dn(ia,ja,ka)
     & ,db(ia,ja,ka),dt(ia,ja,ka)
     & ,aw(ia,ja,ka),ae(ia,ja,ka)
     & ,as(ia,ja,ka),an(ia,ja,ka)
     & ,ab(ia,ja,ka),at(ia,ja,ka)
     & ,sp(ia,ja,ka),su(ia,ja,ka)
     & ,ap(ia,ja,ka),flux(ia,ja,ka)
     & ,apu(ia,ja,ka),apv(ia,ja,ka),apw(ia,ja,ka)
     & ,dpx(ia,ja,ka),dpy(ia,ja,ka),dpz(ia,ja,ka)
     & ,p(ia,ja,ka),pc(ia,ja,ka)

      common/vel/
     &  ug(ia,ja,ka),ug0(ia,ja,ka)
     & ,vg(ia,ja,ka),vg0(ia,ja,ka)
     & ,wg(ia,ja,ka),wg0(ia,ja,ka)
     & ,ul1(ia,ja,ka),ul10(ia,ja,ka)
     & ,vl1(ia,ja,ka),vl10(ia,ja,ka)
     & ,wl1(ia,ja,ka),wl10(ia,ja,ka)
     & ,ul2(ia,ja,ka),ul20(ia,ja,ka)
     & ,vl2(ia,ja,ka),vl20(ia,ja,ka)
     & ,wl2(ia,ja,ka),wl20(ia,ja,ka)
     & ,ul3(ia,ja,ka),ul30(ia,ja,ka)
     & ,vl3(ia,ja,ka),vl30(ia,ja,ka)
     & ,wl3(ia,ja,ka),wl30(ia,ja,ka)

      common/mom/
     &  qm1(ia,ja,ka),q0(ia,ja,ka)
     & ,r21(ia,ja,ka),r32(ia,ja,ka)
     & ,pk(ia,ja,ka)
     & ,bq1(ia,ja,ka),bq2(ia,ja,ka)
     & ,q1(ia,ja,ka),q10(ia,ja,ka)
     & ,q2(ia,ja,ka),q20(ia,ja,ka)
     & ,q3(ia,ja,ka),q30(ia,ja,ka)

      common/prop/
     &  ff(ia,ja,ka),ff0(ia,ja,ka),dd(ia,ja,ka)
     & ,dng(ia,ja,ka),dng0(ia,ja,ka)
     & ,dnl(ia,ja,ka),dnl0(ia,ja,ka)
     & ,den(ia,ja,ka),den0(ia,ja,ka)
     & ,dvsg(ia,ja,ka),dvsl(ia,ja,ka)
     & ,hcvg(ia,ja,ka),hcvl(ia,ja,ka)
     & ,hcpg(ia,ja,ka),hcpl(ia,ja,ka)
     & ,tcg(ia,ja,ka),tcl(ia,ja,ka)
     & ,tg(ia,ja,ka),tl(ia,ja,ka)
     & ,st(ia,ja,ka)

      common/icnt/
     &  itimx,ito,itomx,nts,ntsmx,npt,nradic
     & ,ivel,iph,iq,iequ,idrag,ibreak,icoll,ii
     & ,iquick

      common/rcnt/
     &  ulin,vlin,wlin
     & ,ugin,vgin,wgin
     & ,q1in,q2in,q3in
     & ,r32in
     & ,pkin,pkmin,pkmax
     & ,tgin,tlin
     & ,delt,tmax,time
     & ,order,ratio
     & ,cd,dinj,xl,yl,zl
     & ,beta,tdcoef
     & ,tiny,great,pi,cnt
     & ,rinj,rs,xinj,pinj
     & ,urf(20),res(20),gam(20)
     & ,sgas,sliq,sormx,valrat
     & ,gsmax0,gsmax,gdelta
