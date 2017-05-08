!
!=====COMMON BLOCKS
!
      parameter(ia=35,ja=20,ka=20)

      common/grdcb/
     &  xc(ia),x(ia),dx(ia)
     & ,yc(ja),y(ja),dy(ja)
     & ,zc(ka),z(ka),dz(ka)
!
     & ,fxw(ia),fxe(ia)
     & ,fys(ja),fyn(ja)
     & ,fzb(ka),fzt(ka)
!
     & ,area(6,ia,ja,ka)
     & ,vol(ia,ja,ka)
!
     & ,imn,imx,jmn,jmx,kmn,kmx
     & ,iinj,jinj,kinj
     & ,imon,jmon,kmon

      common/coeffcb/
     &  fw(ia,ja,ka),fe(ia,ja,ka)
     & ,fs(ia,ja,ka),fn(ia,ja,ka)
     & ,fb(ia,ja,ka),ft(ia,ja,ka)
!
     & ,dw(ia,ja,ka),de(ia,ja,ka)
     & ,ds(ia,ja,ka),dn(ia,ja,ka)
     & ,db(ia,ja,ka),dt(ia,ja,ka)
!
     & ,aw(ia,ja,ka),ae(ia,ja,ka)
     & ,as(ia,ja,ka),an(ia,ja,ka)
     & ,ab(ia,ja,ka),at(ia,ja,ka)
     & ,ap(ia,ja,ka)
!
     & ,sp(ia,ja,ka),su(ia,ja,ka)
!
     & ,dpx(ia,ja,ka),dpy(ia,ja,ka),dpz(ia,ja,ka)
     & ,pc(ia,ja,ka),p(ia,ja,ka),flux(ia,ja,ka)
     
      common/vcoeffcb/
     &  awu(ia,ja,ka),aeu(ia,ja,ka)
     & ,asu(ia,ja,ka),anu(ia,ja,ka)
     & ,abu(ia,ja,ka),atu(ia,ja,ka)
     & ,apu(ia,ja,ka),suu(ia,ja,ka)
!
     & ,awv(ia,ja,ka),aev(ia,ja,ka)
     & ,asv(ia,ja,ka),anv(ia,ja,ka)
     & ,abv(ia,ja,ka),atv(ia,ja,ka)
     & ,apv(ia,ja,ka),suv(ia,ja,ka)
!
     & ,aww(ia,ja,ka),aew(ia,ja,ka)
     & ,asw(ia,ja,ka),anw(ia,ja,ka)
     & ,abw(ia,ja,ka),atw(ia,ja,ka)
     & ,apw(ia,ja,ka),suw(ia,ja,ka)

      common/velcb/
     &  ug(ia,ja,ka),uug(ia,ja,ka),ug0(ia,ja,ka)
     & ,vg(ia,ja,ka),vvg(ia,ja,ka),vg0(ia,ja,ka)
     & ,wg(ia,ja,ka),wwg(ia,ja,ka),wg0(ia,ja,ka)
!
     & ,uq0(ia,ja,ka),uq00(ia,ja,ka)
     & ,vq0(ia,ja,ka),vq00(ia,ja,ka)
     & ,wq0(ia,ja,ka),wq00(ia,ja,ka)
!
     & ,uq1(ia,ja,ka),uq10(ia,ja,ka)
     & ,vq1(ia,ja,ka),vq10(ia,ja,ka)
     & ,wq1(ia,ja,ka),wq10(ia,ja,ka)
!
     & ,uq2(ia,ja,ka),uq20(ia,ja,ka)
     & ,vq2(ia,ja,ka),vq20(ia,ja,ka)
     & ,wq2(ia,ja,ka),wq20(ia,ja,ka)
!
     & ,uq3(ia,ja,ka),uq30(ia,ja,ka)
     & ,vq3(ia,ja,ka),vq30(ia,ja,ka)
     & ,wq3(ia,ja,ka),wq30(ia,ja,ka)
!
     & ,uq4(ia,ja,ka),uq40(ia,ja,ka)
     & ,vq4(ia,ja,ka),vq40(ia,ja,ka)
     & ,wq4(ia,ja,ka),wq40(ia,ja,ka)
!
     & ,uq5(ia,ja,ka),uq50(ia,ja,ka)
     & ,vq5(ia,ja,ka),vq50(ia,ja,ka)
     & ,wq5(ia,ja,ka),wq50(ia,ja,ka)
!
     & ,uq6(ia,ja,ka),uq60(ia,ja,ka)
     & ,vq6(ia,ja,ka),vq60(ia,ja,ka)
     & ,wq6(ia,ja,ka),wq60(ia,ja,ka)
!
     & ,uq7(ia,ja,ka),uq70(ia,ja,ka)
     & ,vq7(ia,ja,ka),vq70(ia,ja,ka)
     & ,wq7(ia,ja,ka),wq70(ia,ja,ka)

      common/momcb/
     &  bq0(ia,ja,ka),bq1(ia,ja,ka),bq2(ia,ja,ka)
     & ,bq4(ia,ja,ka),bq5(ia,ja,ka),bq6(ia,ja,ka)
     & ,bq7(ia,ja,ka)
!
     & ,r32(ia,ja,ka)
!
     & ,qm1(ia,ja,ka),qm2(ia,ja,ka)
!
     & ,q0(ia,ja,ka),q00(ia,ja,ka)
     & ,q1(ia,ja,ka),q10(ia,ja,ka)
     & ,q2(ia,ja,ka),q20(ia,ja,ka)
     & ,q3(ia,ja,ka),q30(ia,ja,ka)
     & ,q4(ia,ja,ka),q40(ia,ja,ka)
     & ,q5(ia,ja,ka),q50(ia,ja,ka)
     & ,q6(ia,ja,ka),q60(ia,ja,ka)
     & ,q7(ia,ja,ka),q70(ia,ja,ka)
!
     & ,qa0(ia,ja,ka),qb0(ia,ja,ka)
     & ,qa1(ia,ja,ka),qb1(ia,ja,ka)
     & ,qa2(ia,ja,ka),qb2(ia,ja,ka)
     & ,qa3(ia,ja,ka),qb3(ia,ja,ka)
     & ,qa4(ia,ja,ka),qb4(ia,ja,ka)
     & ,qa5(ia,ja,ka),qb5(ia,ja,ka)
     & ,qa6(ia,ja,ka),qb6(ia,ja,ka)
     & ,qa7(ia,ja,ka),qb7(ia,ja,ka)

      common/propcb/
     &  ff(ia,ja,ka),ff0(ia,ja,ka),dd(ia,ja,ka)
!
     & ,dng(ia,ja,ka),dng0(ia,ja,ka)
     & ,dnl(ia,ja,ka),dnl0(ia,ja,ka)
     & ,den(ia,ja,ka),den0(ia,ja,ka)
!
     & ,dvsg(ia,ja,ka),dvsl(ia,ja,ka)
     & ,hcvg(ia,ja,ka),hcvl(ia,ja,ka)
     & ,hcpg(ia,ja,ka),hcpl(ia,ja,ka)
     & ,tcg(ia,ja,ka),tcl(ia,ja,ka)
     & ,tg(ia,ja,ka),tl(ia,ja,ka)
!
     & ,st(ia,ja,ka)

      common/icntcb/
     &  itimx,ito,itomx,nts,ntsmx,npt,nradic
     & ,ivel,iph,iq,iequ,idrag,ibreak,icoll,ii
     & ,iquick,iqmax

      common/rcntcb/
     &  uqin,vqin,wqin,ruqin
     & ,ugin,vgin,wgin,rugin
     & ,q0in,q1in,q2in,q3in,q4in,q5in,q6in,q7in
!
     & ,r32in,rexpd,rexpm
     & ,pkin,pkmin,pkmax
!
     & ,tgin,tlin
     & ,delt,tmax,time
     & ,order
     & ,cd,dinj,xl,yl,zl
     & ,beta,tdcoef
     & ,tiny,great,pi,cnt
     & ,rinj,rs,xinj,pinj
     & ,urf(20),res(20),resn(20),gam(20),q3tot(ia)
     & ,fluxn1,fluxn2
     & ,sgas,sliq,sormx
     & ,gsmax0,gsmax,gdelta,urfc