!
!=====COMMON
!
      parameter(ia=100,ja=50)

      common/grid/
     & zc(ia),dzs(ia),rc(ja),drs(ja)
     & ,dasw(ia,ja),dase(ia,ja)
     & ,dass(ia,ja),dasn(ia,ja)
     & ,dvs(ia,ja)
     & ,imin,imax,jmin,jmax
     & ,imn,imx,jmn,jmx,iinj,jinj

      common/coeff/
     & fw(ia,ja),fe(ia,ja),fs(ia,ja),fn(ia,ja)
     & ,dw(ia,ja),de(ia,ja),ds(ia,ja),dn(ia,ja)
     & ,aw(ia,ja),ae(ia,ja),as(ia,ja),an(ia,ja)
     & ,awu(ia,ja),aeu(ia,ja),asu(ia,ja),anu(ia,ja)
     & ,awv(ia,ja),aev(ia,ja),asv(ia,ja),anv(ia,ja)
     & ,suu(ia,ja),suv(ia,ja)
     & ,apu(ia,ja),apv(ia,ja)
     & ,sp(ia,ja),su(ia,ja)
     & ,ap(ia,ja),ap0(ia,ja)
     & ,p(ia,ja),pc(ia,ja)
     & ,du(ia,ja),dv(ia,ja)

      common/vel/
     & ug(ia,ja),ugx(ia,ja),ugxx(ia,ja),ug0(ia,ja)
     & ,vg(ia,ja),vgx(ia,ja),vgxx(ia,ja),vg0(ia,ja)
     & ,ul1(ia,ja),ul10(ia,ja)
     & ,vl1(ia,ja),vl10(ia,ja)
     & ,ul2(ia,ja),ul20(ia,ja)
     & ,vl2(ia,ja),vl20(ia,ja)
     & ,ul3(ia,ja),ul30(ia,ja)
     & ,vl3(ia,ja),vl30(ia,ja)
     & ,drg(ia,ja)
     & ,drgu1(ia,ja),drgu2(ia,ja),drgu3(ia,ja)
     & ,drgv1(ia,ja),drgv2(ia,ja),drgv3(ia,ja)
     & ,drg1a(ia,ja),drg2a(ia,ja),drg3a(ia,ja)
     & ,drg1b(ia,ja),drg2b(ia,ja),drg3b(ia,ja)

      common/mom/
     & r21(ia,ja),r32(ia,ja)
     & ,qm1(ia,ja),q0(ia,ja)
     & ,q1(ia,ja),q10(ia,ja)
     & ,q2(ia,ja),q20(ia,ja)
     & ,q3(ia,ja),q30(ia,ja)
     & ,qq(ia,ja),qq0(ia,ja)
     & ,pk(ia,ja)

      common/prop/
     & dng(ia,ja),dng0(ia,ja)
     & ,dnl(ia,ja),dnl0(ia,ja)
     & ,den(ia,ja),den0(ia,ja)
     & ,dvsg(ia,ja),dvsl(ia,ja)
     & ,dvis(ia,ja)
     & ,hcpg(ia,ja),hcvg(ia,ja)
     & ,tg(ia,ja),tl(ia,ja)

      common/icnt/
     & itc,itcmx
     & ,ntsmx,npt,ninjc
     & ,ivel,iph,iq,iap

      common/rcnt/
     & tiny,great,ur,pi,cnt
     & ,delt,tmax
     & ,ulin,vlin
     & ,q1in,q2in,q3in
     & ,pkmin,pkmax,pkin
     & ,dxc,dyc,xl,yl
     & ,dzc,drc,zl,rl

      common/lcnt/
     & qk1,qk2

      logical qk1,qk2