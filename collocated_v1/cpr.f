!
!=====GAS PRESSURE
!
      subroutine cpr(ip)
      include "common.inc"

      ivel=0
      iph=1
      iq=3
      call zro2
!
!-----SOLVE GAS PSEUDO-VELOCITIES
!
      if(ip.eq.2)call gpv
!
!-----SOLVE: SURFACE AREA / VELOCITY MAIN COEFFICIENT
!
      do i=imn,imx
        do j=jmn,jmx
          du(i,j)=vol(i,j)/(dxs(i)*apu(i,j))
          dv(i,j)=vol(i,j)/(dys(j)*apv(i,j))
        enddo
      enddo
      if(ip.eq.1)call zro4(pc)
!
!-----SOLVE COEFFICIENTS
!
      call var
      call fdc(du,dv)
      do i=imn,imx
        do j=jmn,jmx
          aw(i,j)=fw(i,j)
          ae(i,j)=fe(i,j)
          as(i,j)=fs(i,j)
          an(i,j)=fn(i,j)
        enddo
      enddo
      do ib=1,4
        if(ip.eq.1)call ncs(ib,pc)
        if(ip.eq.2)call ncs(ib,p)
      enddo
!
!-----SET PC TO ZERO AT INJECTOR CELL(S)
!
      i=iinj
      do j=jinj,jinj+ninjc-1
        sp(i,j)=-great
      enddo
!
!-----SET ZERO GRADIENT AT CELLS SURROUNDING INJECTOR
!
      do j=jinj,jinj+ninjc-1
        aw(i+1,j)=0.
        ae(i-1,j)=0.
      enddo
      j=jinj
      an(i,j-1)=0.
      j=jinj+ninjc-1
      as(i,j+1)=0.
!
!-----MAIN COEFFICIENT
!
      do i=imn,imx
        do j=jmn,jmx
          ap(i,j)=fw(i,j)+fe(i,j)+fs(i,j)+fn(i,j)-sp(i,j)
        enddo
      enddo
!
!-----CONITINUTY IMBALANCE TERM USING RHIE-CHOW VELOCITIES
!
      if(ip.eq.1)call rcc(ugx,vgx)
      if(ip.eq.2)call rcc(ugxx,vgxx)

      if(ip.eq.1)call sol1(pc)
      if(ip.eq.2)call sol1(p)

      return
      end