!
!=====MAIN ARRAY ZERO
!
      subroutine zeromain
      include "common.inc"
      
      do i=1,ia
        do j=1,ja
          do k=1,ka
            unity(i,j,k)=1
            
            awu(i,j,k)=0
            aeu(i,j,k)=0
            asu(i,j,k)=0
            anu(i,j,k)=0
            abu(i,j,k)=0
            atu(i,j,k)=0
            apu(i,j,k)=0
            suu(i,j,k)=0
            uug(i,j,k)=0
            
            awv(i,j,k)=0
            aev(i,j,k)=0
            asv(i,j,k)=0
            anv(i,j,k)=0
            abv(i,j,k)=0
            atv(i,j,k)=0
            apv(i,j,k)=0
            suv(i,j,k)=0
            vvg(i,j,k)=0
            
            aww(i,j,k)=0
            aew(i,j,k)=0
            asw(i,j,k)=0
            anw(i,j,k)=0
            abw(i,j,k)=0
            atw(i,j,k)=0
            apw(i,j,k)=0
            suw(i,j,k)=0
            wwg(i,j,k)=0
          enddo
        enddo
      enddo
!
      do i=1,ia
        do j=1,ja
          do k=1,ka
            qa0(i,j,k)=0
            qa1(i,j,k)=0
            qa2(i,j,k)=0
            qa3(i,j,k)=0
            qa4(i,j,k)=0
            qa5(i,j,k)=0
            qa6(i,j,k)=0
            qa7(i,j,k)=0
!
            qb0(i,j,k)=0
            qb1(i,j,k)=0
            qb2(i,j,k)=0
            qb3(i,j,k)=0
            qb4(i,j,k)=0
            qb5(i,j,k)=0
            qb6(i,j,k)=0
            qb7(i,j,k)=0
!
            qm2(i,j,k)=0
            qm1(i,j,k)=0
            bq0(i,j,k)=0
            bq1(i,j,k)=0
            bq2(i,j,k)=0
            bq4(i,j,k)=0
            bq5(i,j,k)=0
            bq6(i,j,k)=0
            bq7(i,j,k)=0
!
            r32(i,j,k)=0
!
            flux(i,j,k)=0
            p(i,j,k)=0
!
            dw(i,j,k)=0
            de(i,j,k)=0
            ds(i,j,k)=0
            dn(i,j,k)=0
            db(i,j,k)=0
            dt(i,j,k)=0
!
            dpx(i,j,k)=0
            dpy(i,j,k)=0
            dpz(i,j,k)=0
            apu(i,j,k)=great
            apv(i,j,k)=great
            apw(i,j,k)=great
!
!-----CURRENT TIME STEP
!
            ug(i,j,k)=0
            vg(i,j,k)=0
            wg(i,j,k)=0
!
            uq0(i,j,k)=0
            vq0(i,j,k)=0
            wq0(i,j,k)=0
!
            uq1(i,j,k)=0
            vq1(i,j,k)=0
            wq1(i,j,k)=0
!
            uq2(i,j,k)=0
            vq2(i,j,k)=0
            wq2(i,j,k)=0
!
            uq3(i,j,k)=0
            vq3(i,j,k)=0
            wq3(i,j,k)=0
!
            uq4(i,j,k)=0
            vq4(i,j,k)=0
            wq4(i,j,k)=0
!
            uq5(i,j,k)=0
            vq5(i,j,k)=0
            wq5(i,j,k)=0
!
            uq6(i,j,k)=0
            vq6(i,j,k)=0
            wq6(i,j,k)=0
!
            uq7(i,j,k)=0
            vq7(i,j,k)=0
            wq7(i,j,k)=0
!
            q0(i,j,k)=0
            q1(i,j,k)=0
            q2(i,j,k)=0
            q3(i,j,k)=0
            q4(i,j,k)=0
            q5(i,j,k)=0
            q6(i,j,k)=0
            q7(i,j,k)=0
!
            dng(i,j,k)=0
            dnl(i,j,k)=0
!
!-----PREVIOUS TIME STEP
!
            ug0(i,j,k)=0
            vg0(i,j,k)=0
            wg0(i,j,k)=0
!
            uq00(i,j,k)=0
            vq00(i,j,k)=0
            wq00(i,j,k)=0
!
            uq10(i,j,k)=0
            vq10(i,j,k)=0
            wq10(i,j,k)=0
!
            uq20(i,j,k)=0
            vq20(i,j,k)=0
            wq20(i,j,k)=0
!
            uq30(i,j,k)=0
            vq30(i,j,k)=0
            wq30(i,j,k)=0
!
            uq40(i,j,k)=0
            vq40(i,j,k)=0
            wq40(i,j,k)=0
!
            uq50(i,j,k)=0
            vq50(i,j,k)=0
            wq50(i,j,k)=0
!
            uq60(i,j,k)=0
            vq60(i,j,k)=0
            wq60(i,j,k)=0
!
            uq70(i,j,k)=0
            vq70(i,j,k)=0
            wq70(i,j,k)=0
!
            q00(i,j,k)=0
            q10(i,j,k)=0
            q20(i,j,k)=0
            q30(i,j,k)=0
            q40(i,j,k)=0
            q50(i,j,k)=0
            q60(i,j,k)=0
            q70(i,j,k)=0
!
            dng0(i,j,k)=0
            dnl0(i,j,k)=0
          enddo
        enddo
      enddo

      do iequ=1,20
        res(iequ)=0
      enddo
      
      gsmax0=0;gsmax=0;gdelta=0

      return
      end
!
!=====INTERMEDIATE ARRAY ZERO
!
      subroutine zerocoef
      include "common.inc"
      
      do i=1,ia
        do j=1,ja
          do k=1,ka
            aw(i,j,k)=0
            ae(i,j,k)=0
            as(i,j,k)=0
            an(i,j,k)=0
            ab(i,j,k)=0
            at(i,j,k)=0
            ap(i,j,k)=0
            su(i,j,k)=0
            sp(i,j,k)=0
          enddo
        enddo
      enddo

      return
      end
!
!=====ZERO PARTICULAR ARRAY
!
      subroutine zerovar(phi)
      include "common.inc"
      real phi(ia,ja,ka)
      
      do i=1,ia
        do j=1,ja
          do k=1,ka
            phi(i,j,k)=0
          enddo
        enddo
      enddo

      return
      end
!
!=====ZERO BOUNDARY COEFFICIENTS
!
      subroutine zerobcoef
      include "common.inc"
!
!-----WEST
!
        i=imn
        do j=jmn,jmx
          do k=kmn,kmx
            aw(i,j,k)=0
          enddo
        enddo
!
!-----EAST
!
        i=imx
        do j=jmn,jmx
          do k=kmn,kmx
            ae(i,j,k)=0
          enddo
        enddo
!
!-----SOUTH
!
        j=jmn
        do i=imn,imx
          do k=kmn,kmx
            as(i,j,k)=0
          enddo
        enddo
!
!-----NORTH
!
        j=jmx
        do i=imn,imx
          do k=kmn,kmx
           an(i,j,k)=0
          enddo
        enddo
!
!-----BOTTOM
!
        k=kmn
        do i=imn,imx
          do j=jmn,jmx
            ab(i,j,k)=0
          enddo
        enddo
!
!-----TOP
!
        k=kmx
        do i=imn,imx
          do j=jmn,jmx
            at(i,j,k)=0
          enddo
        enddo

      return
      end