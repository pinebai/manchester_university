!
!=====SAVE PREVIOUS TIME STEP VALUES
!
      subroutine store
      include "common.inc"
      
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
!
!-----GAS PHASE MOMENTUM
!
            ug0(i,j,k)=ug(i,j,k)
            vg0(i,j,k)=vg(i,j,k)
            wg0(i,j,k)=wg(i,j,k)
!
!-----Q0 MOMENTUM
!
            uq00(i,j,k)=uq0(i,j,k)
            vq00(i,j,k)=vq0(i,j,k)
            wq00(i,j,k)=wq0(i,j,k)
!
!-----Q1 MOMENTUM
!
            uq10(i,j,k)=uq1(i,j,k)
            vq10(i,j,k)=vq1(i,j,k)
            wq10(i,j,k)=wq1(i,j,k)
!
!-----Q2 MOMENTUM
!
            uq20(i,j,k)=uq2(i,j,k)
            vq20(i,j,k)=vq2(i,j,k)
            wq20(i,j,k)=wq2(i,j,k)
!
!-----Q3 MOMENTUM
!
            uq30(i,j,k)=uq3(i,j,k)
            vq30(i,j,k)=vq3(i,j,k)
            wq30(i,j,k)=wq3(i,j,k)
!
!-----Q4 MOMENTUM
!
            uq40(i,j,k)=uq4(i,j,k)
            vq40(i,j,k)=vq4(i,j,k)
            wq40(i,j,k)=wq4(i,j,k)
!
!-----Q5 MOMENTUM
!
            uq50(i,j,k)=uq5(i,j,k)
            vq50(i,j,k)=vq5(i,j,k)
            wq50(i,j,k)=wq5(i,j,k)
!
!-----Q4 MOMENTUM
!
            uq60(i,j,k)=uq6(i,j,k)
            vq60(i,j,k)=vq6(i,j,k)
            wq60(i,j,k)=wq6(i,j,k)
!
!-----Q5 MOMENTUM
!
            uq70(i,j,k)=uq7(i,j,k)
            vq70(i,j,k)=vq7(i,j,k)
            wq70(i,j,k)=wq7(i,j,k)
!
!-----MOMENTS
!
            q00(i,j,k)=q0(i,j,k)
            q10(i,j,k)=q1(i,j,k)
            q20(i,j,k)=q2(i,j,k)
            q30(i,j,k)=q3(i,j,k)
            q40(i,j,k)=q4(i,j,k)
            q50(i,j,k)=q5(i,j,k)
            q60(i,j,k)=q6(i,j,k)
            q70(i,j,k)=q7(i,j,k)
!
!-----DENSITIES
!
            dng0(i,j,k)=dng(i,j,k)
            dnl0(i,j,k)=dnl(i,j,k)
          enddo
        enddo
      enddo

      return
      end