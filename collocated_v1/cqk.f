!
!=====DERIVED MOMENTS AND GAMMA FUNCTION PARAMETER
!
      subroutine cqk
      include "common.inc"
      
      do i=imn,imx
        do j=jmn,jmx
!
!-----ALL MOMENTS MOST CO-EXIST
!
          if(q1(i,j).gt.tiny
     &      .and.q2(i,j).gt.tiny
     &      .and.q3(i,j).gt.1.e-7)then
!
!-----SOLVE R21, R32, K
!
            r21(i,j)=q2(i,j)/(q1(i,j)+tiny)
            r32(i,j)=q3(i,j)/(q2(i,j)+tiny)
            pk(i,j)=(1-2*r21(i,j)/r32(i,j))
     &              /(r21(i,j)/r32(i,j)-1)
!
!-----APPLY CONDITIONS ON K
!
            if(pk(i,j).gt.1
     &        .and.pk(i,j).lt.pkmin)then
              pk(i,j)=pkmin+tiny
            endif
            if(pk(i,j).gt.pkmax)then
              pk(i,j)=pkmax-tiny
            endif
            if(pk(i,j).gt.pkmin
     &        .and.pk(i,j).lt.pkmax)then
!
!-----SOLVE Q0, QM1
!
              q0(i,j)=q1(i,j)*(pk(i,j)+2)
     &                /(r32(i,j)*pk(i,j))
              qm1(i,j)=q0(i,j)*(pk(i,j)+2)
     &                 /(r32(i,j)*(pk(i,j)-1))
            else
              pk(i,j)=0.
              r21(i,j)=0.
              r32(i,j)=0.
              q3(i,j)=0.
              q2(i,j)=0.
              q1(i,j)=0.
              q0(i,j)=0.
              qm1(i,j)=0.
            endif
          else
            pk(i,j)=0.
            r21(i,j)=0.
            r32(i,j)=0.
            q3(i,j)=0.
            q2(i,j)=0.
            q1(i,j)=0.
            q0(i,j)=0.
            qm1(i,j)=0.
          endif
        enddo
      enddo
      
      return
      end