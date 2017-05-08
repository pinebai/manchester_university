!
!=====GAS PRESSURE AND VELOCITIES
!
      subroutine cuv
      include "common.inc"
      
      do i=imn,imx
        do j=jmn,jmx
!
!-----CORRECT U - COMPONENT
!
          ug(i,j)=ug(i,j)-du(i,j)*(pc(i+1,j)-pc(i-1,j))
          ug(i,j)=ur*ug(i,j)+(1-ur)*ugx(i,j)
          ugx(i,j)=ug(i,j)
!
!-----CORRECT V - COMPONENT
!
          vg(i,j)=vg(i,j)-dv(i,j)*(pc(i,j+1)-pc(i,j-1))
          vg(i,j)=ur*vg(i,j)+(1-ur)*vgx(i,j)
          vgx(i,j)=vg(i,j)
!
!-----CORRECT PRESSURE
!
          pcref=pc(imx,jmx)
          p(i,j)=p(i,j)+ur*(pc(i,j)-pcref)
        enddo
      enddo

      return
      end