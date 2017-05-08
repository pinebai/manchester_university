!
!=====BREAK-UP MODEL
!
      subroutine breakup
      include "common.inc"
      order=1.e-4
!
!-----STRIPPING CONSTANT
!
      cs=1.0
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            term1=0;term2=0;term3=0;term4=0
            call dmom2(i,j,k)
            if(pk.gt.1.and.i.gt.imn.and.time.gt.0.5e-5)then
!
!=====RELATIVE VELOCITY
!
              velu=ul3(i,j,k)-ug(i,j,k)
              velv=vl3(i,j,k)-vg(i,j,k)
              velw=wl3(i,j,k)-wg(i,j,k)
              vrel=(velu**2+velv**2+velw**2)**0.5
!
!=====INTEGRAL LIMITS
!
!
!-----FOR STRIPPING BREAK-UP
!
              rs=st(i,j,k)**2/(2*dng(i,j,k)*vrel**3*dvsg(i,j,k))
              rslb=rs
              rsub=50.e-5
!
              xslb=(pk+2)*(rslb/r32)
              if(xslb.lt.tiny)xslb=0
              if(xslb.gt.10)xslb=10
!
              xsub=(pk+2)*(rsub/r32)
              if(xsub.lt.tiny)xslb=0
              if(xsub.gt.10)xslb=10
!
!-----FOR BAG BREAK-UP
!
              rb=3*st(i,j,k)/(dng(i,j,k)*vrel**2)
              rblb=rb
              rbub=rs
!
              xblb=(pk+2)*(rblb/r32)
              if(xblb.lt.tiny)xblb=0
              if(xblb.gt.10)xblb=10
!
              xbub=(pk+2)*(rbub/r32)
              if(xbub.lt.tiny)xbub=0
              if(xbub.gt.10)xbub=10
!
!=====TRUNCATED MOMENTS
!
!
!-----FOR STRIPPING BREAK-UP
!
              q0s=qt(0,xslb,xsub)
              q1s=qt(2,xslb,xsub)
              q1p5s=qt(3,xslb,xsub)
!
!-----FOR BAG BREAK-UP
!
              qm0p5b=qt(-1,xblb,xbub)
              q0p5b=qt(1,xblb,xbub)
!
!=====COMMON BREAK-UP TERMS
!
              cterm1=6.2
              cterm2=cs/vrel
              cterm3=dnl(i,j,k)/dng(i,j,k)
              cterm4=dvsl(i,j,k)/(2*dnl(i,j,k)*vrel)
              cterm5=pi*(dnl(i,j,k)/(2*st(i,j,k)))**0.5
!
!=====Q1 MOMENT
!
!
!-----STRIPPING BREAK-UP
!
              if(iq.eq.1)then
                term1=q1s/(cterm1**2*cterm2*cterm3*cterm4)
                term2=q0s/(cterm2*cterm3**0.5)
!
!-----BAG BREAK-UP
!
                term3=3*qm0p5b/cterm5
!
!=====Q2 MOMENT
!
!
!-----STRIPPING BREAK-UP
!
              elseif(iq.eq.2)then
                term1=q1p5s/(cterm1*cterm2*cterm3**0.75*cterm4**0.5)
                term2=q1s/(cterm2*cterm3**0.5)
!
!-----BAG BREAK-UP
!
                term3=q0p5b/cterm5
              endif
!
!=====BREAK-UP SOURCE TERM
!
              term4=dnl(i,j,k)*(term1-term2+term3)*vol(i,j,k)
              if(term1.gt.term2)su(i,j,k)=su(i,j,k)+term4
              
              if(iq.eq.1)then
                write(6,'(''t1,t2,t3,su,aw*q'',5(e11.3))')
     &            term1,term2,term3,term4,aw(i,j,k)*q1(i-1,j,k)
              elseif(iq.eq.2)then
                write(6,'(''t1,t2,t3,su,aw*q'',5(e11.3))')
     &            term1,term2,term3,term4,aw(i,j,k)*q2(i-1,j,k)
              endif
              
            endif
          enddo
        enddo
      enddo

      return
      end