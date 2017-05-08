!
!=====OUTPUT: FLOW
!
      subroutine output1
      include "common.inc"
      character*12 datafile
      
15    format('variables=&
        &"X(mm)","Y(mm)","Z(mm)"&
        &,"UG","UQ3"&
        &,"VG","VQ3"&
        &,"WG","WQ3"&
        &,"P","R32","Q3"')
25    format(12(e15.5))
35    format('zone f=point, i=',i4,',j=',i4,',k=',i4)
      if(nts.eq.1)ii=1000
      if(mod(nts,npt).eq.0.or.nts.eq.ntsmx)then
        ii=ii+1
        datafile='tecp0000.dat'
        write(datafile(5:8),'(i4.4)')ii
        open(ii,file=datafile)
        write(ii,15)
        write(ii,35)((imx+1)-(imn-1)+1)&
          ,((jmx+1)-(jmn-1)+1)&
          ,((kmx+1)-(kmn-1)+1)
        do k=kmn-1,kmx+1
          do j=jmn-1,jmx+1
            do i=imn-1,imx+1
              q3log=0
              if(q3(i,j,k).gt.tiny)q3log=-log10(q3(i,j,k)/q3in)
              write(ii,25)&
                xc(i)*1000,yc(j)*1000,zc(k)*1000&
                ,ug(i,j,k),uq3(i,j,k)&
                ,vg(i,j,k),vq3(i,j,k)&
                ,wg(i,j,k),wq3(i,j,k)&
                ,p(i,j,k),r32(i,j,k),q3log
            enddo
          enddo
        enddo
        close(ii)
      endif
      
      return
      end
!
!=====OUTPUT: PENETRATION
!
      subroutine penetration
      include "common.inc"
      
15    format(2(e15.5))
      do i=1,ia
        q3tot(i)=0.
      enddo
!
      do i=imn,imx
        do j=jmn,jmx
          do k=kmn,kmx
            q3tot(i)=q3(i,j,k)/area(1,i,j,k)+q3tot(i)
          enddo
        enddo
        if(q3tot(i).lt.(q3in/(pi*rs**2))*1.e-4)q3tot(i)=0.
      enddo
!
      do i=imx,imn,-1
        if(q3tot(i).gt.tiny)then
          xpen=x(i)*1000.
          if(xpen.lt.tiny)xpen=0.
          write(6,'('' PENETRATION (mm) = '',F6.2)')xpen
          write(6,*)'-----------------------------------------------'
          goto 100
        endif
      enddo
100   continue
!
      return
      end