!
!=====OUTPUT: FLOW
!
      subroutine output1
      include "common.inc"
      character*11 datafile
      
15    format('variables=
     &  "X(mm)","Y(mm)","Z(mm)"
     &  ,"UG","UL3"
     &  ,"VG","VL3"
     &  ,"WG","WL3"
     &  ,"P","Q3"
     &  ,"R21","R32","PK"')
25    format(14(e15.5))
35    format('zone f=point, i=',i4,',j=',i4,',k=',i4)
      if(nts.eq.1)ii=900
      if(mod(nts,npt).eq.0.or.nts.eq.ntsmx)then
        ii=ii+1
        datafile='tecp000.dat'
        write(datafile(5:7),'(i3.3)')ii
        open(ii,file=datafile)
        write(ii,15)
        write(ii,35)((imx+1)-(imn-1)+1)
     &    ,((jmx+1)-(jmn-1)+1)
     &    ,((kmx+1)-(kmn-1)+1)
        do k=kmn-1,kmx+1
          do j=jmn-1,jmx+1
            do i=imn-1,imx+1
            write(ii,25)
     &        xc(i)*1000,yc(j)*1000,zc(k)*1000
     &        ,ug(i,j,k),ul3(i,j,k)
     &        ,vg(i,j,k),vl3(i,j,k)
     &        ,wg(i,j,k),wl3(i,j,k)
     &        ,p(i,j,k),q3(i,j,k)
     &        ,r21(i,j,k),r32(i,j,k),pk(i,j,k)
            enddo
          enddo
        enddo
        close(ii)
      endif
      
      return
      end