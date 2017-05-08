!
!=====OUTPUT
!
      subroutine opt

      include "common.inc"
      
15    format('variables=
     &  "z","r"
     &  ,"pk","r32"
     &  ,"ug","vg","p","q3","ul3","vl3"
     &  ,"ul12","vl12","ul23","vl23"
     &  ,"drg1a2a","drg1b2b","drg2a3a","drg2b3b"
     &  ,"drgu1u2","drgv1v2","drgu2u3","drgv2v3"')
25    format(22(e15.5))
35    format('zone f=point, i=',i4,',j=',i4)
      open(10,file='tecplot.dat')
      write(10,15)
      write(10,35)imx-imn+1,jmx-jmn+1
      do j=jmn,jmx
        do i=imn,imx
          write(10,25)
     &      zc(i),rc(j)
     &      ,pk(i,j),r32(i,j)
     &      ,ug(i,j),vg(i,j),p(i,j)
     &      ,q3(i,j),ul3(i,j),vl3(i,j)
     &      ,ul1(i,j)/(ul2(i,j)+tiny)
     &      ,vl1(i,j)/(vl2(i,j)+tiny)
     &      ,ul2(i,j)/(ul3(i,j)+tiny)
     &      ,vl2(i,j)/(vl3(i,j)+tiny)
     &      ,drg1a(i,j)/(drg2a(i,j)+1.e+4*tiny)
     &      ,drg1b(i,j)/(drg2b(i,j)+1.e+4*tiny)
     &      ,drg2a(i,j)/(drg3a(i,j)+1.e+4*tiny)
     &      ,drg2b(i,j)/(drg3b(i,j)+1.e+4*tiny)
     &      ,drgu1(i,j)/(drgu2(i,j)+1.e+4*tiny)
     &      ,drgv1(i,j)/(drgv2(i,j)+1.e+4*tiny)
     &      ,drgu2(i,j)/(drgu3(i,j)+1.e+4*tiny)
     &      ,drgv2(i,j)/(drgv3(i,j)+1.e+4*tiny)
        enddo
      enddo
      close(10)
          
      return
      end