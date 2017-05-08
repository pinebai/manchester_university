program numbers

100 continue
read*,r
i=int(abs(r))
j=nint(abs(r))
if(j.eq.i)j=j+1
i=sign(real(i),r)
j=sign(real(j),r)
d=abs(r-i)
print*,i,1-d,j,d
! read*,i
if(abs(r).lt.10)goto 100
end