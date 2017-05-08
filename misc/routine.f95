program routine
real a,b
read*,a
call solve(a,b)
print*,a,b
end

subroutine solve(aa,bb)
real aa,bb
bb=aa**2
end