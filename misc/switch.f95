!
!=INTEGRAL FUNCTION
!
  program qalt

1 read(*,*)power
  
  if(power.lt.-1.999 .and.power.gt.-2.001)then
    s=1
  elseif(power.lt.-0.999 .and.power.gt.-1.001)then
    s=2
  elseif(power.gt.-0.001 .and.power.lt.0.001)then
    s=3
  elseif(power.gt.0.999 .and.power.lt.1.001)then
    s=4
  elseif(power.gt.1.999 .and.power.lt.2.001)then
    s=5
  elseif(power.gt.2.999 .and.power.lt.3.001)then
    s=6
  elseif(power.gt.3.999 .and.power.lt.4.001)then
    s=7
  elseif(power.gt.4.999 .and.power.lt.5.001)then
    s=8
  elseif(power.gt.-0.5)then
    s=9
  else
    s=10
  endif

  print*,s

  if(power.lt.10)goto 1

  end