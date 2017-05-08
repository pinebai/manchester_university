!
!=====INITIALISATION
!
      block data
      include "common.inc"
!
!-----DURATION
!
      data tmax/0.5e-3/
!
!-----TIME STEP
!
      data delt/2.e-6/
!
!-----DOMAIN SIZE
!
      data xl,yl,zl/0.024,0.004,0.004/
!
!-----ITERATION PARAMETERS
!
      data itimx,itomx,sormx,order/2,2,1.e-4,1.e-6/
!
!-----DISTRIBUTION PARAMETER VALUES
!
      data iqmax,r32in,pkin,rexpd/6,15.e-6,2.0,0./
!
!-----INJECTOR PRESSURE
!
      data pinj/9.9e+6/
!
!-----INJECTION TEMPERATURE
!
      data tlin/300./
!
!-----BREAK-UP, COLLISIONS & QUICK SCHEME SWITCHES [0=OFF]
!
      data ibreak,icoll,iquick,idrag/0,0,1,0/
!
!-----CONE HALF-ANGLE AND NOZZLE DRAG COEFFICIENT
!
      data beta,cd/5.,0.7/
!
!-----INJECTOR RADIUS AND LENGTH
!
      data rinj,xinj/1.5e-4,5.e-4/
!
!-----NUMBER OF RADIAL INJECTION CELLS [1,2,3]
!
      data nradic/1/
!
!-----TIME DISCRETISATION COEFFICIENT
!
      data tdcoef/1./
!
!-----CONSTANTS
!
      data npt/1/
      data tiny,great/1.e-30,1.e+30/
!
!-----UNDER RELAXTION FACTORS
!
      data urf(2)/0.8/
      
      end