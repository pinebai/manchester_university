!
!======INITIALISATION
!
      block data
      include "common.inc"
!
!-----DURATION
!
      data tmax/1.e-4/
!
!-----ITERATION PARAMETERS
!
      data itimx,itomx,sormx/2,150,1.e-4/
!
!-----UNDER RELAXTION FACTORS
!
      data urf(1),urf(2),urf(3),urf(4),urf(5)/0.1,0.1,0.1,0.1,0.1/
      data urf(6),urf(7),urf(8),urf(9),urf(10)/0.1,0.1,0.1,0.1,0.1/
      data urf(11),urf(12),urf(13),urf(14),urf(15)/0.1,0.1,0.1,0.1,0.1/
      data urf(16),urf(18),urf(19),urf(20)/0.1,0.05,0.05,0.05/
!
!-----DISTRIBUTION PARAMETER VALUES
!
      data pkin,pkmin,pkmax/3.5,1.05,15./
      data r32in/15.e-6/
!
!-----CONE HALF-ANGLE
!
      data beta/6./
!
!-----INJECTOR RADIUS AND LENGTH
!
      data rinj/1.e-4/
      data xinj/5.e-4/
!
!-----NUMBER OF RADIAL INJECTION CELLS [1,2,3]
!
      data nradic/1/
!
!-----INJECTION TEMPERATURE
!
      data tlin/300./
!
!-----INJECTOR PRESSURE
!
      data pinj/10.e+6/
!
!-----INJECTOR DRAG COEFFICIENT
!
      data cd/0.7/
!
!-----TIME DISCRETISATION COEFFICIENT
!
      data tdcoef/1./
!
!-----CONSTANTS
!
      data ibreak,icoll/0,1/
      data iquick/1/
      data npt/1/
      data tiny,great/1.e-30,1.e+30/
      
      end