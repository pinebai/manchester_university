!------------------------------------------------
! DEFINE ARRAYS
!------------------------------------------------
!
MODULE mod_all
  LOGICAL,PARAMETER::t=.true.,f=.false.
  LOGICAL::lgas=t,lspray=t,lrebound=t,lpolar=t,ldrag=t,lbreakup=f,lcollision=f
  INTEGER,PARAMETER::nnx=25,nny=25
  INTEGER::nx=nnx,ny=nny,nts=300
  REAL::delt=1.5e-6,pi=3.141593,order=1.e-6
  REAL::xinj=5.e-4,yinj=5.e-4
  REAL::xl=5.e-3,yl=5.e-3
  REAL::uin=50.,vin=5.
  REAL::q0in,q1in,q2in,q3in
  REAL,DIMENSION(nnx,nny),TARGET:: &
! GAS
  uo,vo,u,v,q3so,q3s, &
! INCOMING SPRAY
  u0o,v0o,u0,v0,q0o,q0, &
  u1o,v1o,u1,v1,q1o,q1, &
  u2o,v2o,u2,v2,q2o,q2, &
  u3o,v3o,u3,v3,q3o,q3, &
! REBOUNDING SPRAY
  u0io,v0io,u0i,v0i,q0io,q0i, &
  u1io,v1io,u1i,v1i,q1io,q1i, &
  u2io,v2io,u2i,v2i,q2io,q2i, &
  u3io,v3io,u3i,v3i,q3io,q3i
  REAL,DIMENSION(nnx,nny):: &
! COEFFICIENTS
  ffo,ff,dd,vol,apu,apv,dpx,dpy,pp, &
  fw,fe,fs,fn,dw,de,ds,dn,aw,ae,as,an,ap,su,sp, &
! GAS PROPERTIES
  deno_g,den_g,vis_g,t_g,cp_g,cv_g,p, &
! LIQUID PROPERTIES
  deno_l,den_l,vis_l,t_l,cp_l,cv_l,st, &
! SOURCES
  du0,dv0,du1,dv1,du2,dv2,du3,dv3,bq0,bq1,bq2,cq0,cq1,cq2, &
  du0i,dv0i,du1i,dv1i,du2i,dv2i,du3i,dv3i,bq0i,bq1i,bq2i,cq0i,cq1i,cq2i
! GRID
  REAL,DIMENSION(nnx)::x,xc,dx,fxw,fxe
  REAL,DIMENSION(nny)::y,yc,dy,fys,fyn
  REAL,DIMENSION(4,nnx,nny)::area
END MODULE mod_all
!------------------------------------------------
! DEFINE INTERFACE MODULES
!------------------------------------------------
!
MODULE mod_int_f1
  INTERFACE
    FUNCTION f1(ib,i,j,psi,psi_p)
      USE mod_all
      REAL,DIMENSION(:,:),OPTIONAL::psi
      REAL,DIMENSION(:,:),POINTER,OPTIONAL::psi_p
      REAL,DIMENSION(nnx,nny)::phi
    END FUNCTION f1
  END INTERFACE
END MODULE mod_int_f1
MODULE mod_int_f2
  INTERFACE
    FUNCTION f2(ib,i,j,psi1,psi2,psi1_p,psi2_p)
      USE mod_all
      REAL,DIMENSION(:,:),OPTIONAL::psi1,psi2
      REAL,DIMENSION(:,:),POINTER,OPTIONAL::psi1_p,psi2_p
      REAL,DIMENSION(nnx,nny)::phi1,phi2
    END FUNCTION f2
  END INTERFACE
END MODULE mod_int_f2
MODULE mod_int_dbc
  INTERFACE
    SUBROUTINE dbc(ib,var,psi,psi_p)
      USE mod_all
      REAL,DIMENSION(:,:),OPTIONAL::psi
      REAL,DIMENSION(:,:),POINTER,OPTIONAL::psi_p
      REAL,DIMENSION(nnx,nny)::phi
    END SUBROUTINE dbc
  END INTERFACE
END MODULE mod_int_dbc
MODULE mod_int_nbc
  INTERFACE
    SUBROUTINE nbc(ib,psi,psi_p)
      USE mod_all
      REAL,DIMENSION(:,:),OPTIONAL::psi
      REAL,DIMENSION(:,:),POINTER,OPTIONAL::psi_p
      REAL,DIMENSION(nnx,nny)::phi
    END SUBROUTINE nbc
  END INTERFACE
END MODULE mod_int_nbc
MODULE mod_int_coeffs
  INTERFACE
    SUBROUTINE coeffs(qo_p,q_p,phi_p)
      USE mod_all
      REAL,DIMENSION(:,:),POINTER::qo_p,q_p,phi_p
    END SUBROUTINE coeffs
  END INTERFACE
END MODULE mod_int_coeffs
MODULE mod_int_convection
  INTERFACE
    SUBROUTINE convection(u_p,v_p)
      USE mod_all
      USE mod_int_f2
      REAL,DIMENSION(:,:),POINTER::u_p,v_p
    END SUBROUTINE convection
  END INTERFACE
END MODULE mod_int_convection
MODULE mod_int_temporal
  INTERFACE
    SUBROUTINE temporal(phio_p,phi_p)
      USE mod_all
      REAL,DIMENSION(:,:),POINTER::phio_p,phi_p
    END SUBROUTINE temporal
  END INTERFACE
END MODULE mod_int_temporal
MODULE mod_int_solver
  INTERFACE
    SUBROUTINE solver(urf,psi,psi_p)
      USE mod_all
      REAL,DIMENSION(:,:),OPTIONAL::psi
      REAL,DIMENSION(:,:),POINTER,OPTIONAL::psi_p
      REAL,DIMENSION(nnx,nny)::aa,cc,phi
    END SUBROUTINE solver
  END INTERFACE
END MODULE mod_int_solver
MODULE mod_int_boundary
  INTERFACE
    SUBROUTINE boundary(phi_p)
      USE mod_all
      REAL,DIMENSION(:,:),POINTER::phi_p
    END SUBROUTINE boundary
  END INTERFACE
END MODULE mod_int_boundary
MODULE mod_int_source
  INTERFACE
    SUBROUTINE source(phi_p)
      USE mod_all
      REAL,DIMENSION(:,:),POINTER::phi_p
    END SUBROUTINE source
  END INTERFACE
END MODULE mod_int_source
MODULE mod_int_quick
  INTERFACE
    SUBROUTINE quick(phi_p)
      USE mod_all
      REAL,DIMENSION(:,:),POINTER::phi_p
    END SUBROUTINE quick
  END INTERFACE
END MODULE mod_int_quick
MODULE mod_int_drag
  INTERFACE
    FUNCTION drag(i,j,vel_g,vel_l,qa,qb)
      USE mod_all
      INTEGER,INTENT(IN)::i,j
      REAL,DIMENSION(:,:),INTENT(IN)::vel_g,vel_l,qa,qb
    END FUNCTION drag
  END INTERFACE
END MODULE mod_int_drag
MODULE mod_int_moment
  INTERFACE
    SUBROUTINE moment(i,j,power,rlb,rub,p0,p1,p2,p3,q)
      COMMON/lcoeff/b0,b1,b2,b3
      INTEGER,INTENT(IN)::i,j
      REAL,INTENT(IN)::power
      REAL,DIMENSION(:,:),INTENT(IN)::p0,p1,p2,p3
      REAL,INTENT(INOUT)::rlb,rub
      REAL,INTENT(OUT)::q
      REAL,EXTERNAL::flagu
    END SUBROUTINE moment
  END INTERFACE
END MODULE mod_int_moment
MODULE mod_int_breakup
  INTERFACE
    FUNCTION breakup(i,j,iq,rub,ud,vd,p0,p1,p2,p3)
      USE mod_all,ONLY:lbreakup,q2in,q3in,u,v,den_g,den_l,vis_l,st
      REAL,DIMENSION(:,:),INTENT(IN)::ud,vd,p0,p1,p2,p3
      INTEGER,INTENT(IN)::i,j,iq
      REAL,INTENT(INOUT)::rub
    END FUNCTION breakup
  END INTERFACE
END MODULE mod_int_breakup
MODULE mod_int_collision
  INTERFACE
    FUNCTION collision(i,j,iq,rub,ud3,vd3,ud0,vd0,p0,p1,p2,p3)
      USE mod_all,ONLY:lcollision,pi,den_g,st
      USE mod_int_moment
      REAL,DIMENSION(:,:),INTENT(IN)::ud3,vd3,ud0,vd0,p0,p1,p2,p3
      INTEGER,INTENT(IN)::i,j,iq
      REAL,INTENT(INOUT)::rub
    END FUNCTION collision
  END INTERFACE
END MODULE mod_int_collision
MODULE mod_int_wall
  INTERFACE
    FUNCTION wall(i,j,iq,lreb,lspl)
      USE mod_all
      USE mod_int_moment
      REAL,DIMENSION(15)::qs,wes
      INTEGER,INTENT(IN)::i,j,iq
      LOGICAL,INTENT(IN),OPTIONAL::lreb,lspl
    END FUNCTION wall
  END INTERFACE
END MODULE mod_int_wall
!------------------------------------------------
! PROGRAM
!------------------------------------------------
!
PROGRAM main
  USE mod_all
! INITIALISE
  pkin=3.; r32in=15.e-6; q3in=0.145
  q2in=q3in/r32in
  q1in=q2in*(pkin+2)/(r32in*(pkin+1))
  q0in=q1in*(pkin+2)/(r32in*pkin)
  CALL zero; CALL properties; CALL grid
! TIME LOOP
  DO n=1,nts
    PRINT*,n,u(2,2),u(2,3),v(2,2),v(2,3)
    CALL store
! ITERATION LOOP
    DO nn=1,2
      CALL calc
      CALL sources
    ENDDO
    CALL output
  ENDDO
! OUTPUT
  CALL output
END PROGRAM main
!------------------------------------------------
! SOLVE TRANSPORT EQUATIONS
!------------------------------------------------
! U_P    -> X COMPONENT VELOCITY ETC.
! QO_P   -> RELATED OLD MOMENT
! Q_P    -> RELATED MOMENT
! PHIO_P -> OLD VALUE OF ARGUMENT BEING SOLVED
! PHI_P  -> ARGUMENT BEING SOLVED
!
SUBROUTINE calc
  USE mod_all
  REAL,DIMENSION(:,:),POINTER::u_p,v_p,qo_p,q_p,phio_p,phi_p
! GAS PHASE
  IF(lgas)THEN
    q3so=q3o+q3io; q3s=q3+q3i
    NULLIFY(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    u_p=>u; v_p=>v; qo_p=>q3so; q_p=>q3s
    phio_p=>uo; phi_p=>u; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
    phio_p=>vo; phi_p=>v; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
    CALL calcp
  ENDIF
! SPRAY EXITING INJECTOR
! Q0
  IF(lspray)THEN
    NULLIFY(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    u_p=>u0; v_p=>v0; qo_p=>q0o; q_p=>q0
    phio_p=>u0o; phi_p=>u0; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
    phio_p=>v0o; phi_p=>v0; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
    phio_p=>qo_p; phi_p=>q_p; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
! Q1
    NULLIFY(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    u_p=>u1; v_p=>v1; qo_p=>q1o; q_p=>q1
    phio_p=>u1o; phi_p=>u1; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
    phio_p=>v1o; phi_p=>v1; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
    phio_p=>qo_p; phi_p=>q_p; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
! Q2
    NULLIFY(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    u_p=>u2; v_p=>v2; qo_p=>q2o; q_p=>q2
    phio_p=>u2o; phi_p=>u2; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
    phio_p=>v2o; phi_p=>v2; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
    phio_p=>qo_p; phi_p=>q_p; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
! Q3
    NULLIFY(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    u_p=>u3; v_p=>v3; qo_p=>q3o; q_p=>q3
    phio_p=>u3o; phi_p=>u3; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
    phio_p=>v3o; phi_p=>v3; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
    phio_p=>qo_p; phi_p=>q_p; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
  ENDIF
! SPRAY REBOUNDING OF EAST WALL
! Q0i
  IF(lrebound)THEN
    NULLIFY(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    u_p=>u0i; v_p=>v0i; qo_p=>q0io; q_p=>q0i
    phio_p=>u0io; phi_p=>u0i; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
    phio_p=>v0io; phi_p=>v0i; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
    phio_p=>qo_p; phi_p=>q_p; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
! Q1i
    NULLIFY(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    u_p=>u1i; v_p=>v1i; qo_p=>q1io; q_p=>q1i
    phio_p=>u1io; phi_p=>u1i; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
    phio_p=>v1io; phi_p=>v1i; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
    phio_p=>qo_p; phi_p=>q_p; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
! Q2i
    NULLIFY(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    u_p=>u2i; v_p=>v2i; qo_p=>q2io; q_p=>q2i
    phio_p=>u2io; phi_p=>u2i; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
    phio_p=>v2io; phi_p=>v2i; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
    phio_p=>qo_p; phi_p=>q_p; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
! Q3i
    NULLIFY(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    u_p=>u3i; v_p=>v3i; qo_p=>q3io; q_p=>q3i
    phio_p=>u3io; phi_p=>u3i; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
    phio_p=>v3io; phi_p=>v3i; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
    phio_p=>qo_p; phi_p=>q_p; CALL solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
    NULLIFY(phio_p,phi_p)
  ENDIF
  CONTAINS
    SUBROUTINE solve(u_p,v_p,qo_p,q_p,phio_p,phi_p)
      USE mod_all
      USE mod_int_coeffs
      USE mod_int_convection
      USE mod_int_quick
      USE mod_int_boundary
      USE mod_int_source
      USE mod_int_temporal
      USE mod_int_solver
      REAL,DIMENSION(:,:),POINTER::u_p,v_p,qo_p,q_p,phio_p,phi_p
      CALL zeromcoeff
      CALL coeffs(qo_p,q_p,phi_p)
      CALL convection(u_p,v_p)
      CALL diffusion
      CALL hybrid
      IF(ASSOCIATED(phi_p,q0).or.ASSOCIATED(phi_p,q1) &
      .or.ASSOCIATED(phi_p,q2).or.ASSOCIATED(phi_p,q3) &
      .or.ASSOCIATED(phi_p,q0i).or.ASSOCIATED(phi_p,q1i) &
      .or.ASSOCIATED(phi_p,q2i).or.ASSOCIATED(phi_p,q3i))CALL quick(phi_p)
      CALL boundary(phi_p)
      CALL source(phi_p)
      CALL temporal(phio_p,phi_p)
      CALL solver(urf=0.7,psi_p=phi_p)
      DO i=2,nx-1
        DO j=2,ny-1
          IF(ASSOCIATED(phi_p,u))apu(i,j)=1/ap(i,j)
          IF(ASSOCIATED(phi_p,v))apv(i,j)=1/ap(i,j)
        ENDDO
      ENDDO
    END SUBROUTINE solve
END SUBROUTINE calc
!------------------------------------------------
! COEFFICIENTS OF TRANSPORT TERMS
!------------------------------------------------
! FFO -> OLD FLUX COEFFICINET
! FF  -> FLUX COEFFICINET
! DD  -> DIFFUSION COEFFICINET
!
SUBROUTINE coeffs(qo_p,q_p,phi_p)
  USE mod_all
  REAL,DIMENSION(:,:),POINTER::qo_p,q_p,phi_p
  IF(ASSOCIATED(phi_p,u).or.ASSOCIATED(phi_p,v))THEN
    ffo=deno_g*(1.-(4/3.)*pi*qo_p); ff=den_g*(1.-(4/3.)*pi*q_p)
    dd=vis_g*(1.-(4/3.)*pi*q_p)
  ELSEIF(ASSOCIATED(phi_p,u0).or.ASSOCIATED(phi_p,v0) &
  .or.ASSOCIATED(phi_p,u1).or.ASSOCIATED(phi_p,v1) &
  .or.ASSOCIATED(phi_p,u2).or.ASSOCIATED(phi_p,v2) &
  .or.ASSOCIATED(phi_p,u0i).or.ASSOCIATED(phi_p,v0i) &
  .or.ASSOCIATED(phi_p,u1i).or.ASSOCIATED(phi_p,v1i) &
  .or.ASSOCIATED(phi_p,u2i).or.ASSOCIATED(phi_p,v2i))THEN
    ffo=qo_p; ff=q_p
    dd=0.7*vis_l*q_p/den_l
  ELSEIF(ASSOCIATED(phi_p,u3).or.ASSOCIATED(phi_p,v3) &
  .or.ASSOCIATED(phi_p,u3i).or.ASSOCIATED(phi_p,v3i))THEN
    ffo=(4/3.)*pi*deno_l*qo_p; ff=(4/3.)*pi*den_l*q_p
    dd=0.7*vis_l*(4/3.)*pi*q_p
  ELSEIF(ASSOCIATED(phi_p,q0).or.ASSOCIATED(phi_p,q1).or.ASSOCIATED(phi_p,q2) &
  .or.ASSOCIATED(phi_p,q0i).or.ASSOCIATED(phi_p,q1i).or.ASSOCIATED(phi_p,q2i))THEN
    ffo=1.; ff=1.
    dd=0.
  ELSEIF(ASSOCIATED(phi_p,q3).or.ASSOCIATED(phi_p,q3i))THEN
    ffo=deno_l; ff=den_l
    dd=0.
  ENDIF
END SUBROUTINE coeffs
!------------------------------------------------
! CONVECTION COEFFICIENTS
!------------------------------------------------
! FW -> WEST COEFFICIENT
! FE -> EAST COEFFICIENT
! FS -> SOUTH COEFFICIENT
! FN -> NORTH COEFFICIENT
!
SUBROUTINE convection(u_p,v_p)
  USE mod_all
  USE mod_int_f2
  REAL,DIMENSION(:,:),POINTER::u_p,v_p
  DO i=2,nx-1
    DO j=2,ny-1
      fw(i,j)=area(1,i,j)*f2(1,i,j,ff,psi2_p=u_p)
      fe(i,j)=area(2,i,j)*f2(2,i,j,ff,psi2_p=u_p)
      fs(i,j)=area(3,i,j)*f2(3,i,j,ff,psi2_p=v_p)
      fn(i,j)=area(4,i,j)*f2(4,i,j,ff,psi2_p=v_p)
    ENDDO
  ENDDO
END SUBROUTINE convection
!------------------------------------------------
! DIFFUSION COEFFICIENTS
!------------------------------------------------
! DW -> WEST COEFFICIENT
! DE -> EAST COEFFICIENT
! DS -> SOUTH COEFFICIENT
! DN -> NORTH COEFFICIENT
!
SUBROUTINE diffusion
  USE mod_all
  USE mod_int_f1
  DO i=2,nx-1
    DO j=2,ny-1
      dw(i,j)=area(1,i,j)*f1(1,i,j,dd)/(xc(i)-xc(i-1))
      de(i,j)=area(2,i,j)*f1(2,i,j,dd)/(xc(i+1)-xc(i))
      ds(i,j)=area(3,i,j)*f1(3,i,j,dd)/(yc(j)-yc(j-1))
      dn(i,j)=area(4,i,j)*f1(4,i,j,dd)/(yc(j+1)-yc(j))
    ENDDO
  ENDDO
END SUBROUTINE diffusion
!------------------------------------------------
! SPATIAL DIFFERENCING SCHEME: HYBRID
!------------------------------------------------
! AW -> WEST MAIN COEFFICIENT
! AE -> EAST MAIN COEFFICIENT
! AS -> SOUTH MAIN COEFFICIENT
! AN -> NORTH MAIN COEFFICIENT
!
SUBROUTINE hybrid
  USE mod_all
  DO i=2,nx-1
    DO j=2,ny-1
      aw(i,j)=max(fw(i,j),(dw(i,j)+fw(i,j)/2),0.)
      ae(i,j)=max(-fe(i,j),(de(i,j)-fe(i,j)/2),0.)
      as(i,j)=max(fs(i,j),(ds(i,j)+fs(i,j)/2),0.)
      an(i,j)=max(-fn(i,j),(dn(i,j)-fn(i,j)/2),0.)
    ENDDO
  ENDDO
END SUBROUTINE hybrid
!------------------------------------------------
! SPATIAL DIFFERENCING SCHEME: QUICK
!------------------------------------------------
!
SUBROUTINE quick(phi_p)
  USE mod_all
  REAL,DIMENSION(:,:),POINTER::phi_p
  DO i=2,nx-1
    DO j=2,ny-1
      IF(fw(i,j)<1.-30)qw=0.; IF(fw(i,j)>-1.-30)qw=1.
      IF(fe(i,j)<1.-30)qe=0.; IF(fe(i,j)>-1.-30)qe=1.
      IF(fs(i,j)<1.-30)qs=0.; IF(fs(i,j)>-1.-30)qs=1.
      IF(fn(i,j)<1.-30)qn=0.; IF(fn(i,j)>-1.-30)qn=1.
      aw(i,j)=dw(i,j)+qw*fw(i,j)
      ae(i,j)=de(i,j)-(1-qe)*fe(i,j)
      as(i,j)=ds(i,j)+qs*fs(i,j)
      an(i,j)=dn(i,j)-(1-qn)*fn(i,j)
! PHI(WW) AT THE WEST BOUNDARY IS LINEARLY INTERPOLATED FROM PHI(W) & PHI(P), ETC.
! WEST
      t1a=1/8.*(3*phi_p(i,j)-2*phi_p(i-1,j)-phi_p(i-2,j))*qw*fw(i,j)
      IF(i==2)t1a=1/8.*(3*phi_p(i,j)-3*phi_p(i-1,j))*qw*fw(i,j)
      t1b=1/8.*(phi_p(i-1,j)+2*phi_p(i,j)-3*phi_p(i+1,j))*qe*fe(i,j)
! EAST
      t2a=1/8.*(3*phi_p(i-1,j)-2*phi_p(i,j)-phi_p(i+1,j))*(1-qw)*fw(i,j)
      t2b=1/8.*(2*phi_p(i+1,j)+phi_p(i+2,j)-3*phi_p(i,j))*(1-qe)*fe(i,j)
      IF(i==nx-1)t2b=1/8.*(3*phi_p(i+1,j)-3*phi_p(i,j))*(1-qe)*fe(i,j)
! SOUTH
      t3a=1/8.*(3*phi_p(i,j)-2*phi_p(i,j-1)-phi_p(i,j-2))*qs*fs(i,j)
      IF(j==2)t3a=1/8.*(3*phi_p(i,j)-3*phi_p(i,j-1))*qs*fs(i,j)
      t3b=1/8.*(phi_p(i,j-1)+2*phi_p(i,j)-3*phi_p(i,j+1))*qn*fn(i,j)
! NORTH
      t4a=1/8.*(3*phi_p(i,j-1)-2*phi_p(i,j)-phi_p(i,j+1))*(1-qs)*fs(i,j)
      t4b=1/8.*(2*phi_p(i,j+1)+phi_p(i,j+2)-3*phi_p(i,j))*(1-qn)*fn(i,j)
      IF(j==ny-1)t4b=1/8.*(3*phi_p(i,j+1)-3*phi_p(i,j))*(1-qn)*fn(i,j)
! SOURCE TERM
      su(i,j)=su(i,j)+t1a+t1b+t2a+t2b+t3a+t3b+t4a+t4b
    ENDDO
  ENDDO
END SUBROUTINE quick
!------------------------------------------------
! TERMPORAL DIFFERENCING SCHEME: EULER IMPLICIT
!------------------------------------------------
! DELF -> FLUX BALANCE
! RDTO -> OLD TIME STEP COEFFICIENT
! RDT  -> TIME STEP COEFFICIENT
! DELT -> TIME STEP
!
SUBROUTINE temporal(phio_p,phi_p)
  USE mod_all
  REAL,DIMENSION(:,:),POINTER::phio_p,phi_p
  DO i=2,nx-1
    DO j=2,ny-1
      delf=fe(i,j)-fw(i,j)+fn(i,j)-fs(i,j)
      rdto=ffo(i,j)*vol(i,j)/delt
      rdt=ff(i,j)*vol(i,j)/delt
      IF(ASSOCIATED(phi_p,q0).or.ASSOCIATED(phi_p,q1) &
      .or.ASSOCIATED(phi_p,q2).or.ASSOCIATED(phi_p,q3) &
      .or.ASSOCIATED(phi_p,q0i).or.ASSOCIATED(phi_p,q1i) &
      .or.ASSOCIATED(phi_p,q2i).or.ASSOCIATED(phi_p,q3i))THEN
        ap(i,j)=aw(i,j)+ae(i,j)+as(i,j)+an(i,j)-sp(i,j)+rdt+delf
      ELSE
        ap(i,j)=aw(i,j)+ae(i,j)+as(i,j)+an(i,j)-sp(i,j)+rdto
      ENDIF
      su(i,j)=su(i,j)+rdto*phio_p(i,j)
    ENDDO
  ENDDO
END SUBROUTINE temporal
!------------------------------------------------
! SOLUTION ALGORITHM [ A.x = b]: TDMA
!------------------------------------------------
! URF -> UNDER-RELAXATION FACTOR
! AA  -> TEMPORY ARRAY
! CC  -> TEMPORY ARRAY
! D   -> TEMPORY DENOMINATOR
! S   -> TEMPORY SOURCE
!
SUBROUTINE solver(urf,psi,psi_p)
  USE mod_all
  REAL,DIMENSION(:,:),OPTIONAL::psi
  REAL,DIMENSION(:,:),POINTER,OPTIONAL::psi_p
  REAL,DIMENSION(nnx,nny)::aa=0.,cc=0.,phi=0.
  IF(PRESENT(psi))phi=psi
  IF(PRESENT(psi_p))phi=psi_p
! UNDER-RELAX
  DO i=2,nx-1
    DO j=2,ny-1
      IF(abs(sp(i,j))<1.e+30)THEN
        ap(i,j)=ap(i,j)/urf
        su(i,j)=su(i,j)+(1.-urf)*ap(i,j)*phi(i,j)
      ENDIF
    ENDDO
  ENDDO
  DO nn=1,2
! W-E LINES
    DO j=2,ny-1
      DO i=2,nx-1
        d=ap(i,j)-aw(i,j)*aa(i-1,j)+1.e-30
        s=su(i,j)+as(i,j)*phi(i,j-1)+an(i,j)*phi(i,j+1)
        IF(i==2)s=s+aw(i,j)*phi(i-1,j)
        IF(i==nx-1)s=s+ae(i,j)*phi(i+1,j)
        aa(i,j)=ae(i,j)/d
        cc(i,j)=(aw(i,j)*cc(i-1,j)+s)/d
      ENDDO
    ENDDO
    DO j=2,ny-1
      DO i=nx-1,2,-1
        phi(i,j)=phi(i+1,j)*aa(i,j)+cc(i,j)
      ENDDO
    ENDDO
! S-N LINES
    DO i=2,nx-1
      DO j=2,ny-1
        d=ap(i,j)-as(i,j)*aa(i-1,j)+1.e-30
        s=su(i,j)+aw(i,j)*phi(i-1,j)+ae(i,j)*phi(i+1,j)
        IF(j==2)s=s+as(i,j)*phi(i,j-1)
        IF(j==ny-1)s=s+an(i,j)*phi(i,j+1)
        aa(i,j)=an(i,j)/d
        cc(i,j)=(as(i,j)*cc(i,j-1)+s)/d
      ENDDO
    ENDDO
    DO i=2,nx-1
      DO j=ny-1,2,-1
        phi(i,j)=phi(i,j+1)*aa(i,j)+cc(i,j)
      ENDDO
    ENDDO
  ENDDO
  IF(PRESENT(psi))psi=phi
  IF(PRESENT(psi_p))psi_p=phi
END SUBROUTINE solver
!------------------------------------------------
! PRESSURE FIELD: COLLOCATED METHOD
!------------------------------------------------
! PP   -> PRESSURE CORRECTION
! P    -> PRESSURE
! DPX  -> PRESSURE GRADIENT IN X AXIS ETC.
! APU  -> 1/AP ETC.
! VELW -> WEST FACE VELOCITY ETC.
! Q3S  -> SUM OF Q3'S
!
SUBROUTINE calcp
  USE mod_all
  USE mod_int_f1
  USE mod_int_nbc
  USE mod_int_solver
  CALL zeromcoeff; pp=0.
! SOLVE FLUXES
  ffo=deno_g*(1.-(4/3.)*pi*q3so); ff=den_g*(1.-(4/3.)*pi*q3s)
! SOLVE COEFFICIENTS
  DO i=2,nx-1
    DO j=2,ny-1
      aw(i,j)=f1(1,i,j,den_g)*(1.-(4/3.)*pi*f1(1,i,j,q3s))**2 &
      *area(1,i,j)**2*f1(1,i,j,apu)
      ae(i,j)=f1(2,i,j,den_g)*(1.-(4/3.)*pi*f1(2,i,j,q3s))**2 &
      *area(2,i,j)**2*f1(2,i,j,apu)
      as(i,j)=f1(3,i,j,den_g)*(1.-(4/3.)*pi*f1(3,i,j,q3s))**2 &
      *area(3,i,j)**2*f1(3,i,j,apv)
      an(i,j)=f1(4,i,j,den_g)*(1.-(4/3.)*pi*f1(4,i,j,q3s))**2 &
      *area(4,i,j)**2*f1(4,i,j,apv)
      ap(i,j)=aw(i,j)+ae(i,j)+as(i,j)+an(i,j)-sp(i,j)
! CELL FACE VELOCITIES
! WEST
      dxwp=xc(i)-xc(i-1)
      volw=(1.-(4/3.)*pi*f1(1,i,j,q3s))*area(1,i,j)*dxwp
      dpxw=(p(i,j)-p(i-1,j))/dxwp
      dpxwa=(dpx(i-1,j)+dpx(i,j))/2
      velw=f1(1,i,j,u)-volw*f1(1,i,j,apu)*(dpxw-dpxwa)
! EAST
      dxpe=xc(i+1)-xc(i)
      vole=(1.-(4/3.)*pi*f1(2,i,j,q3s))*area(2,i,j)*dxpe
      dpxe=(p(i+1,j)-p(i,j))/dxpe
      dpxea=(dpx(i,j)+dpx(i+1,j))/2
      vele=f1(2,i,j,u)-vole*f1(2,i,j,apu)*(dpxe-dpxea)
! SOUTH
      dysp=yc(j)-yc(j-1)
      vols=(1.-(4/3.)*pi*f1(3,i,j,q3s))*area(3,i,j)*dysp
      dpys=(p(i,j)-p(i,j-1))/dysp
      dpysa=(dpy(i,j-1)+dpy(i,j))/2
      vels=f1(3,i,j,v)-vols*f1(3,i,j,apv)*(dpys-dpysa)
! NORTH
      dypn=yc(j+1)-yc(j)
      voln=(1.-(4/3.)*pi*f1(4,i,j,q3s))*area(4,i,j)*dypn
      dpyn=(p(i,j+1)-p(i,j))/dypn
      dpyna=(dpy(i,j)+dpy(i,j+1))/2
      veln=f1(4,i,j,v)-voln*f1(4,i,j,apv)*(dpyn-dpyna)
! CONTINUITY SOURCE TERM
      fw(i,j)=area(1,i,j)*f1(1,i,j,ff)*velw
      fe(i,j)=area(2,i,j)*f1(2,i,j,ff)*vele
      fs(i,j)=area(3,i,j)*f1(3,i,j,ff)*vels
      fn(i,j)=area(4,i,j)*f1(4,i,j,ff)*veln
      rdt=(ff(i,j)-ffo(i,j))*vol(i,j)/delt
      su(i,j)=fw(i,j)-fe(i,j)+fs(i,j)-fn(i,j)-rdt
    ENDDO
  ENDDO
! SOLVE PRESSURE CORRECTION EQUATION
  CALL solver(urf=1.,psi=pp)
! UPDATE PRESSURE
  urf=0.7
  DO i=2,nx-1
    DO j=2,ny-1
      p(i,j)=p(i,j)+urf*(pp(i,j)-pp(2,2))
    ENDDO
  ENDDO
! ZERO GRADIENT BOUNDARY VALUES
  DO nf=1,4
    CALL nbc(nf,pp); CALL nbc(nf,p)
  ENDDO
! UPDATE VELOCITIES AND PRESSURE GRADIENTS
  DO i=2,nx-1
    DO j=2,ny-1
      dppx=(f1(2,i,j,pp)-f1(1,i,j,pp))/dx(i)
      dpx(i,j)=(f1(2,i,j,p)-f1(1,i,j,p))/dx(i)
      u(i,j)=u(i,j)-(1.-(4/3.)*pi*q3s(i,j))*vol(i,j)*dppx*apu(i,j)
      dppy=(f1(4,i,j,pp)-f1(3,i,j,pp))/dy(j)
      dpy(i,j)=(f1(4,i,j,p)-f1(3,i,j,p))/dy(j)
      v(i,j)=v(i,j)-(1.-(4/3.)*pi*q3s(i,j))*vol(i,j)*dppy*apv(i,j)
    ENDDO
  ENDDO
END SUBROUTINE calcp
!------------------------------------------------
! GRID: DEFINED BY CELL FACE LOCATIONS
!------------------------------------------------
! X -> LENGTH
! Y -> WIDTH (RADIUS)
! DX, DY -> DIMENSIONS OF CELL
! XC, YC -> LOCATION OF CELL CENTRE
! FXW -> DISTANCE RATIO OF WEST FACE TO WEST NODE
! AREA -> SURFACE AREAS OF CELL FACES
! VOL -> VOLUME OF CELL
!
SUBROUTINE grid
  USE mod_all
  REAL,DIMENSION(nny)::r
! CELL FACE LOCATIONS DEFINED BY USER
  x(1)=0.; x(2)=xinj
  DO i=3,nx-1
    x(i)=x(i-1)+xinj
  ENDDO
  y(1)=0.; y(2)=yinj
  DO j=3,ny-1
    y(j)=y(j-1)+yinj
  ENDDO
  r=1.
  IF(lpolar)THEN
    r=y; PRINT*,'polar coordinates'
  ENDIF
! COORDINATES OF CELL CENTERS
  xc(1)=x(1)
  DO i=2,nx-1
    xc(i)=(x(i)+x(i-1))/2
  ENDDO
  xc(nx)=x(nx-1)
  yc(1)=y(1)
  DO j=2,ny-1
    yc(j)=(y(j)+y(j-1))/2
  ENDDO
  yc(ny)=y(ny-1)
! INTERPOLATION FACTORS
  DO i=2,nx-1
    fxw(i)=1-(x(i-1)-xc(i-1))/(xc(i)-xc(i-1))
    fxe(i)=1-(xc(i+1)-x(i))/(xc(i+1)-xc(i))
  ENDDO
  DO j=2,ny-1
    fys(j)=1-(y(j-1)-yc(j-1))/(yc(j)-yc(j-1))
    fyn(j)=1-(yc(j+1)-y(j))/(yc(j+1)-yc(j))
  ENDDO
! LENGTHS
  DO i=2,nx-1
    dx(i)=x(i)-x(i-1)
  ENDDO
  DO j=2,ny-1
    dy(j)=y(j)-y(j-1)
  ENDDO
! AREAS
  DO i=2,nx-1
    DO j=2,ny-1
      area(1,i,j)=dy(j)*(r(j)+r(j-1))/2
      area(2,i,j)=dy(j)*(r(j)+r(j-1))/2
      area(3,i,j)=dx(i)*r(j-1)
      area(4,i,j)=dx(i)*r(j)
    ENDDO
  ENDDO
! VOLUMES
  DO i=2,nx-1
    DO j=2,ny-1
      vol(i,j)=dx(i)*dy(j)*(r(j)+r(j-1))/2
    ENDDO
  ENDDO
END SUBROUTINE grid
!------------------------------------------------
! CELL FACE INTERPOLATED VALUES
!------------------------------------------------
! F1 -> PHI INTERPOLATED AT BOUNDARY IB
! IB -> BOUNDARY INDEX; WEST = 1, EAST =2,...
!
FUNCTION f1(ib,i,j,psi,psi_p)
  USE mod_all
  REAL,DIMENSION(:,:),OPTIONAL::psi
  REAL,DIMENSION(:,:),POINTER,OPTIONAL::psi_p
  REAL,DIMENSION(nnx,nny)::phi
  IF(PRESENT(psi))phi=psi
  IF(PRESENT(psi_p))phi=psi_p
! WEST
  IF(ib==1)THEN
    f1=fxw(i)*phi(i-1,j)+(1-fxw(i))*phi(i,j)
! EAST
  ELSEIF(ib==2)THEN
    f1=(1-fxe(i))*phi(i,j)+fxe(i)*phi(i+1,j)
! SOUTH
  ELSEIF(ib==3)THEN
    f1=fys(j)*phi(i,j-1)+(1-fys(j))*phi(i,j)
! NORTH
  ELSEIF(ib==4)THEN
    f1=(1-fyn(j))*phi(i,j)+fyn(j)*phi(i,j+1)
  ENDIF
END FUNCTION f1
!------------------------------------------------
! CELL FACE TWO INTERPOLATED VALUES
!------------------------------------------------
! F2 -> PHI1*PHI2 INTERPOLATED AT BOUNDARY IB
! IB -> BOUNDARY INDEX; WEST = 1, EAST =2,...
!
FUNCTION f2(ib,i,j,psi1,psi2,psi1_p,psi2_p)
  USE mod_all
  REAL,DIMENSION(:,:),OPTIONAL::psi1,psi2
  REAL,DIMENSION(:,:),POINTER,OPTIONAL::psi1_p,psi2_p
  REAL,DIMENSION(nnx,nny)::phi1,phi2
  IF(PRESENT(psi1))phi1=psi1
  IF(PRESENT(psi2))phi2=psi2
  IF(PRESENT(psi1_p))phi1=psi1_p
  IF(PRESENT(psi2_p))phi2=psi2_p
! WEST
  IF(ib==1)THEN
    f2=fxw(i)*phi1(i-1,j)*phi2(i-1,j)+(1-fxw(i))*phi1(i,j)*phi2(i,j)
! EAST
  ELSEIF(ib==2)THEN
    f2=(1-fxe(i))*phi1(i,j)*phi2(i,j)+fxe(i)*phi1(i+1,j)*phi2(i+1,j)
! SOUTH
  ELSEIF(ib==3)THEN
    f2=fys(j)*phi1(i,j-1)*phi2(i,j-1)+(1-fys(j))*phi1(i,j)*phi2(i,j)
! NORTH
  ELSEIF(ib==4)THEN
    f2=(1-fyn(j))*phi1(i,j)*phi2(i,j)+fyn(j)*phi1(i,j+1)*phi2(i,j+1)
  ENDIF
END FUNCTION f2
!------------------------------------------------
! DIRICHLET CONDITIONS
!------------------------------------------------
! VAR -> BOUNDARY VALUE
! IB -> BOUNDARY INDEX; WEST = 1, EAST =2,...
!
SUBROUTINE dbc(ib,var,psi,psi_p)
  USE mod_all
  REAL,DIMENSION(:,:),OPTIONAL::psi
  REAL,DIMENSION(:,:),POINTER,OPTIONAL::psi_p
  REAL,DIMENSION(nnx,nny)::phi
  IF(PRESENT(psi))phi=psi
  IF(PRESENT(psi_p))phi=psi_p
! WEST
  IF(ib==1)THEN
    DO j=2,ny-1
      phi(1,j)=var
    ENDDO
! EAST
  ELSEIF(ib==2)THEN
    DO j=2,ny-1
      phi(nx,j)=var
    ENDDO
! SOUTH
  ELSEIF(ib==3)THEN
    DO i=2,nx-1
      phi(i,1)=var
    ENDDO
! NORTH
  ELSEIF(ib==4)THEN
    DO i=2,nx-1
      phi(i,ny)=var
    ENDDO
  ENDIF
  IF(PRESENT(psi))psi=phi
  IF(PRESENT(psi_p))psi_p=phi
END SUBROUTINE dbc
!------------------------------------------------
! NEUMANN CONDITIONS
!------------------------------------------------
! IB -> BOUNDARY INDEX; WEST = 1, EAST =2,...
!
SUBROUTINE nbc(ib,psi,psi_p)
  USE mod_all
  REAL,DIMENSION(:,:),OPTIONAL::psi
  REAL,DIMENSION(:,:),POINTER,OPTIONAL::psi_p
  REAL,DIMENSION(nnx,nny)::phi
  IF(PRESENT(psi))phi=psi
  IF(PRESENT(psi_p))phi=psi_p
! WEST
  IF(ib==1)THEN
    DO j=2,ny-1
      phi(1,j)=phi(2,j)
    ENDDO
! EAST
  ELSEIF(ib==2)THEN
    DO j=2,ny-1
      phi(nx,j)=phi(nx-1,j)
    ENDDO
! SOUTH
  ELSEIF(ib==3)THEN
    DO i=2,nx-1
      phi(i,1)=phi(i,2)
    ENDDO
! NORTH
  ELSEIF(ib==4)THEN
    DO i=2,nx-1
      phi(i,ny)=phi(i,ny-1)
    ENDDO
  ENDIF
  IF(PRESENT(psi))psi=phi
  IF(PRESENT(psi_p))psi_p=phi
END SUBROUTINE nbc
!------------------------------------------------
! TRANSPORT AND THERMODYNAMIC PROPERTIES
!------------------------------------------------
! DEN_G -> GAS DENSITY
! VIS_G -> GAS DYNAMIC VISCOSITY
! T_G -> GAS TEMPERATURE
! CP_G -> GAS SPECIFIC HEAT CAPACITY AT CONSTANT PRESSURE
! CV_G -> GAS SPECIFIC HEAT CAPACITY AT CONSTANT VOLUME
! DEN_G -> LIQUID DENSITY
! VIS_G -> LIQUID DYNAMIC VISCOSITY 
! T_L -> LIQUID TEMPERATURE
! CP_L -> LIQUID SPECIFIC HEAT CAPACITY AT CONSTANT PRESSURE
! CV_L -> LIQUID SPECIFIC HEAT CAPACITY AT CONSTANT VOLUME
! ST -> SURFACE TENSION
!
SUBROUTINE properties
  USE mod_all
! GAS PROPERTIES
  den_g=13.2; deno_g=den_g; vis_g=0.00002
  cp_g=1004.9; cv_g=717.8; t_g=300.
  p=den_g*(cp_g-cv_g)*t_g
! LIQUID PROPERTIES
  den_l=840.; deno_l=den_l; vis_l=0.0027
  cp_l=1004.9; cv_l=717.8; t_l=300.
  st=0.0283
END SUBROUTINE properties
!------------------------------------------------
! ZERO MAIN COEFFICIENTS
!------------------------------------------------
!
SUBROUTINE zeromcoeff
  USE mod_all
  aw=0.; ae=0.; as=0.; an=0.; ap=0.; su=0.; sp=0.
END SUBROUTINE zeromcoeff
!------------------------------------------------
! STORE OLD VALUES [FROM PREVIOUS TIME STEP]
!------------------------------------------------
!
SUBROUTINE store
  USE mod_all
  uo=u; vo=v; q3so=q3s
  u0o=u0; v0o=v0; u1o=u1; v1o=v1; u2o=u2; v2o=v2; u3o=u3; v3o=v3
  q0o=q0; q1o=q1; q2o=q2; q3o=q3
  u0io=u0i; v0io=v0i; u1io=u1i; v1io=v1i; u2io=u2i; v2io=v2i; u3io=u3i; v3io=v3i
  q0io=q0i; q1io=q1i; q2io=q2i; q3io=q3i
  deno_g=den_g; deno_l=den_l
END SUBROUTINE store
!------------------------------------------------
! ZERO ARRAYS
!------------------------------------------------
!
SUBROUTINE zero
  USE mod_all
  u=0.; v=0.; q3s=0.; p=0.; dpx=0.; dpy=0.; apu=0.; apv=0.; vol=0.
  u0=0.; v0=0.; u1=0.; v1=0.; u2=0.; v2=0.; u3=0.; v3=0.
  bq0=0.; bq1=0.; bq2=0.; cq0=0.; cq1=0.; cq2=0.
  du0=0.; dv0=0.; du1=0.; dv1=0.; du2=0.; dv2=0.; du3=0.; dv3=0.
  q0=0.; q1=0.; q2=0.; q3=0.
  u0i=0.; v0i=0.; u1i=0.; v1i=0.; u2i=0.; v2i=0.; u3i=0.; v3i=0.
  bq0i=0.; bq1i=0.; bq2i=0.; cq0i=0.; cq1i=0.; cq2i=0.
  du0i=0.; dv0i=0.; du1i=0.; dv1i=0.; du2i=0.; dv2i=0.; du3i=0.; dv3i=0.
  q0i=0.; q1i=0.; q2i=0.; q3i=0.
END SUBROUTINE zero
!------------------------------------------------
! BOUNDARY CONDITIONS
!------------------------------------------------
!
SUBROUTINE boundary(phi_p)
  USE mod_all
  USE mod_int_dbc
  USE mod_int_nbc
  USE mod_int_wall
  REAL,DIMENSION(:,:),POINTER::phi_p
! GAS PHASE
  IF(ASSOCIATED(phi_p,u))THEN
    CALL nbc(1,psi_p=phi_p); CALL dbc(2,0.,psi_p=phi_p)
    CALL nbc(3,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
  ELSEIF(ASSOCIATED(phi_p,v))THEN
    CALL nbc(1,psi_p=phi_p); CALL dbc(2,0.,psi_p=phi_p)
    CALL dbc(3,0.,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    DO j=2,ny-1
      sp(nx-1,j)=sp(nx-1,j)-vis_g(nx-1,j)*area(2,nx-1,j)/(x(nx-1)-xc(nx-1))
    ENDDO
! INCOMING SPRAY
! Q0
  ELSEIF(ASSOCIATED(phi_p,u0))THEN
    CALL dbc(1,0.,psi_p=phi_p); CALL nbc(2,psi_p=phi_p)
    CALL nbc(3,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    phi_p(1,2)=uin; as(2,3)=0.
  ELSEIF(ASSOCIATED(phi_p,v0))THEN
    CALL dbc(1,0.,psi_p=phi_p); CALL nbc(2,psi_p=phi_p)
    CALL dbc(3,0.,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    phi_p(1,2)=vin; as(2,3)=0.
  ELSEIF(ASSOCIATED(phi_p,q0))THEN
    CALL dbc(1,0.,psi_p=phi_p); CALL nbc(2,psi_p=phi_p)
    CALL nbc(3,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    phi_p(1,2)=q0in; as(2,3)=0.
! Q1
  ELSEIF(ASSOCIATED(phi_p,u1))THEN
    CALL dbc(1,0.,psi_p=phi_p); CALL nbc(2,psi_p=phi_p)
    CALL nbc(3,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    phi_p(1,2)=uin; as(2,3)=0.
  ELSEIF(ASSOCIATED(phi_p,v1))THEN
    CALL dbc(1,0.,psi_p=phi_p); CALL nbc(2,psi_p=phi_p)
    CALL dbc(3,0.,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    phi_p(1,2)=vin; as(2,3)=0.
  ELSEIF(ASSOCIATED(phi_p,q1))THEN
    CALL dbc(1,0.,psi_p=phi_p); CALL nbc(2,psi_p=phi_p)
    CALL nbc(3,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    phi_p(1,2)=q1in; as(2,3)=0.
! Q2
  ELSEIF(ASSOCIATED(phi_p,u2))THEN
    CALL dbc(1,0.,psi_p=phi_p); CALL nbc(2,psi_p=phi_p)
    CALL nbc(3,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    phi_p(1,2)=uin; as(2,3)=0.
  ELSEIF(ASSOCIATED(phi_p,v2))THEN
    CALL dbc(1,0.,psi_p=phi_p); CALL nbc(2,psi_p=phi_p)
    CALL dbc(3,0.,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    phi_p(1,2)=vin; as(2,3)=0.
  ELSEIF(ASSOCIATED(phi_p,q2))THEN
    CALL dbc(1,0.,psi_p=phi_p); CALL nbc(2,psi_p=phi_p)
    CALL nbc(3,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    phi_p(1,2)=q2in; as(2,3)=0.
! Q3
  ELSEIF(ASSOCIATED(phi_p,u3))THEN
    CALL dbc(1,0.,psi_p=phi_p); CALL nbc(2,psi_p=phi_p)
    CALL nbc(3,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    phi_p(1,2)=uin; as(2,3)=0.
  ELSEIF(ASSOCIATED(phi_p,v3))THEN
    CALL dbc(1,0.,psi_p=phi_p); CALL nbc(2,psi_p=phi_p)
    CALL dbc(3,0.,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    phi_p(1,2)=vin; as(2,3)=0.
  ELSEIF(ASSOCIATED(phi_p,q3))THEN
    CALL dbc(1,0.,psi_p=phi_p); CALL nbc(2,psi_p=phi_p)
    CALL nbc(3,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    phi_p(1,2)=q3in; as(2,3)=0.
! REBOUNDING SPRAY
! Q0i
  ELSEIF(ASSOCIATED(phi_p,u0i))THEN
    CALL nbc(1,psi_p=phi_p); CALL dbc(2,0.,psi_p=phi_p)
    CALL nbc(3,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    DO j=2,ny-1
      phi_p(nx,j)=-u0(nx-1,j)*0.3
    ENDDO
  ELSEIF(ASSOCIATED(phi_p,v0i))THEN
    CALL nbc(1,psi_p=phi_p); CALL dbc(2,0.,psi_p=phi_p)
    CALL dbc(3,0.,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    DO j=2,ny-1
      phi_p(nx,j)=v0(nx-1,j)*1.4
    ENDDO
  ELSEIF(ASSOCIATED(phi_p,q0i))THEN
    CALL nbc(1,psi_p=phi_p); CALL dbc(2,0.,psi_p=phi_p)
    CALL nbc(3,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    DO j=2,ny-1
      phi_p(nx,j)=q0(nx-1,j)! wall(nx-1,j,0,lreb=t)
    ENDDO
! Q1i
  ELSEIF(ASSOCIATED(phi_p,u1i))THEN
    CALL nbc(1,psi_p=phi_p); CALL dbc(2,0.,psi_p=phi_p)
    CALL nbc(3,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    DO j=2,ny-1
      phi_p(nx,j)=-u1(nx-1,j)*0.3
    ENDDO
  ELSEIF(ASSOCIATED(phi_p,v1i))THEN
    CALL nbc(1,psi_p=phi_p); CALL dbc(2,0.,psi_p=phi_p)
    CALL dbc(3,0.,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    DO j=2,ny-1
      phi_p(nx,j)=v1(nx-1,j)*1.4
    ENDDO
  ELSEIF(ASSOCIATED(phi_p,q1i))THEN
    CALL nbc(1,psi_p=phi_p); CALL dbc(2,0.,psi_p=phi_p)
    CALL nbc(3,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    DO j=2,ny-1
      phi_p(nx,j)=q1(nx-1,j)! wall(nx-1,j,1,lreb=t)
    ENDDO
! Q2i
  ELSEIF(ASSOCIATED(phi_p,u2i))THEN
    CALL nbc(1,psi_p=phi_p); CALL dbc(2,0.,psi_p=phi_p)
    CALL nbc(3,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    DO j=2,ny-1
      phi_p(nx,j)=-u2(nx-1,j)*0.3
    ENDDO
  ELSEIF(ASSOCIATED(phi_p,v2i))THEN
    CALL nbc(1,psi_p=phi_p); CALL dbc(2,0.,psi_p=phi_p)
    CALL dbc(3,0.,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    DO j=2,ny-1
      phi_p(nx,j)=v2(nx-1,j)*1.4
    ENDDO
  ELSEIF(ASSOCIATED(phi_p,q2i))THEN
    CALL nbc(1,psi_p=phi_p); CALL dbc(2,0.,psi_p=phi_p)
    CALL nbc(3,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    DO j=2,ny-1
      phi_p(nx,j)=q2(nx-1,j)! wall(nx-1,j,2,lreb=t)
    ENDDO
! Q3i
  ELSEIF(ASSOCIATED(phi_p,u3i))THEN
    CALL nbc(1,psi_p=phi_p); CALL dbc(2,0.,psi_p=phi_p)
    CALL nbc(3,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    DO j=2,ny-1
      phi_p(nx,j)=-u3(nx-1,j)*0.3
    ENDDO
  ELSEIF(ASSOCIATED(phi_p,v3i))THEN
    CALL nbc(1,psi_p=phi_p); CALL dbc(2,0.,psi_p=phi_p)
    CALL dbc(3,0.,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    DO j=2,ny-1
      phi_p(nx,j)=v3(nx-1,j)*1.4
    ENDDO
  ELSEIF(ASSOCIATED(phi_p,q3i))THEN
    CALL nbc(1,psi_p=phi_p); CALL dbc(2,0.,psi_p=phi_p)
    CALL nbc(3,psi_p=phi_p); CALL nbc(4,psi_p=phi_p)
    DO j=2,ny-1
      phi_p(nx,j)=q3(nx-1,j)! wall(nx-1,j,3,lreb=t)
    ENDDO
  ENDIF
END SUBROUTINE boundary
!------------------------------------------------
! ADDITION OF SOURCE TERMS TO TRANSPORT EQU.
!------------------------------------------------
!
SUBROUTINE source(phi_p)
  USE mod_all
  REAL,DIMENSION(:,:),POINTER::phi_p
  DO i=2,nx-1
    DO j=2,ny-1
! GAS PHASE
      IF(ASSOCIATED(phi_p,u))THEN
        su(i,j)=su(i,j)+du3(i,j)*vol(i,j)+du3i(i,j)*vol(i,j)-(1.-(4/3.)*pi*q3s(i,j))*dpx(i,j)*vol(i,j)
      ELSEIF(ASSOCIATED(phi_p,v))THEN
        su(i,j)=su(i,j)+dv3(i,j)*vol(i,j)+dv3i(i,j)*vol(i,j)-(1.-(4/3.)*pi*q3s(i,j))*dpy(i,j)*vol(i,j)
! INCOMING SPRAY
! Q0
      ELSEIF(ASSOCIATED(phi_p,u0))THEN
        su(i,j)=su(i,j)-du0(i,j)*vol(i,j)+bq0(i,j)*u3(i,j)*vol(i,j)+cq0(i,j)*u0(i,j)*vol(i,j)
      ELSEIF(ASSOCIATED(phi_p,v0))THEN
        su(i,j)=su(i,j)-dv0(i,j)*vol(i,j)+bq0(i,j)*v3(i,j)*vol(i,j)+cq0(i,j)*v0(i,j)*vol(i,j)
      ELSEIF(ASSOCIATED(phi_p,q0))THEN
        su(i,j)=su(i,j)+bq0(i,j)*vol(i,j)+cq0(i,j)*vol(i,j)
! Q1
      ELSEIF(ASSOCIATED(phi_p,u1))THEN
        su(i,j)=su(i,j)-du1(i,j)*vol(i,j)+bq1(i,j)*u3(i,j)*vol(i,j)+cq1(i,j)*u1(i,j)*vol(i,j)
      ELSEIF(ASSOCIATED(phi_p,v1))THEN
        su(i,j)=su(i,j)-dv1(i,j)*vol(i,j)+bq1(i,j)*v3(i,j)*vol(i,j)+cq1(i,j)*v1(i,j)*vol(i,j)
      ELSEIF(ASSOCIATED(phi_p,q1))THEN
        su(i,j)=su(i,j)+bq1(i,j)*vol(i,j)+cq1(i,j)*vol(i,j)
! Q2
      ELSEIF(ASSOCIATED(phi_p,u2))THEN
        su(i,j)=su(i,j)-du2(i,j)*vol(i,j)+bq2(i,j)*u3(i,j)*vol(i,j)+cq2(i,j)*u2(i,j)*vol(i,j)
      ELSEIF(ASSOCIATED(phi_p,v2))THEN
        su(i,j)=su(i,j)-dv2(i,j)*vol(i,j)+bq2(i,j)*v3(i,j)*vol(i,j)+cq2(i,j)*v2(i,j)*vol(i,j)
      ELSEIF(ASSOCIATED(phi_p,q2))THEN
        su(i,j)=su(i,j)+bq2(i,j)*vol(i,j)+cq2(i,j)*vol(i,j)
! Q3
      ELSEIF(ASSOCIATED(phi_p,u3))THEN
        su(i,j)=su(i,j)-du3(i,j)*vol(i,j)
      ELSEIF(ASSOCIATED(phi_p,v3))THEN
        su(i,j)=su(i,j)-dv3(i,j)*vol(i,j)
      ELSEIF(ASSOCIATED(phi_p,q3))THEN
        su(i,j)=su(i,j)
! REBOUNDING SPRAY
! Q0i
      ELSEIF(ASSOCIATED(phi_p,u0i))THEN
        su(i,j)=su(i,j)-du0i(i,j)*vol(i,j)+bq0i(i,j)*u3i(i,j)*vol(i,j)+cq0i(i,j)*u0i(i,j)*vol(i,j)
      ELSEIF(ASSOCIATED(phi_p,v0i))THEN
        su(i,j)=su(i,j)-dv0i(i,j)*vol(i,j)+bq0i(i,j)*v3i(i,j)*vol(i,j)+cq0i(i,j)*v0i(i,j)*vol(i,j)
      ELSEIF(ASSOCIATED(phi_p,q0i))THEN
        su(i,j)=su(i,j)+bq0i(i,j)*vol(i,j)+cq0i(i,j)*vol(i,j)
! Q1i
      ELSEIF(ASSOCIATED(phi_p,u1i))THEN
        su(i,j)=su(i,j)-du1i(i,j)*vol(i,j)+bq1i(i,j)*u3i(i,j)*vol(i,j)+cq1i(i,j)*u1i(i,j)*vol(i,j)
      ELSEIF(ASSOCIATED(phi_p,v1i))THEN
        su(i,j)=su(i,j)-dv1i(i,j)*vol(i,j)+bq1i(i,j)*v3i(i,j)*vol(i,j)+cq1i(i,j)*v1i(i,j)*vol(i,j)
      ELSEIF(ASSOCIATED(phi_p,q1i))THEN
        su(i,j)=su(i,j)+bq1i(i,j)*vol(i,j)+cq1i(i,j)*vol(i,j)
! Q2i
      ELSEIF(ASSOCIATED(phi_p,u2i))THEN
        su(i,j)=su(i,j)-du2i(i,j)*vol(i,j)+bq2i(i,j)*u3i(i,j)*vol(i,j)+cq2i(i,j)*u2i(i,j)*vol(i,j)
      ELSEIF(ASSOCIATED(phi_p,v2i))THEN
        su(i,j)=su(i,j)-dv2i(i,j)*vol(i,j)+bq2i(i,j)*v3i(i,j)*vol(i,j)+cq2i(i,j)*v2i(i,j)*vol(i,j)
      ELSEIF(ASSOCIATED(phi_p,q2i))THEN
        su(i,j)=su(i,j)+bq2i(i,j)*vol(i,j)+cq2i(i,j)*vol(i,j)
! Q3i
      ELSEIF(ASSOCIATED(phi_p,u3i))THEN
        su(i,j)=su(i,j)-du3i(i,j)*vol(i,j)
      ELSEIF(ASSOCIATED(phi_p,v3i))THEN
        su(i,j)=su(i,j)-dv3i(i,j)*vol(i,j)
      ELSEIF(ASSOCIATED(phi_p,q3i))THEN
        su(i,j)=su(i,j)
      ENDIF
    ENDDO
  ENDDO
END SUBROUTINE source
!------------------------------------------------
! SOLUTION OF ALL SOURCE TERMS
!------------------------------------------------
!
SUBROUTINE sources
  USE mod_all
  USE mod_int_drag
  USE mod_int_breakup
  USE mod_int_collision
  USE mod_int_moment
  REAL,DIMENSION(nnx,nny)::qm2,qm1
  DO i=2,nx-1
    DO j=2,ny-1
! INCOMING SPRAY
      IF(q0(i,j)>q0in*order.and.q1(i,j)>q1in*order &
      .and.q2(i,j)>q2in*order.and.q3(i,j)>q3in*order)THEN
! MOMENTS
        rlb=0.; rub=0.
        CALL moment(i,j,-2.,rlb,rub,q0,q1,q2,q3,q); qm2(i,j)=q
        CALL moment(i,j,-1.,rlb,rub,q0,q1,q2,q3,q); qm1(i,j)=q
! DRAG
        du0(i,j)=drag(i,j,u,u0,qm2,qm1)
        dv0(i,j)=drag(i,j,v,v0,qm2,qm1)
        du1(i,j)=drag(i,j,u,u1,qm1,q0)
        dv1(i,j)=drag(i,j,v,v1,qm1,q0)
        du2(i,j)=drag(i,j,u,u2,q0,q1)
        dv2(i,j)=drag(i,j,v,v2,q0,q1)
        du3(i,j)=drag(i,j,u,u3,q1,q2)*(4/3.)*pi*den_l(i,j)
        dv3(i,j)=drag(i,j,v,v3,q1,q2)*(4/3.)*pi*den_l(i,j)
! BREAK-UP
        bq0(i,j)=breakup(i,j,0,rub,u0,v0,q0,q1,q2,q3)
        bq1(i,j)=breakup(i,j,1,rub,u1,v1,q0,q1,q2,q3)
        bq2(i,j)=breakup(i,j,2,rub,u2,v2,q0,q1,q2,q3)
! COLLISION
        cq0(i,j)=collision(i,j,0,rub,u3,v3,u0,v0,q0,q1,q2,q3)
        cq1(i,j)=collision(i,j,1,rub,u3,v3,u0,v0,q0,q1,q2,q3)
        cq2(i,j)=collision(i,j,2,rub,u3,v3,u0,v0,q0,q1,q2,q3)
      ELSE
        q0(i,j)=0.; q1(i,j)=0.; q2(i,j)=0.; q3(i,j)=0.
        du0(i,j)=0.; du1(i,j)=0.; du2(i,j)=0.; du3(i,j)=0.
        dv0(i,j)=0.; dv1(i,j)=0.; dv2(i,j)=0.; dv3(i,j)=0.
        bq0(i,j)=0.; bq1(i,j)=0.; bq2(i,j)=0.
        cq0(i,j)=0.; cq1(i,j)=0.; cq2(i,j)=0.
      ENDIF
! REBOUNDING SPRAY
      IF(q0i(i,j)>q0in*order.and.q1i(i,j)>q1in*order &
      .and.q2i(i,j)>q2in*order.and.q3i(i,j)>q3in*order)THEN
! MOMENTS
        rlb=0.; rub=0.
        CALL moment(i,j,-2.,rlb,rub,q0i,q1i,q2i,q3i,q); qm2(i,j)=q
        CALL moment(i,j,-1.,rlb,rub,q0i,q1i,q2i,q3i,q); qm1(i,j)=q
! DRAG
        du0i(i,j)=drag(i,j,u,u0i,qm2,qm1)
        dv0i(i,j)=drag(i,j,v,v0i,qm2,qm1)
        du1i(i,j)=drag(i,j,u,u1i,qm1,q0i)
        dv1i(i,j)=drag(i,j,v,v1i,qm1,q0i)
        du2i(i,j)=drag(i,j,u,u2i,q0i,q1i)
        dv2i(i,j)=drag(i,j,v,v2i,q0i,q1i)
        du3i(i,j)=drag(i,j,u,u3i,q1i,q2i)*(4/3.)*pi*den_l(i,j)
        dv3i(i,j)=drag(i,j,v,v3i,q1i,q2i)*(4/3.)*pi*den_l(i,j)
! BREAK-UP
        bq0i(i,j)=breakup(i,j,0,rub,u0i,v0i,q0i,q1i,q2i,q3i)
        bq1i(i,j)=breakup(i,j,1,rub,u1i,v1i,q0i,q1i,q2i,q3i)
        bq2i(i,j)=breakup(i,j,2,rub,u2i,v2i,q0i,q1i,q2i,q3i)
! COLLISION
        cq0i(i,j)=collision(i,j,0,rub,u3i,v3i,u0i,v0i,q0i,q1i,q2i,q3i)
        cq1i(i,j)=collision(i,j,1,rub,u3i,v3i,u0i,v0i,q0i,q1i,q2i,q3i)
        cq2i(i,j)=collision(i,j,2,rub,u3i,v3i,u0i,v0i,q0i,q1i,q2i,q3i)
      ELSE
        q0i(i,j)=0.; q1i(i,j)=0.; q2i(i,j)=0.; q3i(i,j)=0.
        du0i(i,j)=0.; du1i(i,j)=0.; du2i(i,j)=0.; du3i(i,j)=0.
        dv0i(i,j)=0.; dv1i(i,j)=0.; dv2i(i,j)=0.; dv3i(i,j)=0.
        bq0i(i,j)=0.; bq1i(i,j)=0.; bq2i(i,j)=0.
        cq0i(i,j)=0.; cq1i(i,j)=0.; cq2i(i,j)=0.
      ENDIF
    ENDDO
  ENDDO
END SUBROUTINE sources
!------------------------------------------------
! DRAG FUNCTION
!------------------------------------------------
!
FUNCTION drag(i,j,vel_g,vel_l,qa,qb)
  USE mod_all
  INTEGER,INTENT(IN)::i,j
  REAL,DIMENSION(:,:),INTENT(IN)::vel_g,vel_l,qa,qb
  drag=0.
  IF(ldrag)THEN
    urel=vel_l(i,j)-vel_g(i,j)
    termu=u3(i,j)-u(i,j); termv=v3(i,j)-v(i,j)
    urelh=sqrt(termu**2+termv**2)
    IF(qa(i,j)>1.e-30)terma=4.5*vis_g(i,j)*qa(i,j)
    IF(qb(i,j)>1.e-30)termb=1.35*(den_g(i,j)*urelh*qb(i,j))**0.687 &
    *((vis_g(i,j)*qa(i,j))/2.)**0.313
    drag=(terma+termb)*urel/den_l(i,j)
  ENDIF
END FUNCTION drag
!------------------------------------------------
! BREAK-UP FUNCTION
!------------------------------------------------
!
FUNCTION breakup(i,j,iq,rub,ud,vd,p0,p1,p2,p3)
  USE mod_all,ONLY:lbreakup,q2in,q3in,u,v,den_g,den_l,vis_l,st
  USE mod_int_moment
  REAL,DIMENSION(:,:),INTENT(IN)::ud,vd,p0,p1,p2,p3
  INTEGER,INTENT(IN)::i,j,iq
  REAL,INTENT(INOUT)::rub
  breakup=0.
  IF(lbreakup)THEN
    velu=ud(i,j)-u(i,j)
    velv=vd(i,j)-v(i,j)
    urelh=SQRT(velu**2+velv**2)
! INTEGRAL LIMITS
    IF(urelh>1.e-2)THEN
      rlb=min(6.*st(i,j)/(den_g(i,j)*urelh**2),rub)
! TRUNCATED MOMENTS
      CALL moment(i,j,2.,rlb,rub,p0,p1,p2,p3,q); q2tr=q
      CALL moment(i,j,3.,rlb,rub,p0,p1,p2,p3,q); q3tr=q
      IF(q2tr>q2in*1.e-10.and.q3tr>q3in*1.e-10)THEN
        r32tr=q3tr/q2tr
! STABLE RADIUS
        rstab=6.2*SQRT(r32tr) &
        *(den_l(i,j)/den_g(i,j))**0.25 &
        *SQRT(vis_l(i,j)/(2.*den_l(i,j)*urelh))
! CHARACTERISTIC TIME SCALE
        tscale=2*5.*r32tr*SQRT(den_l(i,j)/den_g(i,j)) &
        /(urelh*(1.-oh(r32tr,i,j)/7.))
! MOMENT SOURCE TERM
        breakup=(r32tr**3/rstab**(3-iq)-r32tr**(iq))/tscale
      ENDIF
    ENDIF
  ENDIF
  END FUNCTION breakup
!------------------------------------------------
! COLLISION FUNCTION
!------------------------------------------------
!
FUNCTION collision(i,j,iq,rub,ud3,vd3,ud0,vd0,p0,p1,p2,p3)
  USE mod_all,ONLY:lcollision,pi,den_g,st
  USE mod_int_moment
  REAL,DIMENSION(:,:),INTENT(IN)::ud3,vd3,ud0,vd0,p0,p1,p2,p3
  INTEGER,INTENT(IN)::i,j,iq
  REAL,INTENT(INOUT)::rub
  collision=0.
  IF(lcollision)THEN
! RELATIVE VELOCITY COEFFICIENTS
    velu=ud3(i,j)-ud0(i,j)
    velv=vd3(i,j)-vd0(i,j)
    urelh=sqrt(velu**2+velv**2)
! NUMBER OF COLLISIONS PER UNIT VOLUME PER UNIT TIME (*DENSITY)
    collno=0.15*pi*urelh*(p0(i,j)*p2(i,j)+p1(i,j)**2)
! CRITICAL RADII GIVEN BY CRITICAL WEBER NUMBER
! ASSUMING THE CRITICAL WEBER NUMBERS ARE THOSE FOR DODECANE
! NAMELY WEA=2.5, WEB=10, WEC=30
    IF(collno>1.e-30.and.urelh>1.e-4)THEN
      coeff=st(i,j)/(2*den_g(i,j)*urelh**2)
      rlba=2.5*coeff; rlbb=10.*coeff; rlbc=30.*coeff
! AVERAGE RADII FOR EACH REGIME
      ra=rlba; rb=(rlbb+rlbc)/2; rc=rlbc; rd=rlbc
! PROBABILITY OF R<RCRITA
      CALL moment(i,j,0.,rlba,rub,p0,p1,p2,p3,q); proba=q/p0(i,j)
! PROBABILITY OF R<RCRITB
      CALL moment(i,j,0.,rlbb,rub,p0,p1,p2,p3,q); probb=q/p0(i,j)
! PROBABILITY OF R<RCRITC
      CALL moment(i,j,0.,rlbc,rub,p0,p1,p2,p3,q); probc=q/p0(i,j)
! PROBABILITY OF COALESCENCE IN REGION A
      pcoala=0.6*proba*(2-proba)
! PROBABILITY OF COALESCENCE IN REGION B
      pcoalb=0.5*(probc-probb)*(2-probb-probc)
! PROBABILITY OF COALESCENCE IN REGION C
      pcoalc=0.2*(1-probc)**2
! PROBABILITY OF SEPARATION IN REGION C
      psep=0.8*(1-probc)**2
! SOURCE
      dqa=0.;dqb=0.;dqc=0.;dqd=0.
      IF(ra>1.e-6.and.ra<rub)dqa=(2.**(REAL(iq)/3.)-2)*ra**(iq)
      IF(rb>1.e-6.and.rb<rub)dqb=(2.**(REAL(iq)/3.)-2)*rb**(iq)
      IF(rc>1.e-6.and.rc<rub)dqc=(2.**(REAL(iq)/3.)-2)*rc**(iq)
      IF(rd>1.e-6.and.rd<rub)dqd=(5*0.4**(iq/3.)-2)*rd**(iq)
      collision=collno*(pcoala*dqa+pcoalb*dqb+pcoalc*dqc+psep*dqd)
    ENDIF
  ENDIF
END FUNCTION collision
!------------------------------------------------
! WALL MOMENT FUNCTION
!------------------------------------------------
!
FUNCTION wall(i,j,iq,lreb,lspl)
  USE mod_all
  USE mod_int_moment
  REAL,DIMENSION(15)::qs,wes
  INTEGER,INTENT(IN)::i,j,iq
  LOGICAL,INTENT(IN),OPTIONAL::lreb,lspl
  IF(q0(i,j)>q0in*order.and.q1(i,j)>q1in*order &
  .and.q2(i,j)>q2in*order.and.q3(i,j)>q3in*order)THEN
    velu=u3(i,j)-u(i,j)
    velv=v3(i,j)-v(i,j)
    urelh=sqrt(velu**2+velv**2)
    rlb=0.;rub=0.
    CALL moment(i,j,3.,rlb,rub,q0,q1,q2,q3,q)
    ratio=q/q3(i,j)
    IF(ratio>0.5.and.ratio<2.)THEN
      delr=(rub-rlb)/15
      DO ns=1,15
        IF(ns==1)rrlb=rlb ;rub=rrlb+delr
        IF(ns/=1)rlb=rub ;rub=rub+delr
        CALL moment(i,j,REAL(iq),rlb,rub,q0,q1,q2,q3,q)
        qs(ns)=q; wes(ns)=we(i,j,(rlb+rub)/2,urelh)
        IF(lreb.and.wes(ns)>10.)qsi=qsi+qs(ns)
        IF(lspl.and.wes(ns)>5..and.wes(ns)<10.)qsi=qsi+qs(ns)
      ENDDO
      wall=qsi
    ENDIF
  ENDIF
END FUNCTION wall
!------------------------------------------------
! OHNEZORGE NUMBER
!------------------------------------------------
!
FUNCTION oh(rad,i,j)
  USE mod_all,ONLY:den_l,vis_l,st
  oh=vis_l(i,j)/SQRT(den_l(i,j)*rad*st(i,j))
END FUNCTION oh
!------------------------------------------------
! WEBER NUMBER
!------------------------------------------------
!
FUNCTION we(i,j,rad,urelh)
  USE mod_all,ONLY:den_g,st
  we=2*den_g(i,j)*urelh**2*rad/st(i,j)
END FUNCTION we
!------------------------------------------------
! OUTPUT: TECPLOT
!------------------------------------------------
!
SUBROUTINE output
  USE mod_all
15 FORMAT('variables="X mm","Y mm","U","V","P","U3","V3","Q3","R32","U3i","V3i","Q3i","R32i"')
25 FORMAT(13(e15.5))
35 FORMAT('zone f=point, i=',i4,',j=',i4)
  OPEN(10,file='TECPLOT.DAT')
  WRITE(10,15)
  WRITE(10,35)nx,ny
  DO j=1,ny
    DO i=1,nx
      r32=0.; IF(q2(i,j)>1.e-30)r32=q3(i,j)/q2(i,j)
      r32i=0.; IF(q2i(i,j)>1.e-30)r32i=q3i(i,j)/q2i(i,j)
      WRITE(10,25)xc(i)*1000,yc(j)*1000, &
      u(i,j),v(i,j),p(i,j), &
      u3(i,j),v3(i,j),q3(i,j),r32, &
      u3i(i,j),v3i(i,j),q3i(i,j),r32i
    ENDDO
  ENDDO
  CLOSE(10)
END SUBROUTINE output
!------------------------------------------------
! CALCULATION OF MOMENTS: LAGUERRE POLNOMIALS
!------------------------------------------------
! P0 -> Q0 ETC.
!
SUBROUTINE moment(i,j,power,rlb,rub,p0,p1,p2,p3,q)
  COMMON/lcoeff/b0,b1,b2,b3
  INTEGER,INTENT(IN)::i,j
  REAL,INTENT(IN)::power
  REAL,DIMENSION(:,:),INTENT(IN)::p0,p1,p2,p3
  REAL,INTENT(INOUT)::rlb,rub
  REAL,INTENT(OUT)::q
  REAL,EXTERNAL::flagu
! SIZE DISTRIBUTION NORMALISATION RADIUS
  rnorm=(p3(i,j)/p0(i,j))**(1/3.)
! NORMALISED MOMENTS
  a0=(p0(i,j)/p0(i,j))*(1./rnorm)**0; a1=(p1(i,j)/p0(i,j))*(1./rnorm)**1
  a2=(p2(i,j)/p0(i,j))*(1./rnorm)**2; a3=(p3(i,j)/p0(i,j))*(1./rnorm)**3
! LAGUERRE COEFFICIENTS
  b0=plagu(0,a0); b1=plagu(1,a1); b2=plagu(2,a2); b3=plagu(3,a3)
! LIMITS OF DISTRIBUTION FUNCTION
  root1=0; root2=0; delr=0.5e-6
  IF(rlb<1.e-6)THEN
    rrlb=1.e-6
  ELSE
    root1=rlb; rrlb=rlb
  ENDIF
  IF(rub<1.e-6)THEN
    rrub=50.e-6
  ELSE
    root2=rub; rrub=rub
  ENDIF
  IF(abs(root1)<1.e-30)THEN
  pfunc=3.
  DO r=rrlb,rrub,delr
    terml=flagu((r-delr)/rnorm,pfunc)
    termc=flagu(r/rnorm,pfunc)
    termu=flagu((r+delr)/rnorm,pfunc)
    IF(termc>1.e-30.and.terml<termu)THEN
      root1=r; GOTO 1
    ENDIF
  ENDDO
  DO r=rrlb,rrub,delr
    terml=flagu((r-delr)/rnorm,pfunc)
    termc=flagu(r/rnorm,pfunc)
    termu=flagu((r+delr)/rnorm,pfunc)
    IF(terml>1.e-30.and.termc>1.e-30.and.termu>1.e-30)then
      root1=r; GOTO 1
    ENDIF
  ENDDO
1 CONTINUE
  ENDIF
  IF(abs(root2)<1.e-30)THEN
  DO r=rrub,rrlb,-delr
    terml=flagu((r-delr)/rnorm,pfunc)
    termc=flagu(r/rnorm,pfunc)
    termu=flagu((r+delr)/rnorm,pfunc)
    IF(termc>1.e-30.and.terml>termu.and.r>root1)THEN
      root2=r; GOTO 2
    ENDIF
  ENDDO
2 CONTINUE
  ENDIF
! ADJUST UPPER LIMIT
  IF(root1>1.e-30.and.root2>root1)THEN
    rlbn=root1/rnorm
    rubn=root2/rnorm
    icount=0; icountmx=10
  3 CONTINUE
    icount=icount+1
    pfunc=3.
    CALL qalt(flagu,pfunc,rlbn,rubn,qan)
    qnorm=p0(i,j)/(1./rnorm)**(pfunc)
    qa=qan*qnorm
    ratio=qa/p3(i,j)
    IF(ratio<0.98 .and.icount<icountmx)THEN
      IF(ratio>0.8)rubn=rubn+0.02
      IF(ratio<0.8)rubn=rubn+0.2
      GOTO 3
    ELSEIF(ratio>1.02 .and.icount<icountmx)THEN
      IF(ratio<1.2)rubn=rubn-0.02
      IF(ratio>1.2)rubn=rubn-0.2
      GOTO 3
    ENDIF
    IF(rlbn>1.e-30.and.rubn>rlbn)THEN
! INTEGRATION
      CALL qalt(flagu,power,rlbn,rubn,qn)
! UN-NORMALISE MOMENT
      qnorm=p0(i,j)/(1./rnorm)**(power)
      q=qn*qnorm; rlb=root1; rub=root2
    ENDIF
  ENDIF
  IF(q<1.e-30)q=0.
END SUBROUTINE moment
!------------------------------------------------
! LAGUERRE SIZE DISTRIBUTION FUNCTION
!------------------------------------------------
!
FUNCTION flagu(x,power)
  COMMON/lcoeff/b0,b1,b2,b3
  REAL,INTENT(IN)::x,power
  flagu=exp(-x)*x**(power) &
  *(b0*plagu(0,x)+b1*plagu(1,x)+b2*plagu(2,x)+b3*plagu(3,x))
END FUNCTION flagu
!------------------------------------------------
! LAGUERRE POLYNOMIALS
!------------------------------------------------
!
FUNCTION plagu(iq,var)
  INTEGER,INTENT(IN)::iq
  REAL,INTENT(IN)::var
  IF(iq==0)THEN
    plagu=1.
  ELSEIF(iq==1)THEN
    plagu=1.-var
  ELSEIF(iq==2)THEN
    plagu=(1/2.)*(2.-4.*var+var**2)
  ELSEIF(iq==3)THEN
    plagu=(1/6.)*(6.-18.*var+9.*var**2-var**3)
  ELSEIF(iq==4)THEN
    plagu=(1/24.)*(24.-96.*var+72.*var**2-16.*var**3+var**4)
  ENDIF
END FUNCTION plagu
!------------------------------------------------
! INTEGRAL OF LAGUERRE FUNCTION
!------------------------------------------------
!
SUBROUTINE qalt(func,power,a,b,s)
  REAL,INTENT(IN)::power,a,b
  REAL,INTENT(INOUT)::s
  REAL,EXTERNAL::func
  CALL qtrap(func,power,a,b,s)
END SUBROUTINE qalt
!------------------------------------------------
! INTEGRATION REFINEMENT [Num. Res. in F77, Ch 4.2]
!------------------------------------------------
!
SUBROUTINE qtrap(func,power,a,b,s)
  INTEGER::j
  REAL::power,a,b,s,olds
  REAL,PARAMETER::eps=1.e-6
  INTEGER,PARAMETER::jmax=5
  REAL,EXTERNAL::func
  olds=-1.e+30
  DO j=1,jmax
    CALL trapzd(func,power,a,b,s,j)
    olds=s
    IF(j>5.and.abs(s-olds)<eps*abs(olds))EXIT
  ENDDO
END SUBROUTINE qtrap
!------------------------------------------------
! NUMERICAL INTEGRATION
! USING EXTENDED TRAPEZOIDAL RULE [Num. Res. in F77, Ch 4.2]
!------------------------------------------------
!
SUBROUTINE trapzd(func,power,a,b,s,n)
  INTEGER::n,it,j
  REAL::power,a,b,s,del,sum,tnm,x
  REAL,EXTERNAL::func
  IF(n==1)THEN
    s=0.5*(b-a)*(func(a,power)+func(b,power))
  ELSE
    it=2**(n-2)
    tnm=it
    del=(b-a)/tnm
    x=a+0.5*del
    sum=0.
    DO j=1,it
      sum=sum+func(x,power)
      x=x+del
    ENDDO
    s=0.5*(s+(b-a)*sum/tnm)
  ENDIF
END SUBROUTINE trapzd