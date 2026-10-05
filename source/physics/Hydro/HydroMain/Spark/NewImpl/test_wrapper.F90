program test_wrapper
  implicit none
#include "Simulation.h"
#include "constants.h"
#include "Spark.h"
  interface
    subroutine hy_rk_getFaceFlux(stage,starState,flat3d,flx,fly,flz,limits,deltas, &
                                scr_rope,scr_flux,scr_uPlus,scr_uMinus,loGC)
      integer, intent(in) :: stage,loGC(3),limits(LOW:HIGH,MDIM,MAXSTAGE)
      real, intent(in) :: starState(1:,loGC(1):,loGC(2):,loGC(3):)
      real, intent(in) :: flat3d(loGC(1):,loGC(2):,loGC(3):),deltas(MDIM)
      real, intent(out) :: flx(1:,loGC(1):,loGC(2):,loGC(3):)
      real, intent(out) :: fly(1:,loGC(1):,loGC(2):,loGC(3):)
      real, intent(out) :: flz(1:,loGC(1):,loGC(2):,loGC(3):)
      real, intent(out) :: scr_rope(1:,loGC(1):,loGC(2):,loGC(3):)
      real, intent(out) :: scr_flux(1:,loGC(1):,loGC(2):,loGC(3):)
      real, intent(out) :: scr_uPlus(1:,loGC(1):,loGC(2):,loGC(3):)
      real, intent(out) :: scr_uMinus(1:,loGC(1):,loGC(2):,loGC(3):)
    end subroutine
  end interface
  integer :: limits(2,3,MAXSTAGE),loGC(3),stage,v,i,j,k,d
  real :: state(11,-3:5,-3:5,-3:5),flat(-3:5,-3:5,-3:5),deltas(3)
  real :: fx(5,-3:5,-3:5,-3:5),fy(5,-3:5,-3:5,-3:5),fz(5,-3:5,-3:5,-3:5)
  real :: rhoL,rhoR,sl,sr,ul(5),ur(5),fl(5),fr(5),actual(5)
  real :: rope(7,-3:5,-3:5,-3:5),flux(5,-3:5,-3:5,-3:5)
  real :: plus(7,-3:5,-3:5,-3:5),minus(7,-3:5,-3:5,-3:5),q(7),u(5),expected(5),un
  q=[2.,0.4,-0.3,0.2,1.,1.4,5.]
  state=0; flat=0.3; deltas=0.1; loGC=-3
  state(DENS_VAR,:,:,:)=q(1); state(VELX_VAR,:,:,:)=q(2)
  state(VELY_VAR,:,:,:)=q(3); state(VELZ_VAR,:,:,:)=q(4)
  state(PRES_VAR,:,:,:)=q(5); state(GAMC_VAR,:,:,:)=q(6)
  state(EINT_VAR,:,:,:)=q(7)/q(1)
  limits=0; limits(LOW,:,MAXSTAGE)=1; limits(HIGH,:,MAXSTAGE)=2
  u=[q(1),q(1)*q(2),q(1)*q(3),q(1)*q(4),0.5*q(1)*sum(q(2:4)**2)+q(7)]
  call hy_rk_getFaceFlux(MAXSTAGE,state,flat,fx,fy,fz,limits,deltas,rope,flux,plus,minus,loGC)
  if (maxval(abs(rope(:,1,1,1)-q))>1.e-5) stop 1
  if (maxval(abs(plus(:,1,1,1)-q))>1.e-5) stop 2
  if (maxval(abs(minus(:,1,1,1)-q))>1.e-5) stop 3
  do d=1,NDIM
    un=q(1+d); expected=u*un
    expected(1+d)=expected(1+d)+q(5)
    expected(5)=(u(5)+q(5))*un
    select case(d)
    case(1)
      if (maxval(abs(fx(:,1,1,1)-expected))>1.e-5) stop 4
    case(2)
      if (maxval(abs(fy(:,1,1,1)-expected))>1.e-5) stop 5
    case(3)
      if (maxval(abs(fz(:,1,1,1)-expected))>1.e-5) stop 6
    end select
  end do
  if (maxval(abs(flux(:,1,1,1)-expected))>1.e-5) stop 7
  ! Nonuniform density with supplied zero flattening: constant reconstructions,
  ! marker-selected HLLE, and original-cell energy in every scratch/output.
  flat=0.
  do k=-3,5
    do j=-3,5
      do i=-3,5
        rhoR=2.+0.03*i+0.02*j+0.01*k
        state(DENS_VAR,i,j,k)=rhoR
        state(EINT_VAR,i,j,k)=5./rhoR
        state(SHOK_VAR,i,j,k)=1.
      end do
    end do
  end do
  call hy_rk_getFaceFlux(MAXSTAGE,state,flat,fx,fy,fz,limits,deltas,rope,flux,plus,minus,loGC)
  rhoR=state(DENS_VAR,1,1,1)
  if (abs(rope(1,1,1,1)-rhoR)>1.e-5 .or. abs(rope(7,1,1,1)-5.)>1.e-5) stop 8
  if (abs(plus(1,1,1,1)-rhoR)>1.e-5 .or. abs(minus(1,1,1,1)-rhoR)>1.e-5) stop 9
  do d=1,NDIM
    select case(d)
    case(1)
      rhoL=state(DENS_VAR,0,1,1)
    case(2)
      rhoL=state(DENS_VAR,1,0,1)
    case(3)
      rhoL=state(DENS_VAR,1,1,0)
    end select
    un=q(1+d)
    ul=[rhoL,rhoL*q(2),rhoL*q(3),rhoL*q(4),5.+0.5*rhoL*sum(q(2:4)**2)]
    ur=[rhoR,rhoR*q(2),rhoR*q(3),rhoR*q(4),5.+0.5*rhoR*sum(q(2:4)**2)]
    fl=ul*un; fr=ur*un
    fl(1+d)=fl(1+d)+1.; fr(1+d)=fr(1+d)+1.
    fl(5)=(ul(5)+1.)*un; fr(5)=(ur(5)+1.)*un
    sl=min(un-sqrt(1.4/rhoL),un-sqrt(1.4/rhoR))
    sr=max(un+sqrt(1.4/rhoL),un+sqrt(1.4/rhoR))
    expected=(sr*fl-sl*fr+sr*sl*(ur-ul))/(sr-sl)
    select case(d)
    case(1)
      actual=fx(:,1,1,1)
    case(2)
      actual=fy(:,1,1,1)
    case(3)
      actual=fz(:,1,1,1)
    end select
    if (maxval(abs(actual-expected))>1.e-5) stop 10
  end do
  if (maxval(abs(flux(:,1,1,1)-expected))>1.e-5) stop 11
  print *, 'wrapper tests passed' 
end program
