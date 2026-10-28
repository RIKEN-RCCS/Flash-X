program test_update_wrapper
  implicit none
#include "Simulation.h"
#include "constants.h"
#include "Spark.h"
  interface
subroutine hy_rk_updateSoln(stage,starState,tmpState,grav,flx,fly,flz, &
    deltas,fareaX,fareaY,fareaZ,cvol,xCenter,xLeft,xRight,yLeft,yRight, &
    dt,dtOld,limits,lo,loGC)
  integer, intent(in) :: stage,lo(3),loGC(3),limits(LOW:HIGH,MDIM,MAXSTAGE)
  real, intent(inout) :: starState(1:,loGC(1):,loGC(2):,loGC(3):)
  real, intent(in) :: tmpState(1:,loGC(1):,loGC(2):,loGC(3):)
  real, intent(in) :: grav(1:,loGC(1):,loGC(2):,loGC(3):)
  real, intent(in) :: flx(1:,loGC(1):,loGC(2):,loGC(3):)
  real, intent(in) :: fly(1:,loGC(1):,loGC(2):,loGC(3):)
  real, intent(in) :: flz(1:,loGC(1):,loGC(2):,loGC(3):)
  real, intent(in) :: deltas(MDIM),dt,dtOld
  real, intent(in) :: fareaX(loGC(1):,loGC(2):,loGC(3):)
  real, intent(in) :: fareaY(loGC(1):,loGC(2):,loGC(3):),fareaZ(loGC(1):,loGC(2):,loGC(3):)
  real, intent(in) :: cvol(loGC(1):,loGC(2):,loGC(3):),xCenter(loGC(1):)
  real, intent(in) :: xLeft(loGC(1):),xRight(loGC(1):),yLeft(loGC(2):),yRight(loGC(2):)
end subroutine
  end interface
  integer :: limits(2,3,MAXSTAGE),lo(3),loGC(3),d,v,i,j,k
  real :: state(11,0:3,0:3,0:3),ref(11,0:3,0:3,0:3),gravity(3,0:3,0:3,0:3)
  real :: fx(5,0:3,0:3,0:3),fy(5,0:3,0:3,0:3),fz(5,0:3,0:3,0:3)
  real :: area(0:3,0:3,0:3),coords(0:3),spacing(3),us(5),u0(5),source(5),div(5),expected(5)
  state=0; ref=0; area=1; coords=0; lo=1; loGC=0; spacing=[0.5,0.25,2.]
  limits=0; limits(LOW,:,MAXSTAGE)=1; limits(HIGH,:,MAXSTAGE)=2
  state(DENS_VAR,:,:,:)=3.; state(VELX_VAR,:,:,:)=1.; state(VELY_VAR,:,:,:)=-0.5
  state(VELZ_VAR,:,:,:)=0.25; state(ENER_VAR,:,:,:)=10.; state(EINT_VAR,:,:,:)=2000.
  state(PRES_VAR,:,:,:)=99.; state(GPOT_VAR,:,:,:)=3.; state(GPOL_VAR,:,:,:)=9.
  ref(DENS_VAR,:,:,:)=2.; ref(ENER_VAR,:,:,:)=8.
  ref(GPOT_VAR,:,:,:)=4.; ref(GPOL_VAR,:,:,:)=2.
  gravity(1,:,:,:)=1.; gravity(2,:,:,:)=-2.; gravity(3,:,:,:)=3.
  do k=0,3
    do j=0,3
      do i=0,3
        do v=1,5
          fx(v,i,j,k)=0.01*v*i; fy(v,i,j,k)=0.02*v*j; fz(v,i,j,k)=0.03*v*k
        end do
      end do
    end do
  end do
  us=[3.,3.,-1.5,0.75,30.]; u0=[2.,0.,0.,0.,16.]; source=[0.,3.,-6.,9.,8.25]
  do v=1,5
    div(v)=0.01*v/spacing(1)
    if (NDIM>1) div(v)=div(v)+0.02*v/spacing(2)
    if (NDIM>2) div(v)=div(v)+0.03*v/spacing(3)
  end do
  expected=0.25*u0+0.75*us+0.1*(source-div)
  call hy_rk_updateSoln(MAXSTAGE,state,ref,gravity,fx,fy,fz,spacing,area,area,area,area,coords, &
                       coords,coords,coords,coords,0.2,2.,limits,lo,loGC)
  if (abs(state(DENS_VAR,1,1,1)-expected(1))>1.e-5) stop 1
  if (abs(state(VELX_VAR,1,1,1)-expected(2)/expected(1))>1.e-5) stop 2
  if (abs(state(ENER_VAR,1,1,1)-expected(5)/expected(1))>1.e-5) stop 3
  if (abs(state(GPOT_VAR,1,1,1)-3.35)>1.e-5) stop 4
  if (state(GPOL_VAR,1,1,1)/=9. .or. state(PRES_VAR,1,1,1)/=99.) stop 5
  if (state(DENS_VAR,0,0,0)/=3.) stop 6
  print *, 'update wrapper tests passed'
end program
