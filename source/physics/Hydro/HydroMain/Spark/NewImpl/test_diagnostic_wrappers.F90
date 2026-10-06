program test_diagnostic_wrappers
  use Hydro_data
  use Grid_tile
  implicit none
#include "Simulation.h"
#include "constants.h"
  interface
subroutine Hydro_computeDt(blockDesc,x,dx,uxgrid,y,dy,uygrid,z,dz,uzgrid, &
                           blkLimits,blkLimitsGC,U,dtCheck,dtMinLoc,extraInfo)
  use Grid_tile, only: Grid_tile_t
  type(Grid_tile_t), intent(in) :: blockDesc
  integer, intent(in) :: blkLimits(LOW:HIGH,MDIM),blkLimitsGC(LOW:HIGH,MDIM)
  real, intent(in) :: x(blkLimitsGC(LOW,1):blkLimitsGC(HIGH,1)),dx(blkLimitsGC(LOW,1):blkLimitsGC(HIGH,1)),uxgrid( &
  & blkLimitsGC(LOW,1):blkLimitsGC(HIGH,1))
  real, intent(in) :: y(blkLimitsGC(LOW,2):blkLimitsGC(HIGH,2)),dy(blkLimitsGC(LOW,2):blkLimitsGC(HIGH,2)),uygrid( &
  & blkLimitsGC(LOW,2):blkLimitsGC(HIGH,2))
  real, intent(in) :: z(blkLimitsGC(LOW,3):blkLimitsGC(HIGH,3)),dz(blkLimitsGC(LOW,3):blkLimitsGC(HIGH,3)),uzgrid( &
  & blkLimitsGC(LOW,3):blkLimitsGC(HIGH,3))
  real, pointer :: U(:,:,:,:)
  real, intent(inout) :: dtCheck
  integer, intent(inout) :: dtMinLoc(5)
  real, optional, intent(inout) :: extraInfo
end subroutine
subroutine hy_rk_shockDetect(Uin,Vc,blkLimitsGC,loGC)
  use Grid_tile, only: Grid_tile_t
  integer, intent(in) :: loGC(3),blkLimitsGC(LOW:HIGH,MDIM)
  real, intent(inout) :: Uin(1:,loGC(1):,loGC(2):,loGC(3):)
  real, intent(out) :: Vc(loGC(1):,loGC(2):,loGC(3):)
end subroutine
  end interface
  type(Grid_tile_t) :: tile
  real, target :: storage(11,-1:5,-1:5,-1:5)
  real, pointer :: state(:,:,:,:)
  real :: sound(-1:5,-1:5,-1:5),x(-1:5),y(-1:5),z(-1:5)
  real :: dx(-1:5),dy(-1:5),dz(-1:5),wx(-1:5),wy(-1:5),wz(-1:5),dt,extra,expected
  integer :: loGC(3),guard(2,3),bounds(2,3),loc(5),i,j,k
  state=>storage; loGC=-1; guard(LOW,:)=0; guard(HIGH,:)=4
  bounds(LOW,:)=1; bounds(HIGH,:)=3
  storage=0
  storage(DENS_VAR,:,:,:)=1.; storage(GAMC_VAR,:,:,:)=1.; storage(PRES_VAR,:,:,:)=1.
  storage(SHOK_VAR,:,:,:)=9.
  do k=-1,5
    do j=-1,5
      do i=-1,5
        storage(PRES_VAR,i,j,k)=1.+max(0,i)
        storage(VELX_VAR,i,j,k)=-0.5*i
      end do
    end do
  end do
  call hy_rk_shockDetect(state,sound,guard,loGC)
  if (storage(SHOK_VAR,1,1,1)/=1.) stop 1
  if (storage(SHOK_VAR,0,1,1)/=0.) stop 2
  if (storage(SHOK_VAR,-1,1,1)/=9.) stop 3
  if (abs(sound(1,1,1)-sqrt(2.))>1.e-5 .or. sound(-1,1,1)/=0.) stop 4
  storage(PRES_VAR,:,:,:)=1.; storage(VELX_VAR,:,:,:)=0.
  x=2.; y=acos(-1.)/2; z=0.; dx=0.5; dy=0.25; dz=1.; wx=0.; wy=-2.; wz=-3.
  tile%level=5; hy_meshMe=7; hy_lChyp=0.; dt=100.; loc=-1; extra=123.
  hy_useHydro=.false.
  call Hydro_computeDt(tile,x,dx,wx,y,dy,wy,z,dz,wz,bounds,guard,state,dt,loc,extra)
  if (dt/=100. .or. .not.hy_hydroComputeDtFirstCall .or. any(loc/=-1)) stop 5
  hy_useHydro=.true.
  call Hydro_computeDt(tile,x,dx,wx,y,dy,wy,z,dz,wz,bounds,guard,state,dt,loc,extra)
  expected=0.4
  if (NDIM>1) expected=0.8/12.
  if (abs(dt-expected)>1.e-5 .or. any(loc/=[1,1,1,5,7])) stop 6
  if (hy_hydroComputeDtFirstCall .or. extra/=123.) stop 7
  expected=1.
  if (NDIM>1) expected=3.
  if (NDIM>2) expected=4.
  if (abs(hy_lChyp-expected)>1.e-5) stop 8
  hy_geometry=SPHERICAL; dt=100.; hy_lChyp=0.; loc=-1
  call Hydro_computeDt(tile,x,dx,wx,y,dy,wy,z,dz,wz,bounds,guard,state,dt,loc)
  expected=0.4
  if (NDIM>1) expected=0.8/6.
  if (abs(dt-expected)>1.e-5) stop 9
  print *, 'diagnostic wrapper tests passed'
end program
