! Host compatibility entry point for Hydro_computeDt; no magnetic fields.
subroutine Hydro_computeDt(blockDesc,x,dx,uxgrid,y,dy,uygrid,z,dz,uzgrid, &
                           blkLimits,blkLimitsGC,U,dtCheck,dtMinLoc,extraInfo)
  use Hydro_data, only: hy_geometry,hy_lChyp,hy_cfl,hy_meshMe,hy_useHydro, &
                       hy_hydroComputeDtFirstCall,hy_updateHydroFluxes
  use Grid_tile, only: Grid_tile_t
  use newimpl_hydro, only: HY_OK,HY_BAD_INPUT,HY_BAD_STATE
  use newimpl_dt, only: DT_CARTESIAN,DT_CYLINDRICAL,DT_SPHERICAL
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
#include "Simulation.h"
#include "constants.h"
#ifdef SPARK_GLM
#error NewImpl timestep currently supports hydrodynamics only
#endif
  integer, parameter :: rk=kind(1.0)
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
  integer :: block,cellLo(3),cellHi(3),lineLo(3),lineHi(3),ndim,map(6),geometry,level,rank
  integer :: i,j,k,status,failedCell(3),location(5)
  real(rk) :: cfl,sound,inverse(3),rate,speed,maxRate,maxSpeed
real(rk) :: macro_scope_1_local_rho,macro_scope_1_local_p,macro_scope_1_local_g,macro_scope_1_local_a2
real(rk) :: macro_scope_5_local_width(3)
real(rk) :: macro_scope_16_local_velocity(3),macro_scope_16_local_gridVelocity(3),macro_scope_16_local_signal
integer :: macro_scope_16_local_d
real(rk) :: macro_scope_22_local_candidate
  if (.not.hy_useHydro .or. .not.hy_updateHydroFluxes) return
  hy_hydroComputeDtFirstCall=.false.
  block=1; ndim=NDIM; cellLo=blkLimits(LOW,:); cellHi=blkLimits(HIGH,:)
  lineLo=blkLimitsGC(LOW,:); lineHi=blkLimitsGC(HIGH,:)
  map=[DENS_VAR,VELX_VAR,VELY_VAR,VELZ_VAR,PRES_VAR,GAMC_VAR]
  cfl=hy_cfl; level=blockDesc%level; rank=hy_meshMe
  geometry=DT_CARTESIAN
  if (hy_geometry==CYLINDRICAL) geometry=DT_CYLINDRICAL
  if (hy_geometry==SPHERICAL) geometry=DT_SPHERICAL
    status=HY_BAD_INPUT; failedCell=0
    if (ndim<1 .or. ndim>3 .or. block<1 .or. any(cellHi<cellLo)) call Driver_abort( &
  & "NewImpl: invalid timestep calculation")
    if (any(map<1) .or. any(lineLo>cellLo) .or. any(lineHi<cellHi)) call Driver_abort( &
  & "NewImpl: invalid timestep calculation")
    if (geometry<DT_CARTESIAN .or. geometry>DT_SPHERICAL) call Driver_abort("NewImpl: invalid timestep calculation")
    if (.not.ieee_is_finite(cfl) .or. cfl<=0._rk) call Driver_abort("NewImpl: invalid timestep calculation")
    if (.not.all(ieee_is_finite([dtCheck,hy_lChyp]))) call Driver_abort("NewImpl: invalid timestep calculation")
    if (dtCheck<=0._rk .or. hy_lChyp<0._rk) call Driver_abort("NewImpl: invalid timestep calculation")
    maxRate=0._rk; maxSpeed=hy_lChyp; location=0
    do k=cellLo(3),cellHi(3)
      do j=cellLo(2),cellHi(2)
        do i=cellLo(1),cellHi(1)
          failedCell=[i,j,k]
macro_scope_1: do
associate ( &
  macro_scope_1_arg_state => U, &
  macro_scope_1_arg_map => map, &
  macro_scope_1_arg_block => block, &
  macro_scope_1_arg_i => i, &
  macro_scope_1_arg_j => j, &
  macro_scope_1_arg_k => k, &
  macro_scope_1_arg_densityGuard => 0._rk, &
  macro_scope_1_arg_sound => sound, &
  macro_scope_1_arg_status => status)
macro_scope_1_arg_status=HY_BAD_STATE; macro_scope_1_arg_sound=0._rk
macro_scope_1_local_rho=macro_scope_1_arg_state(macro_scope_1_arg_map(1),macro_scope_1_arg_i,macro_scope_1_arg_j, &
  & macro_scope_1_arg_k)
macro_scope_1_local_p=macro_scope_1_arg_state(macro_scope_1_arg_map(5),macro_scope_1_arg_i,macro_scope_1_arg_j, &
  & macro_scope_1_arg_k)
macro_scope_1_local_g=macro_scope_1_arg_state(macro_scope_1_arg_map(6),macro_scope_1_arg_i,macro_scope_1_arg_j, &
  & macro_scope_1_arg_k)
if (.not.all(ieee_is_finite([macro_scope_1_local_rho,macro_scope_1_local_p,macro_scope_1_local_g]))) exit macro_scope_1
macro_scope_1_local_rho=max(macro_scope_1_local_rho,macro_scope_1_arg_densityGuard)
if (macro_scope_1_local_rho<=0._rk) exit macro_scope_1
macro_scope_1_local_a2=macro_scope_1_local_g*macro_scope_1_local_p/macro_scope_1_local_rho
if (.not.ieee_is_finite(macro_scope_1_local_a2)) exit macro_scope_1
if (macro_scope_1_local_a2<0._rk) exit macro_scope_1
macro_scope_1_arg_sound=sqrt(macro_scope_1_local_a2); macro_scope_1_arg_status=HY_OK
end associate
exit macro_scope_1
end do macro_scope_1
          if (status/=HY_OK) call Driver_abort("NewImpl: invalid timestep calculation")
macro_scope_5: do
associate ( &
  macro_scope_5_arg_x => x, &
  macro_scope_5_arg_y => y, &
  macro_scope_5_arg_dx => dx, &
  macro_scope_5_arg_dy => dy, &
  macro_scope_5_arg_dz => dz, &
  macro_scope_5_arg_cellLo => cellLo, &
  macro_scope_5_arg_i => i, &
  macro_scope_5_arg_j => j, &
  macro_scope_5_arg_k => k, &
  macro_scope_5_arg_ndim => ndim, &
  macro_scope_5_arg_geometry => geometry, &
  macro_scope_5_arg_inverse => inverse, &
  macro_scope_5_arg_status => status)
macro_scope_5_arg_status=HY_BAD_INPUT; macro_scope_5_arg_inverse=0._rk; macro_scope_5_local_width=1._rk
macro_scope_5_local_width(1)=macro_scope_5_arg_dx(macro_scope_5_arg_cellLo(1))
if (macro_scope_5_arg_ndim>1) macro_scope_5_local_width(2)=macro_scope_5_arg_dy(macro_scope_5_arg_cellLo(2))
if (macro_scope_5_arg_ndim>2) macro_scope_5_local_width(3)=macro_scope_5_arg_dz(macro_scope_5_arg_cellLo(3))
if (macro_scope_5_arg_geometry==DT_CYLINDRICAL .and. macro_scope_5_arg_ndim>2) macro_scope_5_local_width(3) &
  & =macro_scope_5_arg_x(macro_scope_5_arg_i)*macro_scope_5_arg_dz(macro_scope_5_arg_k)
if (macro_scope_5_arg_geometry==DT_SPHERICAL) then
  if (macro_scope_5_arg_ndim>1) macro_scope_5_local_width(2)=macro_scope_5_arg_x(macro_scope_5_arg_i) &
  & *macro_scope_5_arg_dy(macro_scope_5_arg_j)
  if (macro_scope_5_arg_ndim>2) macro_scope_5_local_width(3)=macro_scope_5_arg_x(macro_scope_5_arg_i)*sin( &
  & macro_scope_5_arg_y(macro_scope_5_arg_j))*macro_scope_5_arg_dz(macro_scope_5_arg_k)
end if
if (.not.all(ieee_is_finite(macro_scope_5_local_width(1:macro_scope_5_arg_ndim)))) exit macro_scope_5
if (any(macro_scope_5_local_width(1:macro_scope_5_arg_ndim)<=0._rk)) exit macro_scope_5
macro_scope_5_arg_inverse(1:macro_scope_5_arg_ndim)=1._rk/macro_scope_5_local_width(1:macro_scope_5_arg_ndim)
macro_scope_5_arg_status=HY_OK
end associate
exit macro_scope_5
end do macro_scope_5
          if (status/=HY_OK) call Driver_abort("NewImpl: invalid timestep calculation")
macro_scope_16: do
associate ( &
  macro_scope_16_arg_state => U, &
  macro_scope_16_arg_map => map, &
  macro_scope_16_arg_block => block, &
  macro_scope_16_arg_i => i, &
  macro_scope_16_arg_j => j, &
  macro_scope_16_arg_k => k, &
  macro_scope_16_arg_ndim => ndim, &
  macro_scope_16_arg_uxgrid => uxgrid, &
  macro_scope_16_arg_uygrid => uygrid, &
  macro_scope_16_arg_uzgrid => uzgrid, &
  macro_scope_16_arg_sound => sound, &
  macro_scope_16_arg_inverse => inverse, &
  macro_scope_16_arg_rate => rate, &
  macro_scope_16_arg_speed => speed, &
  macro_scope_16_arg_status => status)
macro_scope_16_arg_status=HY_BAD_STATE; macro_scope_16_arg_rate=0._rk; macro_scope_16_arg_speed=0._rk
macro_scope_16_local_gridVelocity=0._rk
macro_scope_16_local_gridVelocity(1)=macro_scope_16_arg_uxgrid(macro_scope_16_arg_i)
if (macro_scope_16_arg_ndim>1) macro_scope_16_local_gridVelocity(2)=macro_scope_16_arg_uygrid(macro_scope_16_arg_j)
if (macro_scope_16_arg_ndim>2) macro_scope_16_local_gridVelocity(3)=macro_scope_16_arg_uzgrid(macro_scope_16_arg_k)
do macro_scope_16_local_d=1,macro_scope_16_arg_ndim
  macro_scope_16_local_velocity(macro_scope_16_local_d)=macro_scope_16_arg_state(macro_scope_16_arg_map( &
  & 1+macro_scope_16_local_d),macro_scope_16_arg_i,macro_scope_16_arg_j,macro_scope_16_arg_k)
end do
if (.not.all(ieee_is_finite(macro_scope_16_local_velocity(1:macro_scope_16_arg_ndim)))) exit macro_scope_16
if (.not.all(ieee_is_finite(macro_scope_16_local_gridVelocity(1:macro_scope_16_arg_ndim)))) exit macro_scope_16
do macro_scope_16_local_d=1,macro_scope_16_arg_ndim
  macro_scope_16_local_signal=abs(macro_scope_16_local_velocity(macro_scope_16_local_d) &
  & -macro_scope_16_local_gridVelocity(macro_scope_16_local_d))+macro_scope_16_arg_sound
  macro_scope_16_arg_speed=max(macro_scope_16_arg_speed,macro_scope_16_local_signal)
  macro_scope_16_arg_rate=max(macro_scope_16_arg_rate,macro_scope_16_local_signal*macro_scope_16_arg_inverse( &
  & macro_scope_16_local_d))
end do
if (.not.ieee_is_finite(macro_scope_16_arg_rate) .or. .not.ieee_is_finite(macro_scope_16_arg_speed)) exit macro_scope_16
macro_scope_16_arg_status=HY_OK
end associate
exit macro_scope_16
end do macro_scope_16
          if (status/=HY_OK) call Driver_abort("NewImpl: invalid timestep calculation")
macro_scope_21: do
associate ( &
  macro_scope_21_arg_rate => rate, &
  macro_scope_21_arg_speed => speed, &
  macro_scope_21_arg_i => i, &
  macro_scope_21_arg_j => j, &
  macro_scope_21_arg_k => k, &
  macro_scope_21_arg_level => level, &
  macro_scope_21_arg_rank => rank, &
  macro_scope_21_arg_maxRate => maxRate, &
  macro_scope_21_arg_maxSpeed => maxSpeed, &
  macro_scope_21_arg_location => location)
macro_scope_21_arg_maxSpeed=max(macro_scope_21_arg_maxSpeed,macro_scope_21_arg_speed)
if (macro_scope_21_arg_rate>macro_scope_21_arg_maxRate) then
  macro_scope_21_arg_maxRate=macro_scope_21_arg_rate; macro_scope_21_arg_location=[macro_scope_21_arg_i, &
  & macro_scope_21_arg_j,macro_scope_21_arg_k,macro_scope_21_arg_level,macro_scope_21_arg_rank]
end if
end associate
exit macro_scope_21
end do macro_scope_21
        end do
      end do
    end do
macro_scope_22: do
associate ( &
  macro_scope_22_arg_maxRate => maxRate, &
  macro_scope_22_arg_cfl => cfl, &
  macro_scope_22_arg_location => location, &
  macro_scope_22_arg_dtCheck => dtCheck, &
  macro_scope_22_arg_dtMinLoc => dtMinLoc)
macro_scope_22_local_candidate=huge(1._rk)
if (macro_scope_22_arg_maxRate>0._rk) macro_scope_22_local_candidate=macro_scope_22_arg_cfl/macro_scope_22_arg_maxRate
if (macro_scope_22_local_candidate<macro_scope_22_arg_dtCheck) then
  macro_scope_22_arg_dtCheck=macro_scope_22_local_candidate; macro_scope_22_arg_dtMinLoc=macro_scope_22_arg_location
end if
end associate
exit macro_scope_22
end do macro_scope_22
    hy_lChyp=maxSpeed; status=HY_OK; failedCell=0
end subroutine
