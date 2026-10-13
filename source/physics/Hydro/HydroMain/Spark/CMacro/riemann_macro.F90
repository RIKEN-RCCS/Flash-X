!!****if* source/physics/Hydro/HydroMain/Spark/riemann
!!
!! NAME
!!
!!  hy_rk_riemann
!!
!! SYNOPSIS
!!
!!  hy_rk_riemann( integer(IN) :: dir,
!!               real(IN)    :: VL(HY_VARINUMMAX),
!!               real(IN)    :: VR(HY_VARINUMMAX),
!!               real(OUT)   :: Fstar(HY_VARINUM1),
!!               real(OUT)   :: speed,
!!               integer(OUT):: ierr)
!!
!! ARGUMENTS
!!
!!  dir    - a spatial direction for which the flux is being considered and computed
!!  VL     - primitive variable for left state
!!            (DENS,VELX,VELY,VELZ,PRES,MAGX,MAGY,MAGZ + GAMC,GAME,EINT,TEMP)
!!  VR     - primitive variable for right state
!!            (DENS,VELX,VELY,VELZ,PRES,MAGX,MAGY,MAGZ + GAMC,GAME,EINT,TEMP)
!!  Fstar  - computed flux data
!!            (includes face pressure at the end)
!!  speed  - fastest signal velocity to compute dt
!!  ierr   - a flag to check unphysical negative state (0 is ok; 1 is bad)
!!
!! DESCRIPTION
!!
!!   This routine computes high-order Godunov fluxes based on the left and right Riemann states.
!!
!!   The HLLC Riemann fan:
!!
!!            SL      SM=qStar    SR
!!             \        |        /
!!              \       |       /
!!               \  UL* | UR*  /
!!                \     |     /
!!                 \    |    /
!!                  \   |   /
!!           UL      \  |  /       UR
!!                    \ | /
!!                     \|/
!!   --------------------------------------
!!
!! REFERENCES
!!
!!  * S. Li, JCP, 203:344-357, 2005
!!  * Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics, Springer, 1997
!!
!!***
!!! A note: the speed dummy var seems to be unused. Might be deleted?


subroutine riemann(i1,i2,i3,dir,uPlusArray,uMinusArray,shck,flux, &
     hy_tiny,hy_hybridRiemann,hy_C_hyp)
  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Spark.h"

  !! Arguments type declaration -----------
  real, dimension(HY_NUM_VARS) :: VL, VR
  logical :: inShock
  real, pointer:: Fstar(:)
  real:: speed
  integer :: ierr
  integer, intent(IN) :: i1,i2,i3
  real, dimension(:,:,:,:), intent(INOUT) :: uPlusArray, uMinusArray
  real, dimension(:,:,:), intent(INOUT) :: shck
  real, dimension(:,:,:,:), intent(INOUT), target :: flux
  real, intent(IN) :: hy_tiny, hy_C_hyp
  logical, intent(IN) :: hy_hybridRiemann
  !! --------------------------------------

  real :: SL,SR,cfL,cfR,aL2,aR2,velNL,velNR
  real :: dStarL,dStarR,totalPresL,totalPresR
  real :: BxStar,ByStar,BzStar,Bn_hll,pStar,qStar
  real :: denomL,denomR,numerL,numerR
  real, dimension(HY_NUM_FLUX) :: UL,UR,Uhll,UCstarR,UCstarL
  real, dimension(HY_NUM_FLUX) :: FL,FR
  real :: magBL2,magBR2,magNL,magNR
  real :: Bn_glm, Psi_glm, u2, B2, E, Ptot, UB
  integer, intent(IN) :: dir
!$omp declare target

  inShock = any(shck(@spc(i1-1:i1,i2,i3)@) /= 0.0)
  ! Alias things
! ----------------------------- Make this more effective than a copy! ----------------------------------------------!
  VL = uPlusArray(@blkoff(1:HY_NUM_VARS,i1,-1,i2,+0,i3,+0)@)
  VR = uMinusArray(@blk(1:HY_NUM_VARS,i1,i2,i3)@)
  Fstar => flux(@blk(1:HY_NUM_FLUX,i1,i2,i3)@)
  ! Set no error to begin with
  ierr = 0

  ! Normal velocity
  velNL = VL(HY_VELX+dir-1)
  velNR = VR(HY_VELX+dir-1)

  ! Set sound speed
  aL2   = VL(HY_GAMC)*VL(HY_PRES)/VL(HY_DENS)
  aR2   = VR(HY_GAMC)*VR(HY_PRES)/VR(HY_DENS)

  ! Set zero magnetic quantities by default for hydro
  magNL = 0.
  magNR = 0.
#ifdef SPARK_GLM /* compute additional MHD waves */
  ! Reset the left and right states for the GLM evolution equations
  ! This yields the exact solution for the fluxes of the 2x2 lienar GLM system
  ! See Mignone & Tzeferacos 2010, sec. 3.2
  Bn_glm  = 0.5*(VL(HY_MAGX+dir-1)+VR(HY_MAGX+dir-1)) - 0.5/hy_C_hyp*(VR(HY_PSIB) - VL(HY_PSIB))
  Psi_glm = 0.5*(VL(HY_PSIB)+VR(HY_PSIB)) - 0.5*hy_C_hyp*(VR(HY_MAGX+dir-1)-VL(HY_MAGX+dir-1))
  magNL = VL(HY_MAGX+dir-1)
  magNR = VR(HY_MAGX+dir-1)
  magBL2= dot_product(VL(HY_MAGX:HY_MAGZ),VL(HY_MAGX:HY_MAGZ))/VL(HY_DENS)
  magBR2= dot_product(VR(HY_MAGX:HY_MAGZ),VR(HY_MAGX:HY_MAGZ))/VR(HY_DENS)
#endif

  ! Check unphysical negativity
  if ((VL(HY_DENS) < hy_tiny .and. VL(HY_DENS) > 0.) .or. &
        (VR(HY_DENS) < hy_tiny .and. VR(HY_DENS) > 0.) .or. &
        (VL(HY_PRES) < hy_tiny .and. VL(HY_PRES) > 0.) .or. &
        (VR(HY_PRES) < hy_tiny .and. VR(HY_PRES) > 0.)) then
      ! This could be vacuum limit. We return with zero flux.
      Fstar = 0.
      return
  elseif (aL2 < 0. .or. aR2 < 0.) then
      ierr = 1
      return
  endif

  cfL = sqrt(aL2)
  cfR = sqrt(aR2)
#ifdef SPARK_GLM
  cfL = sqrt(0.5*(aL2 + magBL2 + sqrt((aL2 + magBL2 )**2 - 4.*aL2*magNL*magNL/VL(HY_DENS))))
  cfR = sqrt(0.5*(aR2 + magBR2 + sqrt((aR2 + magBR2 )**2 - 4.*aR2*magNR*magNR/VR(HY_DENS))))
#endif

  ! Get left/right going fastest wave speeds SL & SR for the left and right states
  ! by S. F. Davis, SIAM J. Sci. Stat, Comput., 9(1988) 445.
  ! Also see Miyoshi, Kusano, JCP, 208 (2005)
  SL = min(velNL - cfL, velNR - cfR)
  SR = max(velNL + cfL, velNR + cfR)

  ! Output maximum local wave speed for dt calculation
  speed = max(abs(SL),abs(SR))

  ! Total pressure
  totalPresL = VL(HY_PRES)
  totalPresR = VR(HY_PRES)
#ifdef SPARK_GLM /* for MHD */
  totalPresL = totalPresL + 0.5*dot_product(VL(HY_MAGX:HY_MAGZ),VL(HY_MAGX:HY_MAGZ))
  totalPresR = totalPresR + 0.5*dot_product(VR(HY_MAGX:HY_MAGZ),VR(HY_MAGX:HY_MAGZ))
#endif


! call prim2con(VL,UL)
! call prim2con(VR,UR)
! call prim2flx(dir,VL,FL)
! call prim2flx(dir,VR,FR)
  ! Convert primitive variables to conservative variables
  ! Prim2con L
UL = 0.0
u2 = dot_product(VL(HY_VELX:HY_VELZ),VL(HY_VELX:HY_VELZ))
B2 = 0.
#ifdef SPARK_GLM
B2 = dot_product(VL(HY_MAGX:HY_MAGZ),VL(HY_MAGX:HY_MAGZ))
UL(HY_FMGX:HY_FMGZ) = VL(HY_MAGX:HY_MAGZ)
UL(HY_FPSI) = VL(HY_PSIB)
#endif
UL(HY_MASS) = VL(HY_DENS)
UL(HY_XMOM:HY_ZMOM) = VL(HY_DENS)*VL(HY_VELX:HY_VELZ)
UL(HY_ENER) = 0.5*VL(HY_DENS)*u2 + VL(HY_RHOE) + 0.5*B2

  ! Prim2con R
u2 = dot_product(VR(HY_VELX:HY_VELZ),VR(HY_VELX:HY_VELZ))
B2 = 0.
#ifdef SPARK_GLM
B2 = dot_product(VR(HY_MAGX:HY_MAGZ),VR(HY_MAGX:HY_MAGZ))
UR(HY_FMGX:HY_FMGZ) = VR(HY_MAGX:HY_MAGZ)
UR(HY_FPSI) = VR(HY_PSIB)
#endif
UR(HY_MASS) = VR(HY_DENS)
UR(HY_XMOM:HY_ZMOM) = VR(HY_DENS)*VR(HY_VELX:HY_VELZ)
UR(HY_ENER) = 0.5*VR(HY_DENS)*u2 + VR(HY_RHOE) + 0.5*B2

  ! Prim2flx L
FL = 0.0
u2 = dot_product(VL(HY_VELX:HY_VELZ),VL(HY_VELX:HY_VELZ))
E   = 0.5*VL(HY_DENS)*u2 + VL(HY_RHOE)
Ptot = VL(HY_PRES)

#ifdef SPARK_GLM
B2 = dot_product(VL(HY_MAGX:HY_MAGZ),VL(HY_MAGX:HY_MAGZ))
UB = dot_product(VL(HY_VELX:HY_VELZ),VL(HY_MAGX:HY_MAGZ))
! We will NEED to check units. That could be a pain. #MHDbeNatural
Ptot= Ptot + 0.5*B2
E   = E + 0.5*B2

select case(dir)
case (IAXIS)
   FL(HY_MASS) = VL(HY_DENS)*VL(HY_VELX)
   FL(HY_XMOM) = FL(HY_MASS)*VL(HY_VELX) - VL(HY_MAGX)*VL(HY_MAGX) + Ptot
   FL(HY_YMOM) = FL(HY_MASS)*VL(HY_VELY) - VL(HY_MAGX)*VL(HY_MAGY)
   FL(HY_ZMOM) = FL(HY_MASS)*VL(HY_VELZ) - VL(HY_MAGX)*VL(HY_MAGZ)
   FL(HY_ENER) = (E + Ptot)*VL(HY_VELX) - VL(HY_MAGX)*UB
   FL(HY_FMGX) = 0.
   FL(HY_FMGY) = VL(HY_VELX)*VL(HY_MAGY)-VL(HY_VELY)*VL(HY_MAGX)
   FL(HY_FMGZ) = VL(HY_VELX)*VL(HY_MAGZ)-VL(HY_VELZ)*VL(HY_MAGX)
case (JAXIS)
   FL(HY_MASS) = VL(HY_DENS)*VL(HY_VELY)
   FL(HY_XMOM) = FL(HY_MASS)*VL(HY_VELX) - VL(HY_MAGY)*VL(HY_MAGX)
   FL(HY_YMOM) = FL(HY_MASS)*VL(HY_VELY) - VL(HY_MAGY)*VL(HY_MAGY) + Ptot
   FL(HY_ZMOM) = FL(HY_MASS)*VL(HY_VELZ) - VL(HY_MAGY)*VL(HY_MAGZ)
   FL(HY_ENER) = (E + Ptot)*VL(HY_VELY) - VL(HY_MAGY)*UB
   FL(HY_FMGX) = VL(HY_VELY)*VL(HY_MAGX) - VL(HY_VELX)*VL(HY_MAGY)
   FL(HY_FMGY) = 0.
   FL(HY_FMGZ) = VL(HY_VELY)*VL(HY_MAGZ)-VL(HY_VELZ)*VL(HY_MAGY)
case (KAXIS)
   FL(HY_MASS) = VL(HY_DENS)*VL(HY_VELZ)
   FL(HY_XMOM) = FL(HY_MASS)*VL(HY_VELX) - VL(HY_MAGZ)*VL(HY_MAGX)
   FL(HY_YMOM) = FL(HY_MASS)*VL(HY_VELY) - VL(HY_MAGZ)*VL(HY_MAGY)
   FL(HY_ZMOM) = FL(HY_MASS)*VL(HY_VELZ) - VL(HY_MAGZ)*VL(HY_MAGZ) + Ptot
   FL(HY_ENER) = (E + Ptot)*VL(HY_VELZ) - VL(HY_MAGZ)*UB
   FL(HY_FMGX) = VL(HY_VELZ)*VL(HY_MAGX) - VL(HY_VELX)*VL(HY_MAGZ)
   FL(HY_FMGY) = VL(HY_VELZ)*VL(HY_MAGY) - VL(HY_VELY)*VL(HY_MAGZ)
   FL(HY_FMGZ) = 0.
end select

#else

select case(dir)
case (IAXIS)
   FL(HY_MASS) = VL(HY_DENS)*VL(HY_VELX)
   FL(HY_XMOM) = FL(HY_MASS)*VL(HY_VELX) + Ptot
   FL(HY_YMOM) = FL(HY_MASS)*VL(HY_VELY)
   FL(HY_ZMOM) = FL(HY_MASS)*VL(HY_VELZ)
   FL(HY_ENER) = (E + Ptot)*VL(HY_VELX)
case (JAXIS)
   FL(HY_MASS) = VL(HY_DENS)*VL(HY_VELY)
   FL(HY_XMOM) = FL(HY_MASS)*VL(HY_VELX)
   FL(HY_YMOM) = FL(HY_MASS)*VL(HY_VELY) + Ptot
   FL(HY_ZMOM) = FL(HY_MASS)*VL(HY_VELZ)
   FL(HY_ENER) = (E + Ptot)*VL(HY_VELY)
case (KAXIS)
   FL(HY_MASS) = VL(HY_DENS)*VL(HY_VELZ)
   FL(HY_XMOM) = FL(HY_MASS)*VL(HY_VELX)
   FL(HY_YMOM) = FL(HY_MASS)*VL(HY_VELY)
   FL(HY_ZMOM) = FL(HY_MASS)*VL(HY_VELZ) + Ptot
   FL(HY_ENER) = (E + Ptot)*VL(HY_VELZ)
end select

#endif

    ! Prim2flx R
FR = 0.0
u2 = dot_product(VR(HY_VELX:HY_VELZ),VR(HY_VELX:HY_VELZ))
E   = 0.5*VR(HY_DENS)*u2 + VR(HY_RHOE)
Ptot = VR(HY_PRES)

#ifdef SPARK_GLM
B2 = dot_product(VR(HY_MAGX:HY_MAGZ),VR(HY_MAGX:HY_MAGZ))
UB = dot_product(VR(HY_VELX:HY_VELZ),VR(HY_MAGX:HY_MAGZ))
! We will NEED to check units. That could be a pain. #MHDbeNatural
Ptot= Ptot + 0.5*B2
E   = E + 0.5*B2

select case(dir)
case (IAXIS)
   FR(HY_MASS) = VR(HY_DENS)*VR(HY_VELX)
   FR(HY_XMOM) = FR(HY_MASS)*VR(HY_VELX) - VR(HY_MAGX)*VR(HY_MAGX) + Ptot
   FR(HY_YMOM) = FR(HY_MASS)*VR(HY_VELY) - VR(HY_MAGX)*VR(HY_MAGY)
   FR(HY_ZMOM) = FR(HY_MASS)*VR(HY_VELZ) - VR(HY_MAGX)*VR(HY_MAGZ)
   FR(HY_ENER) = (E + Ptot)*VR(HY_VELX) - VR(HY_MAGX)*UB
   FR(HY_FMGX) = 0.
   FR(HY_FMGY) = VR(HY_VELX)*VR(HY_MAGY)-VR(HY_VELY)*VR(HY_MAGX)
   FR(HY_FMGZ) = VR(HY_VELX)*VR(HY_MAGZ)-VR(HY_VELZ)*VR(HY_MAGX)
case (JAXIS)
   FR(HY_MASS) = VR(HY_DENS)*VR(HY_VELY)
   FR(HY_XMOM) = FR(HY_MASS)*VR(HY_VELX) - VR(HY_MAGY)*VR(HY_MAGX)
   FR(HY_YMOM) = FR(HY_MASS)*VR(HY_VELY) - VR(HY_MAGY)*VR(HY_MAGY) + Ptot
   FR(HY_ZMOM) = FR(HY_MASS)*VR(HY_VELZ) - VR(HY_MAGY)*VR(HY_MAGZ)
   FR(HY_ENER) = (E + Ptot)*VR(HY_VELY) - VR(HY_MAGY)*UB
   FR(HY_FMGX) = VR(HY_VELY)*VR(HY_MAGX) - VR(HY_VELX)*VR(HY_MAGY)
   FR(HY_FMGY) = 0.
   FR(HY_FMGZ) = VR(HY_VELY)*VR(HY_MAGZ)-VR(HY_VELZ)*VR(HY_MAGY)
case (KAXIS)
   FR(HY_MASS) = VR(HY_DENS)*VR(HY_VELZ)
   FR(HY_XMOM) = FR(HY_MASS)*VR(HY_VELX) - VR(HY_MAGZ)*VR(HY_MAGX)
   FR(HY_YMOM) = FR(HY_MASS)*VR(HY_VELY) - VR(HY_MAGZ)*VR(HY_MAGY)
   FR(HY_ZMOM) = FR(HY_MASS)*VR(HY_VELZ) - VR(HY_MAGZ)*VR(HY_MAGZ) + Ptot
   FR(HY_ENER) = (E + Ptot)*VR(HY_VELZ) - VR(HY_MAGZ)*UB
   FR(HY_FMGX) = VR(HY_VELZ)*VR(HY_MAGX) - VR(HY_VELX)*VR(HY_MAGZ)
   FR(HY_FMGY) = VR(HY_VELZ)*VR(HY_MAGY) - VR(HY_VELY)*VR(HY_MAGZ)
   FR(HY_FMGZ) = 0.
end select

#else

select case(dir)
case (IAXIS)
   FR(HY_MASS) = VR(HY_DENS)*VR(HY_VELX)
   FR(HY_XMOM) = FR(HY_MASS)*VR(HY_VELX) + Ptot
   FR(HY_YMOM) = FR(HY_MASS)*VR(HY_VELY)
   FR(HY_ZMOM) = FR(HY_MASS)*VR(HY_VELZ)
   FR(HY_ENER) = (E + Ptot)*VR(HY_VELX)
case (JAXIS)
   FR(HY_MASS) = VR(HY_DENS)*VR(HY_VELY)
   FR(HY_XMOM) = FR(HY_MASS)*VR(HY_VELX)
   FR(HY_YMOM) = FR(HY_MASS)*VR(HY_VELY) + Ptot
   FR(HY_ZMOM) = FR(HY_MASS)*VR(HY_VELZ)
   FR(HY_ENER) = (E + Ptot)*VR(HY_VELY)
case (KAXIS)
   FR(HY_MASS) = VR(HY_DENS)*VR(HY_VELZ)
   FR(HY_XMOM) = FR(HY_MASS)*VR(HY_VELX)
   FR(HY_YMOM) = FR(HY_MASS)*VR(HY_VELY)
   FR(HY_ZMOM) = FR(HY_MASS)*VR(HY_VELZ) + Ptot
   FR(HY_ENER) = (E + Ptot)*VR(HY_VELZ)
end select

#endif




  ! Get HLL states for later use
  if (SL > 0.) then
      Uhll = UL
  elseif ((SL <= 0.) .and. (SR >= 0.)) then
      Uhll = (SR*UR - SL*UL - FR + FL)/(SR - SL)
      !! Pretty sure the following is not need for pressure already in momentum
      !! fluxes (done in prim2flx). I.e., fP = 1.0
      !  Uhll(HY_DENS+dir) = Uhll(HY_DENS+dir) + &
      !       (  totalPresL - totalPresR  )/(SR - SL) * (1.0-fP)
  else
      Uhll = UR
  endif

#ifdef SPARK_GLM /* for MHD */
  ! Calculate intermediate states --------------
  Bn_hll = Uhll(HY_FMGX+dir-1) !=(SR*magNR-SL*magNL)/(SR-SL)
  BxStar = Uhll(HY_FMGX)       !BxStarL = BxStarR = BxHLL
  ByStar = Uhll(HY_FMGY)       !ByStarL = ByStarR = ByHLL
  BzStar = Uhll(HY_FMGZ)       !BzStarL = BzStarR = BzHLL
#endif

  ! (1) Normal velocity component
  ! qStarL = qStarR = qStar
  qStar=(VR(HY_DENS)*velNR*(SR-velNR) &
        - VL(HY_DENS)*velNL*(SL-velNL)  &
        + totalPresL - totalPresR       &
        - magNL**2   + magNR**2)
  qStar = qStar/(VR(HY_DENS)*(SR-velNR) - VL(HY_DENS)*(SL-velNL))

  ! Convenient parameters
  numerL = SL-velNL
  denomL = SL-qStar ! I believe this could be made an inverse to save FLOPS
  numerR = SR-velNR
  denomR = SR-qStar

  ! (2) Total pressure in the star region
  ! pStarL = pStarR = pStar
  pStar = VL(HY_DENS)*numerL*(qStar-velNL)+totalPresL
#ifdef SPARK_GLM /* for MHD */
  pStar = pStar - magNL**2 + Bn_hll**2
#endif

  ! (3) Density
  dStarL = UL(HY_MASS)*numerL/denomL
  dStarR = UR(HY_MASS)*numerR/denomR

  ! (4) Conserved variables in the two-state (left & right) star regions
  UCstarL(HY_MASS)   = dStarL
  UCstarL(HY_ENER)   = UL(HY_ENER)*numerL/denomL + &
        ((pStar*qStar - totalPresL*velNL))/denomL

  UCstarR(HY_MASS)   = dStarR
  UCstarR(HY_ENER)   = UR(HY_ENER)*numerR/denomR + &
        ((pStar*qStar - totalPresR*velNR))/denomR

#ifdef SPARK_GLM /* for MHD */
  UCstarL(HY_FMGX:HY_FMGZ)= Uhll(HY_FMGX:HY_FMGZ)
  UCstarL(HY_ENER) = UCstarL(HY_ENER) &
        -(Bn_hll*dot_product(Uhll(HY_FMGX:HY_FMGZ),Uhll(HY_XMOM:HY_ZMOM))/Uhll(HY_MASS) &
        -  magNL*dot_product(VL(HY_MAGX:HY_MAGZ),  VL(HY_VELX:HY_VELZ)))/denomL

  UCstarR(HY_FMGX:HY_FMGZ)= Uhll(HY_FMGX:HY_FMGZ)
  UCstarR(HY_ENER) = UCstarR(HY_ENER) &
        -(Bn_hll*dot_product(Uhll(HY_FMGX:HY_FMGZ),Uhll(HY_XMOM:HY_ZMOM))/Uhll(HY_MASS) &
        -  magNR*dot_product(VR(HY_MAGX:HY_MAGZ),  VR(HY_VELX:HY_VELZ)))/denomR
#endif

  select case (dir)
  case (IAXIS)
      UCstarL(HY_XMOM) = dStarL*qStar
      UCstarL(HY_YMOM) = UL(HY_YMOM)*numerL/denomL
      UCstarL(HY_ZMOM) = UL(HY_ZMOM)*numerL/denomL

      UCstarR(HY_XMOM) = dStarR*qStar
      UCstarR(HY_YMOM) = UR(HY_YMOM)*numerR/denomR
      UCstarR(HY_ZMOM) = UR(HY_ZMOM)*numerR/denomR

#ifdef SPARK_GLM /* for MHD */
      UCstarL(HY_YMOM) = UCstarL(HY_YMOM) - (BxStar*ByStar-VL(HY_MAGX)*VL(HY_MAGY))/denomL
      UCstarL(HY_ZMOM) = UCstarL(HY_ZMOM) - (BxStar*BzStar-VL(HY_MAGX)*VL(HY_MAGZ))/denomL

      UCstarR(HY_YMOM) = UCstarR(HY_YMOM) - (BxStar*ByStar-VR(HY_MAGX)*VR(HY_MAGY))/denomR
      UCstarR(HY_ZMOM) = UCstarR(HY_ZMOM) - (BxStar*BzStar-VR(HY_MAGX)*VR(HY_MAGZ))/denomR
#endif

  case (JAXIS)
      UCstarL(HY_XMOM) = UL(HY_XMOM)*numerL/denomL
      UCstarL(HY_YMOM) = dStarL*qStar
      UCstarL(HY_ZMOM) = UL(HY_ZMOM)*numerL/denomL

      UCstarR(HY_XMOM) = UR(HY_XMOM)*numerR/denomR
      UCstarR(HY_YMOM) = dStarR*qStar
      UCstarR(HY_ZMOM) = UR(HY_ZMOM)*numerR/denomR

#ifdef SPARK_GLM /* for MHD */
      UCstarL(HY_XMOM) = UCstarL(HY_XMOM) - (ByStar*BxStar-VL(HY_MAGY)*VL(HY_MAGX))/denomL
      UCstarL(HY_ZMOM) = UCstarL(HY_ZMOM) - (ByStar*BzStar-VL(HY_MAGY)*VL(HY_MAGZ))/denomL

      UCstarR(HY_XMOM) = UCstarR(HY_XMOM) - (ByStar*BxStar-VR(HY_MAGY)*VR(HY_MAGX))/denomR
      UCstarR(HY_ZMOM) = UCstarR(HY_ZMOM) - (ByStar*BzStar-VR(HY_MAGY)*VR(HY_MAGZ))/denomR
#endif

  case (KAXIS)
      UCstarL(HY_XMOM) = UL(HY_XMOM)*numerL/denomL
      UCstarL(HY_YMOM) = UL(HY_YMOM)*numerL/denomL
      UCstarL(HY_ZMOM) = dStarL*qStar

      UCstarR(HY_XMOM) = UR(HY_XMOM)*numerR/denomR
      UCstarR(HY_YMOM) = UR(HY_YMOM)*numerR/denomR
      UCstarR(HY_ZMOM) = dStarR*qStar

#ifdef SPARK_GLM /* for MHD */
      UCstarL(HY_XMOM) = UCstarL(HY_XMOM) - (BzStar*BxStar-VL(HY_MAGZ)*VL(HY_MAGX))/denomL
      UCstarL(HY_YMOM) = UCstarL(HY_YMOM) - (BzStar*ByStar-VL(HY_MAGZ)*VL(HY_MAGY))/denomL

      UCstarR(HY_XMOM) = UCstarR(HY_XMOM) - (BzStar*BxStar-VR(HY_MAGZ)*VR(HY_MAGX))/denomR
      UCstarR(HY_YMOM) = UCstarR(HY_YMOM) - (BzStar*ByStar-VR(HY_MAGZ)*VR(HY_MAGY))/denomR
#endif
  end select


  ! End of calculating HLLC intermediate states ---------------------------

  ! (5) Finally, calculate HLLC fluxes
  if (SL >= 0.) then
      Fstar = FL
  elseif ((SL < 0.).and. (qStar >= 0.)) then
      Fstar = FL + SL*(UCstarL - UL)
  elseif ((qStar <0.) .and. (SR >= 0.)) then
      Fstar = FR + SR*(UCstarR - UR)
  else
      Fstar = FR
  endif

  !! The following are the HLLE fluxes.
  if (hy_hybridRiemann .AND. inShock) then
      if (SL > 0.) then
        Fstar = FL
      elseif (SR < 0.) then
        Fstar = FR
      else !if ((SL <= 0.) .and. (SR >= 0.)) then
        Fstar = (SR*FL - SL*FR + SR*SL*(UR - UL))/(SR - SL)
      endif
  end if

#ifdef SPARK_GLM
  ! The exact fluxes for the 2x2 GLM sub-system
  Fstar(HY_FMGX+dir-1) = Psi_glm
  Fstar(HY_FPSI) = hy_C_hyp*hy_C_hyp*Bn_glm
#endif
! print *, Fstar
! print *, flux(1:HY_NUM_FLUX,i1,i2,i3) 
! if (i1==14 .AND. i2==14) then
!   print *, "GPU: (",i1,i2,")", Fstar, new_line('a'), VL, new_line('a'), VR
! endif
end Subroutine riemann

#include "prim2con.F90"
#include "prim2flx.F90"
  

! subroutine riemannOffloaded_inline(i1,i2,i3,dir,dirLims,guardCells)
! use Hydro_data, only: uPlusArray, uMinusArray, shck, hy_tiny, hy_hybridRiemann, hy_C_hyp, flux
!   implicit none

! #include "constants.h"
! #include "Simulation.h"
! #include "Spark.h"

!   !! Arguments type declaration -----------
!   real, dimension(HY_NUM_VARS) :: VL, VR
!   logical :: inShock
!   real, pointer:: Fstar(:)
!   real:: speed
!   integer, dimension(LOW:HIGH,MDIM),  intent(in) :: dirLims
!   integer, dimension(3), intent(in) :: guardCells
!   integer :: ierr, i1,i2,i3
!   !! --------------------------------------

!   real :: SL,SR,cfL,cfR,aL2,aR2,velNL,velNR
!   real :: dStarL,dStarR,totalPresL,totalPresR
!   real :: BxStar,ByStar,BzStar,Bn_hll,pStar,qStar
!   real :: denomL,denomR,numerL,numerR
!   real, dimension(HY_NUM_FLUX) :: UL,UR,Uhll,UCstarR,UCstarL
!   real, dimension(HY_NUM_FLUX) :: FL,FR
!   real :: magBL2,magBR2,magNL,magNR
!   real :: Bn_glm, Psi_glm, u2, B2, E, Ptot, UB
!   integer :: dir

!   !$omp target teams distribute parallel do collapse(3) &
!   !$omp private(i1,i2,i3) shared(dirLims,dir,guardCells,uPlusArray, uMinusArray, shck, hy_tiny, hy_hybridRiemann, hy_C_hyp, flux) default(private)
!   do i3 = 1 + guardCells(3) ,1 + guardCells(3) + dirLims(HIGH,3) - dirLims(LOW,3)
!     do i2 = 1 + guardCells(2) , 1 + guardCells(2) + dirLims(HIGH,2) - dirLims(LOW,2)
!        do i1 = 1 + guardCells(1) , 1 + guardCells(1) + dirLims(HIGH,1) - dirLims(LOW,1) + 1

!   inShock = any(shck(i1-1:i1,i2,i3) /= 0.0)
!   ! Alias things
! ! ----------------------------- Make this more effective than a copy! ----------------------------------------------!
!   VL = uPlusArray(1:HY_NUM_VARS,i1-1,i2,i3)
!   VR = uMinusArray(1:HY_NUM_VARS,i1,i2,i3)
!   Fstar => flux(1:HY_NUM_FLUX,i1,i2,i3)
!   ! Set no error to begin with
!   ierr = 0

!   ! Normal velocity
!   velNL = VL(HY_VELX+dir-1)
!   velNR = VR(HY_VELX+dir-1)

!   ! Set sound speed
!   aL2   = VL(HY_GAMC)*VL(HY_PRES)/VL(HY_DENS)
!   aR2   = VR(HY_GAMC)*VR(HY_PRES)/VR(HY_DENS)

!   ! Set zero magnetic quantities by default for hydro
!   magNL = 0.
!   magNR = 0.
! #ifdef SPARK_GLM /* compute additional MHD waves */
!   ! Reset the left and right states for the GLM evolution equations
!   ! This yields the exact solution for the fluxes of the 2x2 lienar GLM system
!   ! See Mignone & Tzeferacos 2010, sec. 3.2
!   Bn_glm  = 0.5*(VL(HY_MAGX+dir-1)+VR(HY_MAGX+dir-1)) - 0.5/hy_C_hyp*(VR(HY_PSIB) - VL(HY_PSIB))
!   Psi_glm = 0.5*(VL(HY_PSIB)+VR(HY_PSIB)) - 0.5*hy_C_hyp*(VR(HY_MAGX+dir-1)-VL(HY_MAGX+dir-1))
!   magNL = VL(HY_MAGX+dir-1)
!   magNR = VR(HY_MAGX+dir-1)
!   magBL2= dot_product(VL(HY_MAGX:HY_MAGZ),VL(HY_MAGX:HY_MAGZ))/VL(HY_DENS)
!   magBR2= dot_product(VR(HY_MAGX:HY_MAGZ),VR(HY_MAGX:HY_MAGZ))/VR(HY_DENS)
! #endif

!   ! Check unphysical negativity
!   if ((VL(HY_DENS) < hy_tiny .and. VL(HY_DENS) > 0.) .or. &
!         (VR(HY_DENS) < hy_tiny .and. VR(HY_DENS) > 0.) .or. &
!         (VL(HY_PRES) < hy_tiny .and. VL(HY_PRES) > 0.) .or. &
!         (VR(HY_PRES) < hy_tiny .and. VR(HY_PRES) > 0.)) then
!       ! This could be vacuum limit. We return with zero flux.
!       Fstar = 0.
!       return
!   elseif (aL2 < 0. .or. aR2 < 0.) then
!       ierr = 1
!       return
!   endif

!   cfL = sqrt(aL2)
!   cfR = sqrt(aR2)
! #ifdef SPARK_GLM
!   cfL = sqrt(0.5*(aL2 + magBL2 + sqrt((aL2 + magBL2 )**2 - 4.*aL2*magNL*magNL/VL(HY_DENS))))
!   cfR = sqrt(0.5*(aR2 + magBR2 + sqrt((aR2 + magBR2 )**2 - 4.*aR2*magNR*magNR/VR(HY_DENS))))
! #endif

!   ! Get left/right going fastest wave speeds SL & SR for the left and right states
!   ! by S. F. Davis, SIAM J. Sci. Stat, Comput., 9(1988) 445.
!   ! Also see Miyoshi, Kusano, JCP, 208 (2005)
!   SL = min(velNL - cfL, velNR - cfR)
!   SR = max(velNL + cfL, velNR + cfR)

!   ! Output maximum local wave speed for dt calculation
!   speed = max(abs(SL),abs(SR))

!   ! Total pressure
!   totalPresL = VL(HY_PRES)
!   totalPresR = VR(HY_PRES)
! #ifdef SPARK_GLM /* for MHD */
!   totalPresL = totalPresL + 0.5*dot_product(VL(HY_MAGX:HY_MAGZ),VL(HY_MAGX:HY_MAGZ))
!   totalPresR = totalPresR + 0.5*dot_product(VR(HY_MAGX:HY_MAGZ),VR(HY_MAGX:HY_MAGZ))
! #endif


! ! call prim2con(VL,UL)
! ! call prim2con(VR,UR)
! ! call prim2flx(dir,VL,FL)
! ! call prim2flx(dir,VR,FR)
!   ! Convert primitive variables to conservative variables
!   ! Prim2con L
! UL = 0.0
! u2 = dot_product(VL(HY_VELX:HY_VELZ),VL(HY_VELX:HY_VELZ))
! B2 = 0.
! #ifdef SPARK_GLM
! B2 = dot_product(VL(HY_MAGX:HY_MAGZ),VL(HY_MAGX:HY_MAGZ))
! UL(HY_FMGX:HY_FMGZ) = VL(HY_MAGX:HY_MAGZ)
! UL(HY_FPSI) = VL(HY_PSIB)
! #endif
! UL(HY_MASS) = VL(HY_DENS)
! UL(HY_XMOM:HY_ZMOM) = VL(HY_DENS)*VL(HY_VELX:HY_VELZ)
! UL(HY_ENER) = 0.5*VL(HY_DENS)*u2 + VL(HY_RHOE) + 0.5*B2

!   ! Prim2con R
! u2 = dot_product(VR(HY_VELX:HY_VELZ),VR(HY_VELX:HY_VELZ))
! B2 = 0.
! #ifdef SPARK_GLM
! B2 = dot_product(VR(HY_MAGX:HY_MAGZ),VR(HY_MAGX:HY_MAGZ))
! UR(HY_FMGX:HY_FMGZ) = VR(HY_MAGX:HY_MAGZ)
! UR(HY_FPSI) = VR(HY_PSIB)
! #endif
! UR(HY_MASS) = VR(HY_DENS)
! UR(HY_XMOM:HY_ZMOM) = VR(HY_DENS)*VR(HY_VELX:HY_VELZ)
! UR(HY_ENER) = 0.5*VR(HY_DENS)*u2 + VR(HY_RHOE) + 0.5*B2

!   ! Prim2flx L
! FL = 0.0
! u2 = dot_product(VL(HY_VELX:HY_VELZ),VL(HY_VELX:HY_VELZ))
! E   = 0.5*VL(HY_DENS)*u2 + VL(HY_RHOE)
! Ptot = VL(HY_PRES)

! #ifdef SPARK_GLM
! B2 = dot_product(VL(HY_MAGX:HY_MAGZ),VL(HY_MAGX:HY_MAGZ))
! UB = dot_product(VL(HY_VELX:HY_VELZ),VL(HY_MAGX:HY_MAGZ))
! ! We will NEED to check units. That could be a pain. #MHDbeNatural
! Ptot= Ptot + 0.5*B2
! E   = E + 0.5*B2

! select case(dir)
! case (IAXIS)
!    FL(HY_MASS) = VL(HY_DENS)*VL(HY_VELX)
!    FL(HY_XMOM) = FL(HY_MASS)*VL(HY_VELX) - VL(HY_MAGX)*VL(HY_MAGX) + Ptot
!    FL(HY_YMOM) = FL(HY_MASS)*VL(HY_VELY) - VL(HY_MAGX)*VL(HY_MAGY)
!    FL(HY_ZMOM) = FL(HY_MASS)*VL(HY_VELZ) - VL(HY_MAGX)*VL(HY_MAGZ)
!    FL(HY_ENER) = (E + Ptot)*VL(HY_VELX) - VL(HY_MAGX)*UB
!    FL(HY_FMGX) = 0.
!    FL(HY_FMGY) = VL(HY_VELX)*VL(HY_MAGY)-VL(HY_VELY)*VL(HY_MAGX)
!    FL(HY_FMGZ) = VL(HY_VELX)*VL(HY_MAGZ)-VL(HY_VELZ)*VL(HY_MAGX)
! case (JAXIS)
!    FL(HY_MASS) = VL(HY_DENS)*VL(HY_VELY)
!    FL(HY_XMOM) = FL(HY_MASS)*VL(HY_VELX) - VL(HY_MAGY)*VL(HY_MAGX)
!    FL(HY_YMOM) = FL(HY_MASS)*VL(HY_VELY) - VL(HY_MAGY)*VL(HY_MAGY) + Ptot
!    FL(HY_ZMOM) = FL(HY_MASS)*VL(HY_VELZ) - VL(HY_MAGY)*VL(HY_MAGZ)
!    FL(HY_ENER) = (E + Ptot)*VL(HY_VELY) - VL(HY_MAGY)*UB
!    FL(HY_FMGX) = VL(HY_VELY)*VL(HY_MAGX) - VL(HY_VELX)*VL(HY_MAGY)
!    FL(HY_FMGY) = 0.
!    FL(HY_FMGZ) = VL(HY_VELY)*VL(HY_MAGZ)-VL(HY_VELZ)*VL(HY_MAGY)
! case (KAXIS)
!    FL(HY_MASS) = VL(HY_DENS)*VL(HY_VELZ)
!    FL(HY_XMOM) = FL(HY_MASS)*VL(HY_VELX) - VL(HY_MAGZ)*VL(HY_MAGX)
!    FL(HY_YMOM) = FL(HY_MASS)*VL(HY_VELY) - VL(HY_MAGZ)*VL(HY_MAGY)
!    FL(HY_ZMOM) = FL(HY_MASS)*VL(HY_VELZ) - VL(HY_MAGZ)*VL(HY_MAGZ) + Ptot
!    FL(HY_ENER) = (E + Ptot)*VL(HY_VELZ) - VL(HY_MAGZ)*UB
!    FL(HY_FMGX) = VL(HY_VELZ)*VL(HY_MAGX) - VL(HY_VELX)*VL(HY_MAGZ)
!    FL(HY_FMGY) = VL(HY_VELZ)*VL(HY_MAGY) - VL(HY_VELY)*VL(HY_MAGZ)
!    FL(HY_FMGZ) = 0.
! end select

! #else

! select case(dir)
! case (IAXIS)
!    FL(HY_MASS) = VL(HY_DENS)*VL(HY_VELX)
!    FL(HY_XMOM) = FL(HY_MASS)*VL(HY_VELX) + Ptot
!    FL(HY_YMOM) = FL(HY_MASS)*VL(HY_VELY)
!    FL(HY_ZMOM) = FL(HY_MASS)*VL(HY_VELZ)
!    FL(HY_ENER) = (E + Ptot)*VL(HY_VELX)
! case (JAXIS)
!    FL(HY_MASS) = VL(HY_DENS)*VL(HY_VELY)
!    FL(HY_XMOM) = FL(HY_MASS)*VL(HY_VELX)
!    FL(HY_YMOM) = FL(HY_MASS)*VL(HY_VELY) + Ptot
!    FL(HY_ZMOM) = FL(HY_MASS)*VL(HY_VELZ)
!    FL(HY_ENER) = (E + Ptot)*VL(HY_VELY)
! case (KAXIS)
!    FL(HY_MASS) = VL(HY_DENS)*VL(HY_VELZ)
!    FL(HY_XMOM) = FL(HY_MASS)*VL(HY_VELX)
!    FL(HY_YMOM) = FL(HY_MASS)*VL(HY_VELY)
!    FL(HY_ZMOM) = FL(HY_MASS)*VL(HY_VELZ) + Ptot
!    FL(HY_ENER) = (E + Ptot)*VL(HY_VELZ)
! end select

! #endif

!     ! Prim2flx R
! FR = 0.0
! u2 = dot_product(VR(HY_VELX:HY_VELZ),VR(HY_VELX:HY_VELZ))
! E   = 0.5*VR(HY_DENS)*u2 + VR(HY_RHOE)
! Ptot = VR(HY_PRES)

! #ifdef SPARK_GLM
! B2 = dot_product(VR(HY_MAGX:HY_MAGZ),VR(HY_MAGX:HY_MAGZ))
! UB = dot_product(VR(HY_VELX:HY_VELZ),VR(HY_MAGX:HY_MAGZ))
! ! We will NEED to check units. That could be a pain. #MHDbeNatural
! Ptot= Ptot + 0.5*B2
! E   = E + 0.5*B2

! select case(dir)
! case (IAXIS)
!    FR(HY_MASS) = VR(HY_DENS)*VR(HY_VELX)
!    FR(HY_XMOM) = FR(HY_MASS)*VR(HY_VELX) - VR(HY_MAGX)*VR(HY_MAGX) + Ptot
!    FR(HY_YMOM) = FR(HY_MASS)*VR(HY_VELY) - VR(HY_MAGX)*VR(HY_MAGY)
!    FR(HY_ZMOM) = FR(HY_MASS)*VR(HY_VELZ) - VR(HY_MAGX)*VR(HY_MAGZ)
!    FR(HY_ENER) = (E + Ptot)*VR(HY_VELX) - VR(HY_MAGX)*UB
!    FR(HY_FMGX) = 0.
!    FR(HY_FMGY) = VR(HY_VELX)*VR(HY_MAGY)-VR(HY_VELY)*VR(HY_MAGX)
!    FR(HY_FMGZ) = VR(HY_VELX)*VR(HY_MAGZ)-VR(HY_VELZ)*VR(HY_MAGX)
! case (JAXIS)
!    FR(HY_MASS) = VR(HY_DENS)*VR(HY_VELY)
!    FR(HY_XMOM) = FR(HY_MASS)*VR(HY_VELX) - VR(HY_MAGY)*VR(HY_MAGX)
!    FR(HY_YMOM) = FR(HY_MASS)*VR(HY_VELY) - VR(HY_MAGY)*VR(HY_MAGY) + Ptot
!    FR(HY_ZMOM) = FR(HY_MASS)*VR(HY_VELZ) - VR(HY_MAGY)*VR(HY_MAGZ)
!    FR(HY_ENER) = (E + Ptot)*VR(HY_VELY) - VR(HY_MAGY)*UB
!    FR(HY_FMGX) = VR(HY_VELY)*VR(HY_MAGX) - VR(HY_VELX)*VR(HY_MAGY)
!    FR(HY_FMGY) = 0.
!    FR(HY_FMGZ) = VR(HY_VELY)*VR(HY_MAGZ)-VR(HY_VELZ)*VR(HY_MAGY)
! case (KAXIS)
!    FR(HY_MASS) = VR(HY_DENS)*VR(HY_VELZ)
!    FR(HY_XMOM) = FR(HY_MASS)*VR(HY_VELX) - VR(HY_MAGZ)*VR(HY_MAGX)
!    FR(HY_YMOM) = FR(HY_MASS)*VR(HY_VELY) - VR(HY_MAGZ)*VR(HY_MAGY)
!    FR(HY_ZMOM) = FR(HY_MASS)*VR(HY_VELZ) - VR(HY_MAGZ)*VR(HY_MAGZ) + Ptot
!    FR(HY_ENER) = (E + Ptot)*VR(HY_VELZ) - VR(HY_MAGZ)*UB
!    FR(HY_FMGX) = VR(HY_VELZ)*VR(HY_MAGX) - VR(HY_VELX)*VR(HY_MAGZ)
!    FR(HY_FMGY) = VR(HY_VELZ)*VR(HY_MAGY) - VR(HY_VELY)*VR(HY_MAGZ)
!    FR(HY_FMGZ) = 0.
! end select

! #else

! select case(dir)
! case (IAXIS)
!    FR(HY_MASS) = VR(HY_DENS)*VR(HY_VELX)
!    FR(HY_XMOM) = FR(HY_MASS)*VR(HY_VELX) + Ptot
!    FR(HY_YMOM) = FR(HY_MASS)*VR(HY_VELY)
!    FR(HY_ZMOM) = FR(HY_MASS)*VR(HY_VELZ)
!    FR(HY_ENER) = (E + Ptot)*VR(HY_VELX)
! case (JAXIS)
!    FR(HY_MASS) = VR(HY_DENS)*VR(HY_VELY)
!    FR(HY_XMOM) = FR(HY_MASS)*VR(HY_VELX)
!    FR(HY_YMOM) = FR(HY_MASS)*VR(HY_VELY) + Ptot
!    FR(HY_ZMOM) = FR(HY_MASS)*VR(HY_VELZ)
!    FR(HY_ENER) = (E + Ptot)*VR(HY_VELY)
! case (KAXIS)
!    FR(HY_MASS) = VR(HY_DENS)*VR(HY_VELZ)
!    FR(HY_XMOM) = FR(HY_MASS)*VR(HY_VELX)
!    FR(HY_YMOM) = FR(HY_MASS)*VR(HY_VELY)
!    FR(HY_ZMOM) = FR(HY_MASS)*VR(HY_VELZ) + Ptot
!    FR(HY_ENER) = (E + Ptot)*VR(HY_VELZ)
! end select

! #endif




!   ! Get HLL states for later use
!   if (SL > 0.) then
!       Uhll = UL
!   elseif ((SL <= 0.) .and. (SR >= 0.)) then
!       Uhll = (SR*UR - SL*UL - FR + FL)/(SR - SL)
!       !! Pretty sure the following is not need for pressure already in momentum
!       !! fluxes (done in prim2flx). I.e., fP = 1.0
!       !  Uhll(HY_DENS+dir) = Uhll(HY_DENS+dir) + &
!       !       (  totalPresL - totalPresR  )/(SR - SL) * (1.0-fP)
!   else
!       Uhll = UR
!   endif

! #ifdef SPARK_GLM /* for MHD */
!   ! Calculate intermediate states --------------
!   Bn_hll = Uhll(HY_FMGX+dir-1) !=(SR*magNR-SL*magNL)/(SR-SL)
!   BxStar = Uhll(HY_FMGX)       !BxStarL = BxStarR = BxHLL
!   ByStar = Uhll(HY_FMGY)       !ByStarL = ByStarR = ByHLL
!   BzStar = Uhll(HY_FMGZ)       !BzStarL = BzStarR = BzHLL
! #endif

!   ! (1) Normal velocity component
!   ! qStarL = qStarR = qStar
!   qStar=(VR(HY_DENS)*velNR*(SR-velNR) &
!         - VL(HY_DENS)*velNL*(SL-velNL)  &
!         + totalPresL - totalPresR       &
!         - magNL**2   + magNR**2)
!   qStar = qStar/(VR(HY_DENS)*(SR-velNR) - VL(HY_DENS)*(SL-velNL))

!   ! Convenient parameters
!   numerL = SL-velNL
!   denomL = SL-qStar ! I believe this could be made an inverse to save FLOPS
!   numerR = SR-velNR
!   denomR = SR-qStar

!   ! (2) Total pressure in the star region
!   ! pStarL = pStarR = pStar
!   pStar = VL(HY_DENS)*numerL*(qStar-velNL)+totalPresL
! #ifdef SPARK_GLM /* for MHD */
!   pStar = pStar - magNL**2 + Bn_hll**2
! #endif

!   ! (3) Density
!   dStarL = UL(HY_MASS)*numerL/denomL
!   dStarR = UR(HY_MASS)*numerR/denomR

!   ! (4) Conserved variables in the two-state (left & right) star regions
!   UCstarL(HY_MASS)   = dStarL
!   UCstarL(HY_ENER)   = UL(HY_ENER)*numerL/denomL + &
!         ((pStar*qStar - totalPresL*velNL))/denomL

!   UCstarR(HY_MASS)   = dStarR
!   UCstarR(HY_ENER)   = UR(HY_ENER)*numerR/denomR + &
!         ((pStar*qStar - totalPresR*velNR))/denomR

! #ifdef SPARK_GLM /* for MHD */
!   UCstarL(HY_FMGX:HY_FMGZ)= Uhll(HY_FMGX:HY_FMGZ)
!   UCstarL(HY_ENER) = UCstarL(HY_ENER) &
!         -(Bn_hll*dot_product(Uhll(HY_FMGX:HY_FMGZ),Uhll(HY_XMOM:HY_ZMOM))/Uhll(HY_MASS) &
!         -  magNL*dot_product(VL(HY_MAGX:HY_MAGZ),  VL(HY_VELX:HY_VELZ)))/denomL

!   UCstarR(HY_FMGX:HY_FMGZ)= Uhll(HY_FMGX:HY_FMGZ)
!   UCstarR(HY_ENER) = UCstarR(HY_ENER) &
!         -(Bn_hll*dot_product(Uhll(HY_FMGX:HY_FMGZ),Uhll(HY_XMOM:HY_ZMOM))/Uhll(HY_MASS) &
!         -  magNR*dot_product(VR(HY_MAGX:HY_MAGZ),  VR(HY_VELX:HY_VELZ)))/denomR
! #endif

!   select case (dir)
!   case (IAXIS)
!       UCstarL(HY_XMOM) = dStarL*qStar
!       UCstarL(HY_YMOM) = UL(HY_YMOM)*numerL/denomL
!       UCstarL(HY_ZMOM) = UL(HY_ZMOM)*numerL/denomL

!       UCstarR(HY_XMOM) = dStarR*qStar
!       UCstarR(HY_YMOM) = UR(HY_YMOM)*numerR/denomR
!       UCstarR(HY_ZMOM) = UR(HY_ZMOM)*numerR/denomR

! #ifdef SPARK_GLM /* for MHD */
!       UCstarL(HY_YMOM) = UCstarL(HY_YMOM) - (BxStar*ByStar-VL(HY_MAGX)*VL(HY_MAGY))/denomL
!       UCstarL(HY_ZMOM) = UCstarL(HY_ZMOM) - (BxStar*BzStar-VL(HY_MAGX)*VL(HY_MAGZ))/denomL

!       UCstarR(HY_YMOM) = UCstarR(HY_YMOM) - (BxStar*ByStar-VR(HY_MAGX)*VR(HY_MAGY))/denomR
!       UCstarR(HY_ZMOM) = UCstarR(HY_ZMOM) - (BxStar*BzStar-VR(HY_MAGX)*VR(HY_MAGZ))/denomR
! #endif

!   case (JAXIS)
!       UCstarL(HY_XMOM) = UL(HY_XMOM)*numerL/denomL
!       UCstarL(HY_YMOM) = dStarL*qStar
!       UCstarL(HY_ZMOM) = UL(HY_ZMOM)*numerL/denomL

!       UCstarR(HY_XMOM) = UR(HY_XMOM)*numerR/denomR
!       UCstarR(HY_YMOM) = dStarR*qStar
!       UCstarR(HY_ZMOM) = UR(HY_ZMOM)*numerR/denomR

! #ifdef SPARK_GLM /* for MHD */
!       UCstarL(HY_XMOM) = UCstarL(HY_XMOM) - (ByStar*BxStar-VL(HY_MAGY)*VL(HY_MAGX))/denomL
!       UCstarL(HY_ZMOM) = UCstarL(HY_ZMOM) - (ByStar*BzStar-VL(HY_MAGY)*VL(HY_MAGZ))/denomL

!       UCstarR(HY_XMOM) = UCstarR(HY_XMOM) - (ByStar*BxStar-VR(HY_MAGY)*VR(HY_MAGX))/denomR
!       UCstarR(HY_ZMOM) = UCstarR(HY_ZMOM) - (ByStar*BzStar-VR(HY_MAGY)*VR(HY_MAGZ))/denomR
! #endif

!   case (KAXIS)
!       UCstarL(HY_XMOM) = UL(HY_XMOM)*numerL/denomL
!       UCstarL(HY_YMOM) = UL(HY_YMOM)*numerL/denomL
!       UCstarL(HY_ZMOM) = dStarL*qStar

!       UCstarR(HY_XMOM) = UR(HY_XMOM)*numerR/denomR
!       UCstarR(HY_YMOM) = UR(HY_YMOM)*numerR/denomR
!       UCstarR(HY_ZMOM) = dStarR*qStar

! #ifdef SPARK_GLM /* for MHD */
!       UCstarL(HY_XMOM) = UCstarL(HY_XMOM) - (BzStar*BxStar-VL(HY_MAGZ)*VL(HY_MAGX))/denomL
!       UCstarL(HY_YMOM) = UCstarL(HY_YMOM) - (BzStar*ByStar-VL(HY_MAGZ)*VL(HY_MAGY))/denomL

!       UCstarR(HY_XMOM) = UCstarR(HY_XMOM) - (BzStar*BxStar-VR(HY_MAGZ)*VR(HY_MAGX))/denomR
!       UCstarR(HY_YMOM) = UCstarR(HY_YMOM) - (BzStar*ByStar-VR(HY_MAGZ)*VR(HY_MAGY))/denomR
! #endif
!   end select


!   ! End of calculating HLLC intermediate states ---------------------------

!   ! (5) Finally, calculate HLLC fluxes
!   if (SL >= 0.) then
!       Fstar = FL
!   elseif ((SL < 0.).and. (qStar >= 0.)) then
!       Fstar = FL + SL*(UCstarL - UL)
!   elseif ((qStar <0.) .and. (SR >= 0.)) then
!       Fstar = FR + SR*(UCstarR - UR)
!   else
!       Fstar = FR
!   endif

!   !! The following are the HLLE fluxes.
!   if (hy_hybridRiemann .AND. inShock) then
!       if (SL > 0.) then
!         Fstar = FL
!       elseif (SR < 0.) then
!         Fstar = FR
!       else !if ((SL <= 0.) .and. (SR >= 0.)) then
!         Fstar = (SR*FL - SL*FR + SR*SL*(UR - UL))/(SR - SL)
!       endif
!   end if

! #ifdef SPARK_GLM
!   ! The exact fluxes for the 2x2 GLM sub-system
!   Fstar(HY_FMGX+dir-1) = Psi_glm
!   Fstar(HY_FPSI) = hy_C_hyp*hy_C_hyp*Bn_glm
! #endif
! ! print *, Fstar
! ! print *, flux(1:HY_NUM_FLUX,i1,i2,i3) 
! ! if (i1==14 .AND. i2==14) then
! !   print *, "GPU: (",i1,i2,")", Fstar, new_line('a'), VL, new_line('a'), VR
! ! endif
! enddo
! enddo
! enddo
! end Subroutine riemannOffloaded_inline
