!!****if* source/physics/Hydro/HydroMain/Spark/prim2flx
!!
!! NAME
!!
!!  hy_rk_prim2flx
!!
!! SYNOPSIS
!!
!!  hy_rk_prim2flx( integer(IN) :: dir,
!!                   real(IN)    :: V(HY_VARINUM3),
!!                   real(OUT)   :: F(HY_VARINUM1))
!!
!! ARGUMENTS
!!
!! dir- directional index
!! V  - primitive variables  + GAMC,GAME
!! F  - flux
!!
!! DESCRIPTION
!!
!!  This routine calculates conversion from primitive variables to fluxes.
!!
!!***

subroutine prim2flx(dir,V,F)
    implicit none
  
#include "constants.h"
#include "Simulation.h"
#include "Spark.h"
!$omp declare target
   !! Arguments type declaration -----------
   integer, intent(IN)  :: dir
   real, dimension(HY_NUM_VARS), intent(IN) :: V
   real, dimension(HY_NUM_FLUX), intent(OUT) :: F
   !! --------------------------------------

   real  :: u2,E,B2,UB,Ptot

   F = 0.0
   u2 = dot_product(V(HY_VELX:HY_VELZ),V(HY_VELX:HY_VELZ))
   E   = 0.5*V(HY_DENS)*u2 + V(HY_RHOE)
   Ptot = V(HY_PRES)

#ifdef SPARK_GLM
   B2 = dot_product(V(HY_MAGX:HY_MAGZ),V(HY_MAGX:HY_MAGZ))
   UB = dot_product(V(HY_VELX:HY_VELZ),V(HY_MAGX:HY_MAGZ))
   ! We will NEED to check units. That could be a pain. #MHDbeNatural
   Ptot= Ptot + 0.5*B2
   E   = E + 0.5*B2

   select case(dir)
   case (IAXIS)
      F(HY_MASS) = V(HY_DENS)*V(HY_VELX)
      F(HY_XMOM) = F(HY_MASS)*V(HY_VELX) - V(HY_MAGX)*V(HY_MAGX) + Ptot
      F(HY_YMOM) = F(HY_MASS)*V(HY_VELY) - V(HY_MAGX)*V(HY_MAGY)
      F(HY_ZMOM) = F(HY_MASS)*V(HY_VELZ) - V(HY_MAGX)*V(HY_MAGZ)
      F(HY_ENER) = (E + Ptot)*V(HY_VELX) - V(HY_MAGX)*UB
      F(HY_FMGX) = 0.
      F(HY_FMGY) = V(HY_VELX)*V(HY_MAGY)-V(HY_VELY)*V(HY_MAGX)
      F(HY_FMGZ) = V(HY_VELX)*V(HY_MAGZ)-V(HY_VELZ)*V(HY_MAGX)
   case (JAXIS)
      F(HY_MASS) = V(HY_DENS)*V(HY_VELY)
      F(HY_XMOM) = F(HY_MASS)*V(HY_VELX) - V(HY_MAGY)*V(HY_MAGX)
      F(HY_YMOM) = F(HY_MASS)*V(HY_VELY) - V(HY_MAGY)*V(HY_MAGY) + Ptot
      F(HY_ZMOM) = F(HY_MASS)*V(HY_VELZ) - V(HY_MAGY)*V(HY_MAGZ)
      F(HY_ENER) = (E + Ptot)*V(HY_VELY) - V(HY_MAGY)*UB
      F(HY_FMGX) = V(HY_VELY)*V(HY_MAGX) - V(HY_VELX)*V(HY_MAGY)
      F(HY_FMGY) = 0.
      F(HY_FMGZ) = V(HY_VELY)*V(HY_MAGZ)-V(HY_VELZ)*V(HY_MAGY)
   case (KAXIS)
      F(HY_MASS) = V(HY_DENS)*V(HY_VELZ)
      F(HY_XMOM) = F(HY_MASS)*V(HY_VELX) - V(HY_MAGZ)*V(HY_MAGX)
      F(HY_YMOM) = F(HY_MASS)*V(HY_VELY) - V(HY_MAGZ)*V(HY_MAGY)
      F(HY_ZMOM) = F(HY_MASS)*V(HY_VELZ) - V(HY_MAGZ)*V(HY_MAGZ) + Ptot
      F(HY_ENER) = (E + Ptot)*V(HY_VELZ) - V(HY_MAGZ)*UB
      F(HY_FMGX) = V(HY_VELZ)*V(HY_MAGX) - V(HY_VELX)*V(HY_MAGZ)
      F(HY_FMGY) = V(HY_VELZ)*V(HY_MAGY) - V(HY_VELY)*V(HY_MAGZ)
      F(HY_FMGZ) = 0.
   end select

#else

   select case(dir)
   case (IAXIS)
      F(HY_MASS) = V(HY_DENS)*V(HY_VELX)
      F(HY_XMOM) = F(HY_MASS)*V(HY_VELX) + Ptot
      F(HY_YMOM) = F(HY_MASS)*V(HY_VELY)
      F(HY_ZMOM) = F(HY_MASS)*V(HY_VELZ)
      F(HY_ENER) = (E + Ptot)*V(HY_VELX)
   case (JAXIS)
      F(HY_MASS) = V(HY_DENS)*V(HY_VELY)
      F(HY_XMOM) = F(HY_MASS)*V(HY_VELX)
      F(HY_YMOM) = F(HY_MASS)*V(HY_VELY) + Ptot
      F(HY_ZMOM) = F(HY_MASS)*V(HY_VELZ)
      F(HY_ENER) = (E + Ptot)*V(HY_VELY)
   case (KAXIS)
      F(HY_MASS) = V(HY_DENS)*V(HY_VELZ)
      F(HY_XMOM) = F(HY_MASS)*V(HY_VELX)
      F(HY_YMOM) = F(HY_MASS)*V(HY_VELY)
      F(HY_ZMOM) = F(HY_MASS)*V(HY_VELZ) + Ptot
      F(HY_ENER) = (E + Ptot)*V(HY_VELZ)
   end select

#endif
  
  end subroutine prim2flx
  