!!****if* source/physics/Hydro/HydroMain/unsplit/hy_getRiemannState
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!! NAME
!!
!!  hy_getRiemannState
!!
!! SYNOPSIS
!!
!!  hy_getRiemannState( integer(IN) :: blockID,
!!                          integer(IN) :: blkLimits,
!!                          integer(IN) :: blkLimitsGC,
!!                          real(IN)    :: dt,
!!                          integer(IN) :: del(MDIM),
!!                          real(IN)    :: ogravX(:,:,:),
!!                          real(IN)    :: ogravY(:,:,:),
!!                          real(IN)    :: ogravZ(:,:,:),
!!                          real, pointer, dimension (:,:,:,:) :: scrchFaceXPtr
!!                          real, pointer, dimension (:,:,:,:) :: scrchFaceYPtr
!!                          real, pointer, dimension (:,:,:,:) :: scrchFaceZPtr
!!                          real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig,
!!                          logical(IN), optional :: normalFieldUpdate)
!!
!! DESCRIPTION
!!
!!  This routine computes the Riemann state values at cell interfaces using
!!  the cell centered variables and store them in the scratch arrays.
!!
!!  A 2D Cartesian configuration of a single cell is shown:
!!
!!
!!           ---------------------
!!          |          yp         |
!!          |                     |
!!          |                     |
!!          |                     |
!!          |                     |
!!          |xm      (i,j)      xp|
!!          |                     |
!!          |                     |
!!          |                     |
!!  y       |                     |
!!  |       |          ym         |
!!  |        ---------------------
!!  |______x
!!
!!
!! ARGUMENTS
!!
!!  blockID     - local block ID
!!  blkLimits   - an array that holds the lower and upper indices of the section
!!                of block without the guard cells
!!  blkLimitsGC - an array that holds the lower and upper indices of the section
!!                of block with the guard cells
!!  dt          - a current time step
!!  del         - deltas in each direction
!!  ogravX       - gravity component in x-direction at n step
!!  ogravY       - gravity component in y-direction at n step
!!  ogravZ       - gravity component in z-direction at n step
!!  scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr - Pointers to the scrch array (for left/right states)
!!  hy_SpcR,hy_SpcL,hy_SpcSig - Pointers for Species and mass scalar recon.
!!  normalFieldUpdate - a logical switch to choose normal magnetic fields updates only
!!                      (needed for MHD only)
!!
!!***

!!REORDER(4):U, V0, scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr, B[xyz]

Subroutine hy_getRiemannState(tileDesc,U,blkLimits,loGC,hiGC,dt,del,&
                                  ogravX,ogravY,ogravZ,&
                                  scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr,&
                                  hy_SpcR,hy_SpcL,hy_SpcSig,&
                                  normalFieldUpdate)
#include "Simulation.h"
  use Hydro_data,           ONLY : hy_shockDetectOn,     &
                                   hy_meshMe,             &
                                   hy_order,             &
                                   hy_flattening,        &
                                   hy_useGravity,        &
                                   hy_useGravHalfUpdate, &
                                   hy_fallbackLowerCFL,  &
                                   hy_use3dFullCTU,      &
                                   hy_geometry,          &
                                   hy_numXN,             &
                                   hy_cfl,               &
                                   hy_cfl_original,      &
                                   hy_cflFallbackFactor, &
                                   hy_fullSpecMsFluxHandling, &
                                   hy_useHybridOrder,    &
                                   hy_eswitch,           &
                                   hy_upwindTVD

#ifdef FLASH_USM_MHD
  use Hydro_data,           ONLY : hy_killDivB,   &
                                   hy_forceHydroLimit
#endif
#ifdef FLASH_UGLM_MHD
  use Hydro_data,           ONLY : hy_C_hyp, hy_C_par
#endif

  use hy_slopeLimiters, ONLY : mc
  use hy_uhdInterface,  ONLY :     hy_dataReconstOnestep, &
                                   hy_shockDetect,        &
                                   hy_eigenParameters,    &
                                   hy_eigenValue,         &
                                   hy_eigenVector,        &
                                   hy_upwindTransverseFlux
  use Grid_interface,   ONLY : Grid_getCellCoords
  use Grid_tile,        ONLY : Grid_tile_t

  implicit none

#include "constants.h"
#include "UHD.h"

  !! Arguments type declaration ------------------------------------------------------------
  type(Grid_tile_t), intent(IN)   :: tileDesc
  integer, intent(IN),dimension(LOW:HIGH,MDIM):: blkLimits !, blkLimitsGC
  integer, intent(IN),dimension(MDIM)  :: loGC, hiGC
  real,    intent(IN)   :: dt
  real,    intent(IN),dimension(MDIM) :: del
  real, dimension(loGC(IAXIS):,loGC(JAXIS):,loGC(KAXIS):), intent(IN) :: ogravX,ogravY,ogravZ
  real, pointer, dimension(:,:,:,:) :: U
  real, pointer, dimension(:,:,:,:) :: scrchFaceXPtr, scrchFaceYPtr, scrchFaceZPtr
  real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig
  logical, intent(IN), optional :: normalFieldUpdate
  !! ---------------------------------------------------------------------------------------

  integer :: i0,imax,j0,jmax,k0,kmax, i, j, k
  integer,dimension(MDIM) :: dataSize
  real    :: cellCfl,minCfl
  logical :: lowerCflAtBdry
  integer :: dir
  integer,parameter :: dirLast=DIR_X+(NDIM-1)

! MHD only-------------------------------------------------------------------------------
#ifdef FLASH_USM_MHD
#if NFACE_VARS > 0
#if NDIM > 1
  real, pointer, dimension(:,:,:,:) :: Bx,By,Bz
#endif
#endif
  logical :: normalFieldUpdateOnly
#endif
! end of ifdef FLASH_USM_MHD
! MHD only-------------------------------------------------------------------------------

#ifdef FLASH_UHD_HYDRO
  logical, parameter :: normalFieldUpdateOnly = .FALSE.
#endif

  real, dimension(NDIM,loGC(IAXIS):hiGC(IAXIS),&
                       loGC(JAXIS):hiGC(JAXIS),&
                       loGC(KAXIS):hiGC(KAXIS)) :: FlatCoeff, FlatTilde
  real, dimension(     loGC(IAXIS):hiGC(IAXIS),&
                       loGC(JAXIS):hiGC(JAXIS),&
                       loGC(KAXIS):hiGC(KAXIS)) :: DivU

  real :: Sp, dv1, dp1, dp2, presL,presR,hdt

  real, dimension(loGC(IAXIS):hiGC(IAXIS)) :: xCenter
  real, dimension(loGC(JAXIS):hiGC(JAXIS)) :: yCenter

  integer :: k2,k3,kGrav,kHydro,kUSM,order
  integer :: k4,im2,ip2,jm2,jp2,km2,kp2

  ! cylindrical geometry
  integer :: velPhi, velTht, magPhi, magZ
  integer :: HY_velPhi, HY_velTht, H_magPhi, H_magZ

  real :: Rinv, geoFac, eta, enth
  real :: sGeo_dens, sGeo_velx, sGeo_velp, sGeo_velt, sGeo_pres
  real :: sGeo_trans, sGeo_eint, sGeo_magz,sGeo_magp

  ! for 3d only --------------
  real, dimension(HY_SPEC_END) :: TransFluxXY,TransFluxYZ,TransFluxZX,&
                                  TransFluxYX,TransFluxZY,TransFluxXZ
  real :: dt2dxdy6,dt2dydz6,dt2dzdx6,hdtdx,hdtdy,hdtdz
  logical :: cons=.false.
  real    :: cs,ca,cf,as,af,uN
  integer :: transOrder3D
  real, dimension(MDIM) :: beta
  ! for 3d only --------------
  logical :: TransX_updateOnly,TransY_updateOnly,TransZ_updateOnly
  logical :: allTransUpdateOnly

  integer,dimension(4) :: lbx,ubx,lby,uby,lbz,ubz
  integer :: iDim
  real, dimension(HY_VARINUMMAX,NDIM) :: Wp, Wn
  real, dimension(HY_END_VARS) :: Vc
  real, pointer, dimension(:)   :: SigmPtr,SigcPtr,SigpPtr
  !real, allocatable,dimension(:,:,:) :: DivU

  real, allocatable,dimension(:,:,:,:,:),target :: sig
  real, allocatable,dimension(:,:,:,:,:) :: lambda
  real, allocatable,dimension(:,:,:,:,:,:) :: leftEig
  real, allocatable,dimension(:,:,:,:,:,:) :: rghtEig

  ! GLM fluxes & state updates
#ifdef FLASH_UGLM_MHD
  real, dimension(loGC(IAXIS):hiGC(IAXIS),&
                  loGC(JAXIS):hiGC(JAXIS),&
                  loGC(KAXIS):hiGC(KAXIS)) :: &
       GLMxStar,GLMyStar,GLMzStar,BxStar,ByStar,BzStar
#endif /* FLASH_UGLM_MHD */




  !! Set transverse flux interpolation order
  transOrder3D = 1

  !! index for gravity
  kGrav = 0
#ifdef GRAVITY
  kGrav = 1
#endif

  !! index for pure Hydro
  kHydro = 1
#ifdef FLASH_USM_MHD
  if (.NOT. hy_forceHydroLimit) kHydro = 0

  normalFieldUpdateOnly = normalFieldUpdate ! the dummy argument is required to be present for MHD.
#endif
  kUSM = 1 - kHydro

  !! indices for various purposes
  k2=0
  k3=0
#if NDIM == 1
  i0   = @blk_ref(LOW, IAXIS)@
  imax = @blk_ref(HIGH,IAXIS)@
  j0   = 3
  jmax =-1
  k0   = 3
  kmax =-1
#elif NDIM == 2
  i0   = @blk_ref(LOW, IAXIS)@
  imax = @blk_ref(HIGH,IAXIS)@
  j0   = @blk_ref(LOW, JAXIS)@
  jmax = @blk_ref(HIGH,JAXIS)@
  k0   = 3
  kmax =-1
  k2=1
#elif NDIM == 3
  i0   = @blk_ref(LOW, IAXIS)@
  imax = @blk_ref(HIGH,IAXIS)@
  j0   = @blk_ref(LOW, JAXIS)@
  jmax = @blk_ref(HIGH,JAXIS)@
  k0   = @blk_ref(LOW, KAXIS)@
  kmax = @blk_ref(HIGH,KAXIS)@
  k2=1
  k3=1
#endif

  ! half delta t
  hdt = 0.5*dt

!!$#ifdef FLASH_UGLM_MHD
!!$  hy_C_hyp = 0.8/dt*min(del(DIR_X),del(DIR_Y))
!!$#endif

!!$  !Initialize geometric source terms
!!$  sGeo_dens=0.
!!$  sGeo_velx=0.
!!$  sGeo_velp=0.
!!$  sGeo_pres=0.
!!$  sGeo_trans=0.
!!$  sGeo_eint=0.
!!$  sGeo_magz=0.
!!$  sGeo_magp=0.


!!$  ! MHD only-------------------------------------------------------------------------------
!!$#if defined(FLASH_USM_MHD) && (NFACE_VARS > 0) && (NDIM > 1)
!!$  if (hy_order > 1) then
!!$     call Grid_getBlkPtr(blockID,Bx,FACEX)
!!$     call Grid_getBlkPtr(blockID,By,FACEY)
!!$     if (NDIM == 3) call Grid_getBlkPtr(blockID,Bz,FACEZ)
!!$  endif
!!$#endif /* endif of if defined(FLASH_USM_MHD) && NFACE_VARS > 0 && NDIM > 1 */
!!$  ! MHD only-------------------------------------------------------------------------------



  if (.NOT. normalFieldUpdateOnly) then

     @data_size@=@hg_rng@-@lg_rng@+1
     if (hy_geometry /= CARTESIAN) then
        ! Grab cell x-coords for this block
        call Grid_getCellCoords(IAXIS, CENTER, @tile_level@, &
                                loGC, hiGC, xCenter)
#if NDIM > 1
        if (hy_geometry == SPHERICAL) &
             call Grid_getCellCoords(JAXIS, CENTER, @tile_level@, &
                                     loGC, hiGC, yCenter)
#endif
     endif

     allocate( @alloc_sig@)
     allocate( @alloc_lambda@)
     allocate(@alloc_leftEig@)
     allocate(@alloc_rghtEig@)
  end if

  !! -----------------------------------------------------------------------!
  !! (1) Compute eigen structures once and for all -------------------------!
  !! (2) Compute divergence of velocity field      -------------------------!
  !! -----------------------------------------------------------------------!
  !if (hy_upwindTVD) kHydro = -1
!! Note --- fix all these indices
  if (hy_useHybridOrder .AND. .NOT. normalFieldUpdateOnly) then
     k4 = hy_order - 1             !cf. hy_dataReconstOnestep
     if (k4 > 2) k4 = 2 !(i.e., assume order = 3)
     k4 = min(NGUARD-2-kUSM,k4)
     im2=max(@spc1(loGC,IAXIS)@+1 ,i0  -2+kHydro-k4)
     ip2=min(@spc1(hiGC,IAXIS)@-1,imax+2-kHydro+k4)
#if NDIM > 1
     jm2=max(@spc1(loGC,JAXIS)@+1 ,j0  -2+kHydro-k4)
     jp2=min(@spc1(hiGC,JAXIS)@-1,jmax+2-kHydro+k4)
#else
     jm2=1; jp2=1
#endif
#if NDIM > 2
     km2=max(@spc1(loGC,KAXIS)@+1 ,k0  -2+kHydro-k4)
     kp2=min(@spc1(hiGC,KAXIS)@-1,kmax+2-kHydro+k4)
#else
     km2=1; kp2=1
#endif
     do k=km2,kp2
        do j=jm2,jp2
           do i=im2,ip2
              !! We used to compute undivided divergence of velocity fields and (magneto)sonic speed
              !! and store local (magneto)sonic speeds for hybrid order. Now the latter are
              !! computed elsewhere.
              @spc3(DivU,+0,+0,+0)@ = @fld(U,VELX_VAR,+1,+0,+0)@-@fld(U,VELX_VAR,-1,+0,+0)@
              if (NDIM > 1) then
                 @spc3(DivU,+0,+0,+0)@ = @spc3(DivU,+0,+0,+0)@ \
                      +@fld(U,VELY_VAR,+0,+1,+0)@-@fld(U,VELY_VAR,+0,-1,+0)@
                 if (NDIM > 2) then
                    @spc3(DivU,+0,+0,+0)@ = @spc3(DivU,+0,+0,+0)@ \
                         +@fld(U,VELZ_VAR,+0,+0,+1)@-@fld(U,VELZ_VAR,+0,+0,-1)@
                 endif
              endif
              @spc3(DivU,+0,+0,+0)@ = 0.5*@spc3(DivU,+0,+0,+0)@

           enddo ! do i-loop
        enddo ! do j-loop
     enddo ! do k-loop
  endif !end of hybridOrder


  !! -----------------------------------------------------------------------!
  !! (3) Flattening begins here --------------------------------------------!
  !! -----------------------------------------------------------------------!
  if (hy_flattening .AND. .NOT. normalFieldUpdateOnly) then

     ! Initialize with zero
     FlatTilde = 0.
     FlatCoeff = 0.

     ! Flat tilde
     do k=k0-2-min(NGUARD-4,kUSM)*k3,kmax+2+min(NGUARD-4,kUSM)*k3
        do j=j0-2-min(NGUARD-4,kUSM)*k2,jmax+2+min(NGUARD-4,kUSM)*k2
           do i=i0-2-min(NGUARD-4,kUSM),imax+2+min(NGUARD-4,kUSM)
              do dir=1,NDIM
                 select case (dir)
                 case (DIR_X)
                    dp1   = (@fld(U,PRES_VAR,+1,+0,+0)@-@fld(U,PRES_VAR,-1,+0,+0)@)
                    dp2   = (@fld(U,PRES_VAR,+2,+0,+0)@-@fld(U,PRES_VAR,-2,+0,+0)@)
                    dv1   =  @fld(U,VELX_VAR,+1,+0,+0)@-@fld(U,VELX_VAR,-1,+0,+0)@
                    presL =  @fld(U,PRES_VAR,-1,+0,+0)@
                    presR =  @fld(U,PRES_VAR,+1,+0,+0)@
#if NDIM > 1
                 case (DIR_Y)
                    dp1   = (@fld(U,PRES_VAR,+0,+1,+0)@-@fld(U,PRES_VAR,+0,-1,+0)@)
                    dp2   = (@fld(U,PRES_VAR,+0,+2,+0)@-@fld(U,PRES_VAR,+0,-2,+0)@)
                    dv1   =  @fld(U,VELY_VAR,+0,+1,+0)@-@fld(U,VELY_VAR,+0,-1,+0)@
                    presL =  @fld(U,PRES_VAR,+0,-1,+0)@
                    presR =  @fld(U,PRES_VAR,+0,+1,+0)@
#if NDIM > 2
                 case (DIR_Z)
                    dp1   = (@fld(U,PRES_VAR,+0,+0,+1)@-@fld(U,PRES_VAR,+0,+0,-1)@)
                    dp2   = (@fld(U,PRES_VAR,+0,+0,+2)@-@fld(U,PRES_VAR,+0,+0,-2)@)
                    dv1   =  @fld(U,VELZ_VAR,+0,+0,+1)@-@fld(U,VELZ_VAR,+0,+0,-1)@
                    presL =  @fld(U,PRES_VAR,+0,+0,-1)@
                    presR =  @fld(U,PRES_VAR,+0,+0,+1)@
#endif
#endif
                 end select

                 if (abs(dp2) > 1.e-15) then
                    Sp = dp1/dp2 - 0.75
                 else
                    Sp = 0.
                 endif

                 @spc3d(FlatTilde,dir,+0,+0,+0)@ = max(0.0, min(1.0,10.0*Sp))
                 if ((abs(dp1)/min(presL,presR) < 1./3.) .or. dv1 > 0. ) then
                    @spc3d(FlatTilde,dir,+0,+0,+0)@ = 0.
                 endif

              enddo
           enddo
        enddo
     enddo

     ! Flat coefficient
     do k=k0-2+kHydro*k3,kmax+2-kHydro*k3
        do j=j0-2+kHydro*k2,jmax+2-kHydro*k2
           do i=i0-2+kHydro,imax+2-kHydro
              do dir=1,NDIM

                 select case (dir)
                 case (DIR_X)
                    dp1   = (@fld(U,PRES_VAR,+1,+0,+0)@-@fld(U,PRES_VAR,-1,+0,+0)@)

                    if ( dp1 < 0.0 ) then
                       @spc3d(FlatCoeff,dir,+0,+0,+0)@ = max(@spc3d(FlatTilde,dir,+0,+0,+0)@,@spc3d(FlatTilde,dir,+1,+0,+0)@)
                    elseif (dp1 == 0.) then
                       @spc3d(FlatCoeff,dir,+0,+0,+0)@ = @spc3d(FlatTilde,dir,+0,+0,+0)@
                    else
                       @spc3d(FlatCoeff,dir,+0,+0,+0)@ = max(@spc3d(FlatTilde,dir,+0,+0,+0)@,@spc3d(FlatTilde,dir,-1,+0,+0)@)
                    endif
#if NDIM > 1
                 case (DIR_Y)
                    dp1   = (@fld(U,PRES_VAR,+0,+1,+0)@-@fld(U,PRES_VAR,+0,-1,+0)@)

                    if ( dp1 < 0.0 ) then
                       @spc3d(FlatCoeff,dir,+0,+0,+0)@ = max(@spc3d(FlatTilde,dir,+0,+0,+0)@,@spc3d(FlatTilde,dir,+0,+1,+0)@)
                    elseif (dp1 == 0.) then
                       @spc3d(FlatCoeff,dir,+0,+0,+0)@ = @spc3d(FlatTilde,dir,+0,+0,+0)@
                    else
                       @spc3d(FlatCoeff,dir,+0,+0,+0)@ = max(@spc3d(FlatTilde,dir,+0,+0,+0)@,@spc3d(FlatTilde,dir,+0,-1,+0)@)
                    endif
#if NDIM > 2
                 case (DIR_Z)
                    dp1   = (@fld(U,PRES_VAR,+0,+0,+1)@-@fld(U,PRES_VAR,+0,+0,-1)@)

                    if ( dp1 < 0.0 ) then
                       @spc3d(FlatCoeff,dir,+0,+0,+0)@ = max(@spc3d(FlatTilde,dir,+0,+0,+0)@,@spc3d(FlatTilde,dir,+0,+0,+1)@)
                    elseif (dp1 == 0.) then
                       @spc3d(FlatCoeff,dir,+0,+0,+0)@ = @spc3d(FlatTilde,dir,+0,+0,+0)@
                    else
                       @spc3d(FlatCoeff,dir,+0,+0,+0)@ = max(@spc3d(FlatTilde,dir,+0,+0,+0)@,@spc3d(FlatTilde,dir,+0,+0,-1)@)
                    endif
#endif
#endif
                 end select
              enddo
           enddo
        enddo
     enddo
  endif

  !! -----------------------------------------------------------------------!
  !! (4) Start calculating Riemann states ----------------------------------!
  !! -----------------------------------------------------------------------!
  !! Compute Riemann states at each cell
  if (.not. normalFieldUpdateOnly) then

#ifdef CFL_VAR
     minCfl = hy_cfl_original
#else
     minCfl = hy_cfl
#endif

!!$     print*,'RiemannSt lbound(U):',lbound(U)
!!$     print*,'RiemannSt ubound(U):',ubound(U)

     do k=k0-2-(k3*kUSM-kHydro)*k3,kmax+2+(k3*kUSM-kHydro)*k3

        do j=j0-2-(k3*kUSM-kHydro)*k2,jmax+2+(k3*kUSM-kHydro)*k2
           do i=i0-2-(k3*kUSM-kHydro),imax+2+(k3*kUSM-kHydro)
           ! Extra stencil is needed for 3D to correctly calculate transverse fluxes
           !(i.e., cross derivatives in x,y, & z)
!!$              print*,'RiemannSt loop:',i,j,' ...'
              !! save the cell center values for later use
              @vc_rng(HY_DENS,HY_END_VARS-kGrav)@ = &
                      (/@fld(U,DENS_VAR,+0,+0,+0)@&
                       ,@fld_rng(U,VELX_VAR,VELZ_VAR,+0,+0,+0)@&
                       ,@fld(U,PRES_VAR,+0,+0,+0)@&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                       ,@fld_rng(U,MAGX_VAR,MAGZ_VAR,+0,+0,+0)@&
#endif
#ifdef FLASH_UGLM_MHD
                       ,@fld(U,GLMP_VAR,+0,+0,+0)@ &
#endif
                       ,@fld(U,GAMC_VAR,+0,+0,+0)@ &
                       ,@fld(U,GAME_VAR,+0,+0,+0)@ &
                       ,@fld(U,EINT_VAR,+0,+0,+0)@ &
#ifdef FLASH_UHD_3T
                       ,@fld(U,EELE_VAR,+0,+0,+0)@ &
                       ,@fld(U,EION_VAR,+0,+0,+0)@ &
                       ,@fld(U,ERAD_VAR,+0,+0,+0)@ &
#endif
                       /)

              order = hy_order
              lowerCflAtBdry = .FALSE.
#ifdef BDRY_VAR
              !! Reduce order in fluid cell near solid boundary if defined:
              !! Reduce order of spatial reconstruction depending on the distance to the solid boundary
              if (order > 2) then
!!!!! fix this
                 im2=max(@spc1(loGC,IAXIS)@,i-2); ip2=min(@spc1(hiGC,IAXIS)@,i+2)
#if NDIM > 1
                 jm2=max(@spc1(loGC,JAXIS)@,j-2); jp2=min(@spc1(hiGC,JAXIS)@,j+2)
#else
                 jm2 = 1; jp2=1
#endif
#if NDIM > 2
                 km2=max(@spc1(loGC,KAXIS)@,k-2); kp2=min(@spc1(hiGC,KAXIS)@,k+2)
#else
                 km2 = 1; kp2=1
#endif
                 if (maxval(@fld_slice_v(U,BDRY_VAR,im2,ip2,jm2,jp2,km2,kp2)@) .LE. 0.) then !everyone is fluid
                    order = 3
                 else
                    order = 2
                 endif
              endif
              (@fld(U,BDRY_VAR,+0,+0,+0)@*@fld(U,BDRY_VAR,-1,+0,+0)@ < 0.0) .or. &
              (@fld(U,BDRY_VAR,+0,+0,+0)@*@fld(U,BDRY_VAR,+1,+0,+0)@ < 0.0) .or. &
              (@fld(U,BDRY_VAR,+0,+0,+0)@*@fld(U,BDRY_VAR,+0,-k2,+0)@ < 0.0) .or. &
              (@fld(U,BDRY_VAR,+0,+0,+0)@*@fld(U,BDRY_VAR,+0,+k2,+0)@ < 0.0) .or. &
              (@fld(U,BDRY_VAR,+0,+0,+0)@*@fld(U,BDRY_VAR,+0,+0,-k3)@ < 0.0) .or. &
              (@fld(U,BDRY_VAR,+0,+0,+0)@*@fld(U,BDRY_VAR,+0,+0,+k3)@ < 0.0)) then
                 order = 1
                 lowerCflAtBdry = .TRUE.
              endif
#if NDIM > 2
              ! Addtional 3D test whether a solid cell is so close
              ! that we cannot do all the proper transverse computations for this cell.
              if (hy_use3dFullCTU) then
                 if (maxval(@fld_slice(U,BDRY_VAR,-1,+1,-1,+1,-1,+1)@) > 0.) then
                    if (maxval(@fld_slice(U,BDRY_VAR,+0,+0,-1,+1,-1,+1)@) > 0.) lowerCflAtBdry = .TRUE.
                    if (maxval(@fld_slice(U,BDRY_VAR,-1,+1,+0,+0,-1,+1)@) > 0.) lowerCflAtBdry = .TRUE.
                    if (maxval(@fld_slice(U,BDRY_VAR,-1,+1,-1,+1,+0,+0)@) > 0.) lowerCflAtBdry = .TRUE.
                 end if
              end if
#endif
#endif

#ifdef CFL_VAR
              cellCfl = @fld(U,CFL_VAR,+0,+0,+0)@
#else
              cellCfl = hy_cfl
#endif
              if (lowerCflAtBdry) cellCfl = min(cellCfl, hy_cflFallbackFactor / real(NDIM))

              !! Flag for tranverse update
              TransX_updateOnly = .false.
              TransY_updateOnly = .false.
              TransZ_updateOnly = .false.

              if (i > @blk_ref(HIGH,IAXIS)@+1 .or. &
                   i < @blk_ref(LOW, IAXIS)@-1) then
                 TransX_updateOnly = .true.
              endif
              if (i > @blk_ref(HIGH,IAXIS)@+2*kUSM .or. &
                   i < @blk_ref(LOW, IAXIS)@-2*kUSM) then
                 TransY_updateOnly = .true.
                 TransZ_updateOnly = .true.
              endif
#if NDIM > 1
              if (j > @blk_ref(HIGH,JAXIS)@+1 .or. &
                   j < @blk_ref(LOW, JAXIS)@-1) then
                 TransY_updateOnly = .true.
              endif
              if (j > @blk_ref(HIGH,JAXIS)@+2*kUSM .or. &
                   j < @blk_ref(LOW, JAXIS)@-2*kUSM) then
                 TransX_updateOnly = .true.
                 TransZ_updateOnly = .true.
              endif

#if NDIM > 2
              if (k > @blk_ref(HIGH,KAXIS)@+1 .or. &
                   k < @blk_ref(LOW, KAXIS)@-1) then
                 TransZ_updateOnly = .true.
              endif
              if (k > @blk_ref(HIGH,KAXIS)@+2*kUSM .or. &
                   k < @blk_ref(LOW, KAXIS)@-2*kUSM) then
                 TransX_updateOnly = .true.
                 TransY_updateOnly = .true.
              endif
#endif
#endif
              allTransUpdateOnly = TransX_updateOnly
              if (NDIM > 1) allTransUpdateOnly = allTransUpdateOnly .AND. TransY_updateOnly
              if (NDIM > 2) allTransUpdateOnly = allTransUpdateOnly .AND. TransZ_updateOnly

              if (order == 1) then
                 !! DEV: THE FIRST ORDER SHOULD GO INTO THE DATA RECONSTRUCT ONE STEP TO GET
                 !! TRANSVERSE FLUXES
                 if (.NOT. allTransUpdateOnly) then
                    do iDim = 1,NDIM
                       call fallbackToFirstOrder(iDim,@wn_section(iDim)@,@wp_section(iDim)@,Vc,hy_SpcL,hy_SpcR,U,i,j,k)
                    enddo
                 end if

              else ! for high-order schemes


                 !! Left and right Riemann state reconstructions
                 if (hy_fullSpecMsFluxHandling .AND. hy_numXN > 0 &
                      .AND. present(hy_spcR)) then
                    call hy_dataReconstOnestep&
                      (tileDesc,U,loGC,hiGC,    &
                       order,i,j,k,dt,del,     &
                       ogravX,ogravY,ogravZ,   &
                       DivU,FlatCoeff,         &
                       TransX_updateOnly,      &
                       TransY_updateOnly,      &
                       TransZ_updateOnly,      &
                       Wp,Wn,                  &
                       @sig_arg(sig,+0,+0,+0)@,   &
                       @lam_arg(lambda,+0,+0,+0)@,   &
                       @leig_arg(leftEig,+0,+0,+0)@,   &
                       @reig_arg(rghtEig,+0,+0,+0)@,   &
                       cellCfl, &
                       hy_SpcR,hy_SpcL,hy_SpcSig)
                 else
                    call hy_dataReconstOnestep&
                      (tileDesc,U,loGC,hiGC,    &
                       order,i,j,k,dt,del,     &
                       ogravX,ogravY,ogravZ,   &
                       DivU,FlatCoeff,         &
                       TransX_updateOnly,      &
                       TransY_updateOnly,      &
                       TransZ_updateOnly,      &
                       Wp,Wn,                  &
                       @sig_arg(sig,+0,+0,+0)@,   &
                       @lam_arg(lambda,+0,+0,+0)@,   &
                       @leig_arg(leftEig,+0,+0,+0)@,   &
                       @reig_arg(rghtEig,+0,+0,+0)@,   &
                       cellCfl)
                 endif ! if(hy_fullSpecMsFluxHandling ...

              endif ! end of high-order reconstruction schemes


              if (hy_geometry /= CARTESIAN) then
                 !! **************************************************************
                 !! Add geometric source terms in left and Right States          *
                 !! **************************************************************


                 !Initialize geometric source terms
                 sGeo_dens = 0.
                 sGeo_velx = 0.
                 sGeo_velp = 0.
                 sGeo_pres = 0.
                 sGeo_eint = 0.
                 sGeo_magz = 0.
                 sGeo_magp = 0.
                 sGeo_trans= 0.

                 Rinv = 1./@spc1(xCenter,i)@
                 select case (hy_geometry)
                 case (CYLINDRICAL)
                    velPhi    = VELZ_VAR
                    HY_velPhi = HY_VELZ
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                    magPhi    = MAGZ_VAR
                    magZ      = MAGY_VAR
                    H_magPhi  = HY_MAGZ
                    H_magZ    = HY_MAGY
#endif
                    geoFac    = Rinv
                 case (POLAR)
                    velPhi    = VELY_VAR
                    HY_velPhi = HY_VELY
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                    magPhi    = MAGY_VAR
                    magZ      = MAGZ_VAR
                    H_magPhi  = HY_MAGY
                    H_magZ    = HY_MAGZ
#endif
                    geoFac    = Rinv
                 case (SPHERICAL)
                    velPhi    = VELZ_VAR
                    velTht    = VELY_VAR
                    HY_velPhi = HY_VELZ
                    HY_velTht = HY_VELY
                    geoFac    = 2.*Rinv
                 end select

                 cs  = sqrt(@fld(U,GAMC_VAR,+0,+0,+0)@*@fld(U,PRES_VAR,+0,+0,+0)@/@fld(U,DENS_VAR,+0,+0,+0)@)
                 eta = (abs(@fld(U,VELX_VAR,+0,+0,+0)@) + cs) * dt/@del_ref(DIR_X)@
                 eta = (1.-eta) / (cs*dt*abs(geoFac))
                 eta = min(1.,eta)
                 !! comment this line not to use the axis hack
!!$                 geoFac = eta * geoFac
                 !! end of the axis hack
                 enth = @fld(U,EINT_VAR,+0,+0,+0)@ + @fld(U,PRES_VAR,+0,+0,+0)@/@fld(U,DENS_VAR,+0,+0,+0)@

                 !! right/left state source terms
                 sGeo_dens = -@fld(U,DENS_VAR,+0,+0,+0)@ * @fld(U,VELX_VAR,+0,+0,+0)@ * geoFac !src[DN]
#if NDIM > 1
                 if (hy_geometry == SPHERICAL) then
                      sGeo_dens = sGeo_dens &
                      -@fld(U,DENS_VAR,+0,+0,+0)@*@fld(U,velTht,+0,+0,+0)@*cos(@spc1(yCenter,j)@)/sin(@spc1(yCenter,j)@) * 0.5*geoFac
                      sGeo_velt = (@fld(U,velPhi,+0,+0,+0)@**2) * cos(@spc1(yCenter,j)@)/sin(@spc1(yCenter,j)@) * 0.5*geoFac
                   end if
#endif
                 sGeo_velx = (@fld(U,velPhi,+0,+0,+0)@**2) * geoFac                   !src[VR]
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                 sGeo_velx = sGeo_velx - (@fld(U,magPhi,+0,+0,+0)@**2) * geoFac / @fld(U,DENS_VAR,+0,+0,+0)@
                 if (hy_geometry == SPHERICAL) &
                      sGeo_velx = sGeo_velx + @fld(U,velTht,+0,+0,+0)@**2 * geoFac
#endif
                 sGeo_velp = -@fld(U,velPhi,+0,+0,+0)@ * @fld(U,VELX_VAR,+0,+0,+0)@ * geoFac !src[Vphi]
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                 sGeo_velp = sGeo_velp + @fld(U,magPhi,+0,+0,+0)@ * @fld(U,MAGX_VAR,+0,+0,+0)@ * geoFac / @fld(U,DENS_VAR,+0,+0,+0)@
#endif
                 sGeo_pres = sGeo_dens * cs**2                              !src[PR]
                 sGeo_eint = (sGeo_dens*enth)/@fld(U,DENS_VAR,+0,+0,+0)@
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                 sGeo_magp = - @fld(U,velPhi,+0,+0,+0)@ * @fld(U,MAGX_VAR,+0,+0,+0)@ * geoFac
                 sGeo_magz = - @fld(U,VELX_VAR,+0,+0,+0)@ * @fld(U,magZ,+0,+0,+0)@ * geoFac
#endif

#if (0)
#define DIRS DIR_X
#else
#define DIRS DIR_X:dirLast
#endif

                 !! Add sources terms for n+1/2 Left state
                 @wn_ref(HY_DENS,  DIRS)@ = @wn_ref(HY_DENS,  DIRS)@ + hdt * sGeo_dens
                 @wn_ref(HY_VELX,  DIR_X)@ = @wn_ref(HY_VELX,  DIR_X)@ + hdt * sGeo_velx
                 @wn_ref(HY_velPhi,DIR_X)@ = @wn_ref(HY_velPhi,DIR_X)@ + hdt * sGeo_velp
                 @wn_ref(HY_PRES,  DIRS)@ = @wn_ref(HY_PRES,  DIRS)@ + hdt * sGeo_pres
                 @wn_ref(HY_EINT,  DIRS)@ = @wn_ref(HY_EINT,  DIRS)@ + hdt * sGeo_eint
#if NDIM > 1
                 if (hy_geometry == SPHERICAL) &
                      @wn_ref(HY_velTht,DIRS)@ = @wn_ref(HY_velTht,DIRS)@ + hdt * sGeo_velt
#endif
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                 @wn_ref(H_magPhi, DIR_X)@ = @wn_ref(H_magPhi, DIR_X)@ + hdt * sGeo_magp
                 @wn_ref(H_magZ,   DIR_X)@ = @wn_ref(H_magZ,   DIR_X)@ + hdt * sGeo_magz
#endif
                 !! Add source terms for n+1/2 Right state
                 @wp_ref(HY_DENS,  DIRS)@ = @wp_ref(HY_DENS,  DIRS)@ + hdt * sGeo_dens
                 @wp_ref(HY_VELX,  DIR_X)@ = @wp_ref(HY_VELX,  DIR_X)@ + hdt * sGeo_velx
                 @wp_ref(HY_velPhi,DIR_X)@ = @wp_ref(HY_velPhi,DIR_X)@ + hdt * sGeo_velp
                 @wp_ref(HY_PRES,  DIRS)@ = @wp_ref(HY_PRES,  DIRS)@ + hdt * sGeo_pres
                 @wp_ref(HY_EINT,  DIRS)@ = @wp_ref(HY_EINT,  DIRS)@ + hdt * sGeo_eint
#if NDIM > 1
                 if (hy_geometry == SPHERICAL) &
                      @wp_ref(HY_velTht,DIRS)@ = @wp_ref(HY_velTht,DIRS)@ + hdt * sGeo_velt
#endif
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                 @wp_ref(H_magPhi, DIR_X)@ = @wp_ref(H_magPhi, DIR_X)@ + hdt * sGeo_magp
                 @wp_ref(H_magZ,   DIR_X)@ = @wp_ref(H_magZ,   DIR_X)@ + hdt * sGeo_magz
#endif


                 if (@spc1(xCenter,i)@ - 0.5*@del_ref(DIR_X)@ == 0.) then
                    ! the velocity should be zero at r=0.
                    @wn_ref(HY_VELX,  DIR_X)@ = 0.0
                    @wn_ref(HY_velPhi,DIR_X)@ = 0.0
                    if (hy_geometry == SPHERICAL) @wn_ref(HY_velTht,DIR_X)@ = 0.0
#ifdef FLASH_USM_MHD || FLASH_UGLM_MHD
                    @wn_ref(HY_MAGX,  DIR_X)@ = 0.0
                    @wn_ref(H_magPhi, DIR_X)@ = 0.0
#endif
                 elseif (@spc1(xCenter,i)@ + 0.5*@del_ref(DIR_X)@ == 0.) then
                    @wp_ref(HY_VELX,  DIR_X)@ = 0.0
                    @wp_ref(HY_velPhi,DIR_X)@ = 0.0
                    if (hy_geometry == SPHERICAL) @wp_ref(HY_velTht,DIR_X)@ = 0.0
#ifdef FLASH_USM_MHD) || FLASH_UGLM_MHD
                    @wp_ref(HY_MAGX,  DIR_X)@ = 0.0
                    @wp_ref(H_magPhi, DIR_X)@ = 0.0
#endif
                 endif
                 !! Calculate R-momentum geometric source term for transverse fluxes
                 !! We will use the cell-centered states at t^n
                 sGeo_trans = (@fld(U,velPhi,+0,+0,+0)@**2) * geoFac
#ifdef FLASH_USM_MHD || FLASH_UGLM_MHD
                 sGeo_trans = sGeo_trans - (@fld(U,magPhi,+0,+0,+0)@**2) * geoFac/@fld(U,DENS_VAR,+0,+0,+0)@
#endif
                 if (hy_geometry == SPHERICAL) &
                      sGeo_trans = sGeo_trans + @fld(U,velTht,+0,+0,+0)@**2 * geoFac
#if NDIM > 1
                 @wn_ref(HY_VELX,DIR_Y)@ = @wn_ref(HY_VELX,DIR_Y)@ + hdt * sGeo_trans
                 @wp_ref(HY_VELX,DIR_Y)@ = @wp_ref(HY_VELX,DIR_Y)@ + hdt * sGeo_trans
#if NDIM == 3
                 @wn_ref(HY_VELX,DIR_Z)@ = @wn_ref(HY_VELX,DIR_Z)@ + hdt * sGeo_trans
                 @wp_ref(HY_VELX,DIR_Z)@ = @wp_ref(HY_VELX,DIR_Z)@ + hdt * sGeo_trans
#endif
#endif
                 !! Check positivity of density and pressure
                 if (@wn_ref(HY_DENS,DIR_X)@ < 0. .or. @wp_ref(HY_DENS,DIR_X)@ < 0. .or. &
                     @wn_ref(HY_PRES,DIR_X)@ < 0. .or. @wp_ref(HY_PRES,DIR_X)@ < 0.) then
                    if(.NOT.TransX_updateOnly) then
                       print*,'[gRSt] afterGeo fallback to order 1 for DIR_X at i,j=',i,j,&
                            ' in Block ',hy_meshMe
                         print*,'[gRSt] afterGeo',@wn_ref(HY_DENS,DIR_X)@,@wp_ref(HY_DENS,DIR_X)@, &
                              @wn_ref(HY_PRES,DIR_X)@, @wp_ref(HY_PRES,DIR_X)@
                           call fallbackToFirstOrder(DIR_X,@wn_section(DIR_X)@,@wp_section(DIR_X)@,Vc,hy_SpcL,hy_SpcR,U,i,j,k)
                      end if
                   end if
#if NDIM > 1
                   if (@wn_ref(HY_DENS,DIR_Y)@ < 0. .or. @wp_ref(HY_DENS,DIR_Y)@ < 0. .or. &
                       @wn_ref(HY_PRES,DIR_Y)@ < 0. .or. @wp_ref(HY_PRES,DIR_Y)@ < 0.) then
                      if(.NOT.TransY_updateOnly) then
                         print*,'[gRSt] afterGeo fallback to order 1 for DIR_Y at i,j=',i,j,&
                              ' in Block ',hy_meshMe
                         print*,'[gRSt] afterGeo',@wn_ref(HY_DENS,DIR_Y)@,@wp_ref(HY_DENS,DIR_Y)@, &
                              @wn_ref(HY_PRES,DIR_Y)@, @wp_ref(HY_PRES,DIR_Y)@
                           call fallbackToFirstOrder(DIR_Y,@wn_section(DIR_Y)@,@wp_section(DIR_Y)@,Vc,hy_SpcL,hy_SpcR,U,i,j,k)
                      end if
                   end if
#if NDIM ==3
                   if (@wn_ref(HY_DENS,DIR_Z)@ < 0. .or. @wp_ref(HY_DENS,DIR_Z)@ < 0. .or. &
                       @wn_ref(HY_PRES,DIR_Z)@ < 0. .or. @wp_ref(HY_PRES,DIR_Z)@ < 0.) then
                      if(.NOT.TransZ_updateOnly) &
                           call fallbackToFirstOrder(DIR_Z,@wn_section(DIR_Z)@,@wp_section(DIR_Z)@,Vc,hy_SpcL,hy_SpcR,U,i,j,k)
                   end if
#endif
#endif
              endif !end if of if (hy_geometry .ne. CARTESIAN)

#ifdef CFL_VAR
              @fld(U,CFL_VAR,+0,+0,+0)@ = cellCfl
#endif
              minCfl = min(minCfl,cellCfl)


#ifdef GRAVITY
              if (hy_useGravity .and. hy_useGravHalfUpdate)then
                 @wp_rng_ref(HY_VELX,HY_VELZ,DIR_X)@=@wp_rng_ref(HY_VELX,HY_VELZ,DIR_X)@&
                      +hdt*(/@wp_ref(HY_GRAV,DIR_X)@,@spc3(ogravY,+0,+0,+0)@,@spc3(ogravZ,+0,+0,+0)@/)

                 @wn_rng_ref(HY_VELX,HY_VELZ,DIR_X)@=@wn_rng_ref(HY_VELX,HY_VELZ,DIR_X)@&
                      +hdt*(/@wn_ref(HY_GRAV,DIR_X)@,@spc3(ogravY,+0,+0,+0)@,@spc3(ogravZ,+0,+0,+0)@/)
              endif
#endif
              !! Store Riemann states to scratch arrays

              if(.NOT.TransX_updateOnly) @fld_rng(scrchFaceXPtr,HY_P01_FACEXPTR_VAR,HY_P01_FACEXPTR_VAR+HY_SCRATCH_NUM-1,+0,+0,+0)@&
                   = @wp_rng_ref(HY_DENS,HY_END_VARS-kGrav,DIR_X)@

              if(.NOT.TransX_updateOnly) @fld_rng(scrchFaceXPtr,HY_N01_FACEXPTR_VAR,HY_N01_FACEXPTR_VAR+HY_SCRATCH_NUM-1,+0,+0,+0)@&
                   = @wn_rng_ref(HY_DENS,HY_END_VARS-kGrav,DIR_X)@

#if NDIM >= 2
#ifdef GRAVITY
              if (hy_useGravity .and. hy_useGravHalfUpdate)then
                 @wp_rng_ref(HY_VELX,HY_VELZ,DIR_Y)@=@wp_rng_ref(HY_VELX,HY_VELZ,DIR_Y)@&
                      +hdt*(/@spc3(ogravX,+0,+0,+0)@,@wp_ref(HY_GRAV,DIR_Y)@,@spc3(ogravZ,+0,+0,+0)@/)

                 @wn_rng_ref(HY_VELX,HY_VELZ,DIR_Y)@=@wn_rng_ref(HY_VELX,HY_VELZ,DIR_Y)@&
                      +hdt*(/@spc3(ogravX,+0,+0,+0)@,@wn_ref(HY_GRAV,DIR_Y)@,@spc3(ogravZ,+0,+0,+0)@/)
              endif
#endif
              !! Store Riemann states to scratch arrays
              if(.NOT.TransY_updateOnly) @fld_rng(scrchFaceYPtr,HY_P01_FACEYPTR_VAR,HY_P01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,+0,+0,+0)@&
                   = @wp_rng_ref(HY_DENS,HY_END_VARS-kGrav,DIR_Y)@

              if(.NOT.TransY_updateOnly) @fld_rng(scrchFaceYPtr,HY_N01_FACEYPTR_VAR,HY_N01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,+0,+0,+0)@&
                   = @wn_rng_ref(HY_DENS,HY_END_VARS-kGrav,DIR_Y)@
#if NDIM == 3
#ifdef GRAVITY
              if (hy_useGravity .and. hy_useGravHalfUpdate)then
                 @wp_rng_ref(HY_VELX,HY_VELZ,DIR_Z)@=@wp_rng_ref(HY_VELX,HY_VELZ,DIR_Z)@&
                      +hdt*(/@spc3(ogravX,+0,+0,+0)@,@spc3(ogravY,+0,+0,+0)@,@wp_ref(HY_GRAV,DIR_Z)@/)

                 @wn_rng_ref(HY_VELX,HY_VELZ,DIR_Z)@=@wn_rng_ref(HY_VELX,HY_VELZ,DIR_Z)@&
                      +hdt*(/@spc3(ogravX,+0,+0,+0)@,@spc3(ogravY,+0,+0,+0)@,@wn_ref(HY_GRAV,DIR_Z)@/)
              endif
#endif
              !! Store Riemann states to scratch arrays
              if(.NOT.TransZ_updateOnly) @fld_rng(scrchFaceZPtr,HY_P01_FACEZPTR_VAR,HY_P01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,+0,+0,+0)@&
                   = @wp_rng_ref(HY_DENS,HY_END_VARS-kGrav,DIR_Z)@

              if(.NOT.TransZ_updateOnly) @fld_rng(scrchFaceZPtr,HY_N01_FACEZPTR_VAR,HY_N01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,+0,+0,+0)@&
                   = @wn_rng_ref(HY_DENS,HY_END_VARS-kGrav,DIR_Z)@
#endif
#endif

#if NFACE_VARS > 0
#if NDIM > 1
              if (hy_order > 1 .and. hy_killdivB .and. (.not. hy_forceHydroLimit)) then
                 if(.NOT.TransX_updateOnly) @fld(scrchFaceXPtr,HY_P06_FACEXPTR_VAR,+0,+0,+0)@= @fld(Bx,MAG_FACE_VAR,+1,+0,+0)@
                 if(.NOT.TransX_updateOnly) @fld(scrchFaceXPtr,HY_N06_FACEXPTR_VAR,+0,+0,+0)@= @fld(Bx,MAG_FACE_VAR,+0,+0,+0)@

                 if(.NOT.TransY_updateOnly) @fld(scrchFaceYPtr,HY_P07_FACEYPTR_VAR,+0,+0,+0)@= @fld(By,MAG_FACE_VAR,+0,+1,+0)@
                 if(.NOT.TransY_updateOnly) @fld(scrchFaceYPtr,HY_N07_FACEYPTR_VAR,+0,+0,+0)@= @fld(By,MAG_FACE_VAR,+0,+0,+0)@
#if NDIM == 3
                 if(.NOT.TransZ_updateOnly) @fld(scrchFaceZPtr,HY_P08_FACEZPTR_VAR,+0,+0,+0)@= @fld(Bz,MAG_FACE_VAR,+0,+0,+1)@
                 if(.NOT.TransZ_updateOnly) @fld(scrchFaceZPtr,HY_N08_FACEZPTR_VAR,+0,+0,+0)@= @fld(Bz,MAG_FACE_VAR,+0,+0,+0)@
#endif
              endif
#endif
#endif
           enddo
        enddo
     enddo

! MHD only-------------------------------------------------------------------------------
#ifdef FLASH_USM_MHD
  else ! else of if (.not. normalFieldUpdateOnly) then
       !! Updates of magnetic fields in normal direction
#if NFACE_VARS > 0
#if NDIM > 1
     if (hy_order > 1 .and. hy_killdivB .and. (.not. hy_forceHydroLimit)) then
        lbx = lbound(scrchFaceXPtr); ubx = ubound(scrchFaceXPtr)
        lby = lbound(scrchFaceYPtr); uby = ubound(scrchFaceYPtr)
        if (NDIM == 3) then
           lbz = lbound(scrchFaceZPtr); ubz = ubound(scrchFaceZPtr)
        end if
        @scrch_aslice(scrchFaceXPtr,HY_P06_FACEXPTR_VAR)@ = @bx_mslice(Bx,MAGI_FACE_VAR,@lbx_ref(lbx,2)@+1,@lbx_ref(ubx,2)@+1,@lbx_ref(lbx,3)@,@lbx_ref(ubx,3)@,@lbx_ref(lbx,4)@,@lbx_ref(ubx,4)@)@
        @scrch_aslice(scrchFaceXPtr,HY_N06_FACEXPTR_VAR)@ = @bx_mslice(Bx,MAGI_FACE_VAR,@lbx_ref(lbx,2)@,@lbx_ref(ubx,2)@,@lbx_ref(lbx,3)@,@lbx_ref(ubx,3)@,@lbx_ref(lbx,4)@,@lbx_ref(ubx,4)@)@
        @scrch_aslice(scrchFaceYPtr,HY_P07_FACEYPTR_VAR)@ = @bx_mslice(By,MAGI_FACE_VAR,@lbx_ref(lby,2)@,@lbx_ref(uby,2)@,@lbx_ref(lby,3)@+1,@lbx_ref(uby,3)@+1,@lbx_ref(lby,4)@,@lbx_ref(uby,4)@)@
        @scrch_aslice(scrchFaceYPtr,HY_N07_FACEYPTR_VAR)@ = @bx_mslice(By,MAGI_FACE_VAR,@lbx_ref(lby,2)@,@lbx_ref(uby,2)@,@lbx_ref(lby,3)@,@lbx_ref(uby,3)@,@lbx_ref(lby,4)@,@lbx_ref(uby,4)@)@
#if NDIM == 3
        @scrch_aslice(scrchFaceZPtr,HY_P08_FACEZPTR_VAR)@ = @bx_mslice(Bz,MAGI_FACE_VAR,@lbx_ref(lbz,2)@,@lbx_ref(ubz,2)@,@lbx_ref(lbz,3)@,@lbx_ref(ubz,3)@,@lbx_ref(lbz,4)@+1,@lbx_ref(ubz,4)@+1)@
        @scrch_aslice(scrchFaceZPtr,HY_N08_FACEZPTR_VAR)@ = @bx_mslice(Bz,MAGI_FACE_VAR,@lbx_ref(lbz,2)@,@lbx_ref(ubz,2)@,@lbx_ref(lbz,3)@,@lbx_ref(ubz,3)@,@lbx_ref(lbz,4)@,@lbx_ref(ubz,4)@)@
#endif
     endif
#endif
#endif
#endif /* end of ifdef FLASH_USM_MHD */
! MHD only-------------------------------------------------------------------------------
  endif ! end of if (.not. normalFieldUpdateOnly) then

#ifndef CFL_VAR
#if NDIM > 1
  if (.not. normalFieldUpdateOnly) then
#ifdef BDRY_VAR
     !$omp critical (hy_crit_update_hy_cfl)
     hy_cfl = minCfl
     !$omp end critical (hy_crit_update_hy_cfl)
#else
     if (hy_useHybridOrder .OR. hy_fallbackLowerCFL) then
        !$omp critical (hy_crit_update_hy_cfl)
        hy_cfl = minCfl
        !$omp end critical (hy_crit_update_hy_cfl)
     end if
#endif
  end if
#endif
#endif


  !!---------------------------------------------------------------!
  !! (5) Transverse correction terms for 3D -----------------------!
  !!---------------------------------------------------------------!
#if NDIM == 3
  if (.not. normalFieldUpdateOnly) then

     ! hy_use3dFullCTU=.true. will provide CFL <= 1; otherwise, CFL<0.5.
     if (hy_use3dFullCTU) then
        ! Set dt,dx,dy,dz factors
        dt2dxdy6=dt*dt/(6.*@del_ref(DIR_X)@*@del_ref(DIR_Y)@)
        dt2dydz6=dt*dt/(6.*@del_ref(DIR_Y)@*@del_ref(DIR_Z)@)
        dt2dzdx6=dt*dt/(6.*@del_ref(DIR_Z)@*@del_ref(DIR_X)@)


        ! Now let's compute cross derivatives for CTU
        do k=k0-2,kmax+2
           do j=j0-2,jmax+2
              do i=i0-2,imax+2

#ifdef CFL_VAR
                 cellCfl = @fld(U,CFL_VAR,+0,+0,+0)@
#else
                 cellCfl = hy_cfl
#endif
                 !! save the cell center values for later use
                 @vc_rng(HY_DENS,HY_END_VARS-kGrav)@ = &
                      (/@fld(U,DENS_VAR,+0,+0,+0)@&
                       ,@fld_rng(U,VELX_VAR,VELZ_VAR,+0,+0,+0)@&
                       ,@fld(U,PRES_VAR,+0,+0,+0)@&
#ifdef FLASH_USM_MHD || FLASH_UGLM_MHD
                       ,@fld_rng(U,MAGX_VAR,MAGZ_VAR,+0,+0,+0)@&
#endif
#ifdef FLASH_UGLM_MHD
                       ,@fld(U,GLMP_VAR,+0,+0,+0)@ &
#endif
                       ,@fld_rng(U,GAMC_VAR,GAME_VAR,+0,+0,+0)@ &
                       ,@fld(U,EINT_VAR,+0,+0,+0)@ &
#ifdef FLASH_UHD_3T
                       ,@fld(U,EELE_VAR,+0,+0,+0)@ &
                       ,@fld(U,EION_VAR,+0,+0,+0)@ &
                       ,@fld(U,ERAD_VAR,+0,+0,+0)@ &
#endif
                       /)

                 !! ============ X-direction ==================================================================
                 If (i .ge. i0-1 .and. i .le. imax+1) then
#ifdef FLASH_UHD_HYDRO
                 If ((j .ge. j0   .and. j .le. jmax  ) .and. (k .ge. k0 .and. k .le. kmax)) then
#endif
                 ! YZ cross dervatives for X states
                 SigmPtr => @sig_ref(sig,DIR_Z,+0,-1,+0)@
                 SigcPtr => @sig_ref(sig,DIR_Z,+0,+0,+0)@
                 SigpPtr => @sig_ref(sig,DIR_Z,+0,+1,+0)@

                 call  hy_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       @lam_ref(lambda,DIR_Y,+0,+0,+0)@,&
                       @leig_ref(leftEig,DIR_Y,+0,+0,+0)@,&
                       @reig_ref(rghtEig,DIR_Y,+0,+0,+0)@,&
                       HY_END_VARS,@tf1(TransFluxYZ,1)@)

                 ! ZY cross derivatives for X states
                 SigmPtr => @sig_ref(sig,DIR_Y,+0,+0,-1)@
                 SigcPtr => @sig_ref(sig,DIR_Y,+0,+0,+0)@
                 SigpPtr => @sig_ref(sig,DIR_Y,+0,+0,+1)@

                 call hy_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       @lam_ref(lambda,DIR_Z,+0,+0,+0)@,&
                       @leig_ref(leftEig,DIR_Z,+0,+0,+0)@,&
                       @reig_ref(rghtEig,DIR_Z,+0,+0,+0)@,&
                       HY_END_VARS,@tf1(TransFluxZY,1)@)

#ifdef FLASH_USM_MHD || FLASH_UGLM_MHD
                 @tf1(TransFluxYZ,HY_MAGX)@ = 0.
                 @tf1(TransFluxZY,HY_MAGX)@ = 0.
#endif
                 @fld_rng(scrchFaceXPtr,HY_P01_FACEXPTR_VAR,HY_P01_FACEXPTR_VAR+HY_SCRATCH_NUM-1,+0,+0,+0)@=&
                      @fld_rng(scrchFaceXPtr,HY_P01_FACEXPTR_VAR,HY_P01_FACEXPTR_VAR+HY_SCRATCH_NUM-1,+0,+0,+0)@&
                      +(@tf_rng(TransFluxYZ,HY_DENS,HY_END_VARS-kGrav)@+@tf_rng(TransFluxZY,HY_DENS,HY_END_VARS-kGrav)@)*dt2dydz6

                 @fld_rng(scrchFaceXPtr,HY_N01_FACEXPTR_VAR,HY_N01_FACEXPTR_VAR+HY_SCRATCH_NUM-1,+0,+0,+0)@=&
                      @fld_rng(scrchFaceXPtr,HY_N01_FACEXPTR_VAR,HY_N01_FACEXPTR_VAR+HY_SCRATCH_NUM-1,+0,+0,+0)@&
                      +(@tf_rng(TransFluxYZ,HY_DENS,HY_END_VARS-kGrav)@+@tf_rng(TransFluxZY,HY_DENS,HY_END_VARS-kGrav)@)*dt2dydz6

                 !! CHECK FOR NEGATIVITY OF DENSITY AND PRESSURE IN X-DIRECTION
                 IF (@fld(scrchFaceXPtr,HY_P01_FACEXPTR_VAR,+0,+0,+0)@ .le. 0. .or. &
                     @fld(scrchFaceXPtr,HY_N01_FACEXPTR_VAR,+0,+0,+0)@ .le. 0. .or. &
                     @fld(scrchFaceXPtr,HY_P05_FACEXPTR_VAR,+0,+0,+0)@ .le. 0. .or. &
                     @fld(scrchFaceXPtr,HY_N05_FACEXPTR_VAR,+0,+0,+0)@ .le. 0. ) THEN

                    call fallbackToFirstOrder(DIR_X,&
                         @scrch_section(scrchFaceXPtr,HY_N01_FACEXPTR_VAR,+0,+0,+0)@,&
                         @scrch_section(scrchFaceXPtr,HY_P01_FACEXPTR_VAR,+0,+0,+0)@,&
                         Vc, &
                         hy_SpcL,hy_SpcR,U,i,j,k)

#if (NSPECIES+NMASS_SCALARS) > 0
                 ELSE if (hy_fullSpecMsFluxHandling) then
                 ! YZ cross dervatives for X states
                    SigmPtr => @spcsig_ref(hy_SpcSig,DIR_Z,+0,-1,+0)@
                    SigcPtr => @spcsig_ref(hy_SpcSig,DIR_Z,+0,+0,+0)@
                    SigpPtr => @spcsig_ref(hy_SpcSig,DIR_Z,+0,+1,+0)@

                    call hy_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       @lam_ref(lambda,DIR_Y,+0,+0,+0)@,&
                       @leig_ref(leftEig,DIR_Y,+0,+0,+0)@,&
                       @reig_ref(rghtEig,DIR_Y,+0,+0,+0)@,&
                       HY_NSPEC,@tf1(TransFluxYZ,HY_SPEC_BEG)@,&
                       speciesScalar=.true.)

                 ! ZY cross derivatives for X states
                    SigmPtr => @spcsig_ref(hy_SpcSig,DIR_Y,+0,+0,-1)@
                    SigcPtr => @spcsig_ref(hy_SpcSig,DIR_Y,+0,+0,+0)@
                    SigpPtr => @spcsig_ref(hy_SpcSig,DIR_Y,+0,+0,+1)@

                    call hy_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       @lam_ref(lambda,DIR_Z,+0,+0,+0)@,&
                       @leig_ref(leftEig,DIR_Z,+0,+0,+0)@,&
                       @reig_ref(rghtEig,DIR_Z,+0,+0,+0)@,&
                       HY_NSPEC,@tf1(TransFluxZY,HY_SPEC_BEG)@,&
                       speciesScalar=.true.)

                    @spc_rng_ref(hy_SpcR,1,HY_NSPEC,DIR_X,+0,+0,+0)@ = @spc_rng_ref(hy_SpcR,1,HY_NSPEC,DIR_X,+0,+0,+0)@&
                      +(@tf_rng(TransFluxYZ,HY_SPEC_BEG,HY_SPEC_END)@+@tf_rng(TransFluxZY,HY_SPEC_BEG,HY_SPEC_END)@)*dt2dydz6

                    @spc_rng_ref(hy_SpcL,1,HY_NSPEC,DIR_X,+0,+0,+0)@ = @spc_rng_ref(hy_SpcL,1,HY_NSPEC,DIR_X,+0,+0,+0)@&
                      +(@tf_rng(TransFluxYZ,HY_SPEC_BEG,HY_SPEC_END)@+@tf_rng(TransFluxZY,HY_SPEC_BEG,HY_SPEC_END)@)*dt2dydz6
#endif
                 ENDIF
#ifdef FLASH_UHD_HYDRO
                 Endif
#endif
                 Endif

                 !! ============ Y-direction ==================================================================
                                  If (j .ge. j0-1 .and. j .le. jmax+1) then
                 #ifdef FLASH_UHD_HYDRO
                                  If ((i .ge. i0   .and. i .le. imax  ) .and. (k .ge. k0 .and. k .le. kmax)) then
                 #endif
                                  ! ZX cross derivatives for Y states
                                  SigmPtr => @sig_ref(sig,DIR_X,+0,+0,-1)@
                                  SigcPtr => @sig_ref(sig,DIR_X,+0,+0,+0)@
                                  SigpPtr => @sig_ref(sig,DIR_X,+0,+0,+1)@

                                  call  hy_upwindTransverseFlux&
                                       (dir,transOrder3D,&
                                        SigmPtr,SigcPtr,SigpPtr,&
                                        @lam_ref(lambda,DIR_Z,+0,+0,+0)@,&
                                        @leig_ref(leftEig,DIR_Z,+0,+0,+0)@,&
                                        @reig_ref(rghtEig,DIR_Z,+0,+0,+0)@,&
                                        HY_END_VARS,@tf1(TransFluxZX,1)@)

                                  ! XZ cross derivatives for Y states
                                  SigmPtr => @sig_ref(sig,DIR_Z,-1,+0,+0)@
                                  SigcPtr => @sig_ref(sig,DIR_Z,+0,+0,+0)@
                                  SigpPtr => @sig_ref(sig,DIR_Z,+1,+0,+0)@

                                  call  hy_upwindTransverseFlux&
                                       (dir,transOrder3D,&
                                        SigmPtr,SigcPtr,SigpPtr,&
                                        @lam_ref(lambda,DIR_X,+0,+0,+0)@,&
                                        @leig_ref(leftEig,DIR_X,+0,+0,+0)@,&
                                        @reig_ref(rghtEig,DIR_X,+0,+0,+0)@,&
                                        HY_END_VARS,@tf1(TransFluxXZ,1)@)

#ifdef FLASH_USM_MHD || FLASH_UGLM_MHD
                                  @tf1(TransFluxZX,HY_MAGY)@ = 0.
                                  @tf1(TransFluxXZ,HY_MAGY)@ = 0.
#endif
                                  @fld_rng(scrchFaceYPtr,HY_P01_FACEYPTR_VAR,HY_P01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,+0,+0,+0)@=&
                                       @fld_rng(scrchFaceYPtr,HY_P01_FACEYPTR_VAR,HY_P01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,+0,+0,+0)@&
                                       + (@tf_rng(TransFluxZX,HY_DENS,HY_END_VARS-kGrav)@+@tf_rng(TransFluxXZ,HY_DENS,HY_END_VARS-kGrav)@)*dt2dzdx6

                                  @fld_rng(scrchFaceYPtr,HY_N01_FACEYPTR_VAR,HY_N01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,+0,+0,+0)@=&
                                       @fld_rng(scrchFaceYPtr,HY_N01_FACEYPTR_VAR,HY_N01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,+0,+0,+0)@&
                                       + (@tf_rng(TransFluxZX,HY_DENS,HY_END_VARS-kGrav)@+@tf_rng(TransFluxXZ,HY_DENS,HY_END_VARS-kGrav)@)*dt2dzdx6

                                  !! CHECK FOR NEGATIVITY OF DENSITY AND PRESSURE IN Y-DIRECTION
                                  IF (@fld(scrchFaceYPtr,HY_P01_FACEYPTR_VAR,+0,+0,+0)@ .le. 0. .or. &
                                      @fld(scrchFaceYPtr,HY_N01_FACEYPTR_VAR,+0,+0,+0)@ .le. 0. .or. &
                                      @fld(scrchFaceYPtr,HY_P05_FACEYPTR_VAR,+0,+0,+0)@ .le. 0. .or. &
                                      @fld(scrchFaceYPtr,HY_N05_FACEYPTR_VAR,+0,+0,+0)@ .le. 0. ) THEN

                                     call fallbackToFirstOrder(DIR_Y,&
                                          @scrch_section(scrchFaceYPtr,HY_N01_FACEYPTR_VAR,+0,+0,+0)@,&
                                          @scrch_section(scrchFaceYPtr,HY_P01_FACEYPTR_VAR,+0,+0,+0)@,&
                                          Vc, &
                                          hy_SpcL,hy_SpcR,U,i,j,k)

#if (NSPECIES+NMASS_SCALARS) > 0
                                  ELSE if (hy_fullSpecMsFluxHandling) then
                                  !ZX cross dervatives for X states
                                     SigmPtr => @spcsig_ref(hy_SpcSig,DIR_X,+0,+0,-1)@
                                     SigcPtr => @spcsig_ref(hy_SpcSig,DIR_X,+0,+0,+0)@
                                     SigpPtr => @spcsig_ref(hy_SpcSig,DIR_X,+0,+0,+1)@

                                     call hy_upwindTransverseFlux&
                                       (dir,transOrder3D,&
                                        SigmPtr,SigcPtr,SigpPtr,&
                                        @lam_ref(lambda,DIR_Z,+0,+0,+0)@,&
                                        @leig_ref(leftEig,DIR_Z,+0,+0,+0)@,&
                                        @reig_ref(rghtEig,DIR_Z,+0,+0,+0)@,&
                                     @spc_rng_ref(hy_SpcR,1,HY_NSPEC,DIR_Y,+0,+0,+0)@ = @spc_rng_ref(hy_SpcR,1,HY_NSPEC,DIR_Y,+0,+0,+0)@&
                                       +(@tf_rng(TransFluxZX,HY_SPEC_BEG,HY_SPEC_END)@+@tf_rng(TransFluxXZ,HY_SPEC_BEG,HY_SPEC_END)@)*dt2dydz6

                                     @spc_rng_ref(hy_SpcL,1,HY_NSPEC,DIR_Y,+0,+0,+0)@ = @spc_rng_ref(hy_SpcL,1,HY_NSPEC,DIR_Y,+0,+0,+0)@&
                                       +(@tf_rng(TransFluxZX,HY_SPEC_BEG,HY_SPEC_END)@+@tf_rng(TransFluxXZ,HY_SPEC_BEG,HY_SPEC_END)@)*dt2dydz6
#endif
                                  ENDIF
#ifdef FLASH_UHD_HYDRO
                                  Endif
#endif
                                  Endif


                 !! ============ z-direction ==================================================================
                 If (k .ge. k0-1 .and. k .le. kmax+1) then
#ifdef FLASH_UHD_HYDRO
                 If ((i .ge. i0   .and. i .le. imax  ) .and. (j .ge. j0 .and. j .le. jmax)) then
#endif
                 ! XY cross derivatives for Z states
                 SigmPtr => @sig_ref(sig,DIR_Y,-1,+0,+0)@
                 SigcPtr => @sig_ref(sig,DIR_Y,+0,+0,+0)@
                 SigpPtr => @sig_ref(sig,DIR_Y,+1,+0,+0)@

                 call  hy_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       @lam_ref(lambda,DIR_X,+0,+0,+0)@,&
                       @leig_ref(leftEig,DIR_X,+0,+0,+0)@,&
                       @reig_ref(rghtEig,DIR_X,+0,+0,+0)@,&
                       HY_END_VARS,@tf1(TransFluxXY,1)@)

                 ! YX cross derivatives for Z states
                 SigmPtr => @sig_ref(sig,DIR_X,+0,-1,+0)@
                 SigcPtr => @sig_ref(sig,DIR_X,+0,+0,+0)@
                 SigpPtr => @sig_ref(sig,DIR_X,+0,+1,+0)@

                 call  hy_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       @lam_ref(lambda,DIR_Y,+0,+0,+0)@,&
                       @leig_ref(leftEig,DIR_Y,+0,+0,+0)@,&
                       @reig_ref(rghtEig,DIR_Y,+0,+0,+0)@,&
                       HY_END_VARS,@tf1(TransFluxYX,1)@)

#ifdef FLASH_USM_MHD|| FLASH_UGLM_MHD
                 @tf1(TransFluxXY,HY_MAGZ)@ = 0.
                 @tf1(TransFluxYX,HY_MAGZ)@ = 0.
#endif

                 @fld_rng(scrchFaceZPtr,HY_P01_FACEZPTR_VAR,HY_P01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,+0,+0,+0)@=&
                      @fld_rng(scrchFaceZPtr,HY_P01_FACEZPTR_VAR,HY_P01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,+0,+0,+0)@&
                      + (@tf_rng(TransFluxXY,HY_DENS,HY_END_VARS-kGrav)@+@tf_rng(TransFluxYX,HY_DENS,HY_END_VARS-kGrav)@)*dt2dxdy6

                 @fld_rng(scrchFaceZPtr,HY_N01_FACEZPTR_VAR,HY_N01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,+0,+0,+0)@=&
                      @fld_rng(scrchFaceZPtr,HY_N01_FACEZPTR_VAR,HY_N01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,+0,+0,+0)@&
                      + (@tf_rng(TransFluxXY,HY_DENS,HY_END_VARS-kGrav)@+@tf_rng(TransFluxYX,HY_DENS,HY_END_VARS-kGrav)@)*dt2dxdy6


                 !! CHECK FOR NEGATIVITY OF DENSITY AND PRESSURE IN Z-DIRECTION
                 IF (@fld(scrchFaceZPtr,HY_P01_FACEZPTR_VAR,+0,+0,+0)@ .le. 0. .or. &
                     @fld(scrchFaceZPtr,HY_N01_FACEZPTR_VAR,+0,+0,+0)@ .le. 0. .or. &
                     @fld(scrchFaceZPtr,HY_P05_FACEZPTR_VAR,+0,+0,+0)@ .le. 0. .or. &
                     @fld(scrchFaceZPtr,HY_N05_FACEZPTR_VAR,+0,+0,+0)@ .le. 0. ) THEN

                    call fallbackToFirstOrder(DIR_Z,&
                         @scrch_section(scrchFaceZPtr,HY_N01_FACEZPTR_VAR,+0,+0,+0)@,&
                         @scrch_section(scrchFaceZPtr,HY_P01_FACEZPTR_VAR,+0,+0,+0)@,&
                         Vc, &
                         hy_SpcL,hy_SpcR,U,i,j,k)


#if (NSPECIES+NMASS_SCALARS) > 0
                 ELSE if (hy_fullSpecMsFluxHandling) then
                 ! XY cross dervatives for X states
                    SigmPtr => @spcsig_ref(hy_SpcSig,DIR_Y,-1,+0,+0)@
                    SigcPtr => @spcsig_ref(hy_SpcSig,DIR_Y,+0,+0,+0)@
                    SigpPtr => @spcsig_ref(hy_SpcSig,DIR_Y,+1,+0,+0)@

                    call hy_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       @lam_ref(lambda,DIR_X,+0,+0,+0)@,&
                       @leig_ref(leftEig,DIR_X,+0,+0,+0)@,&
                       @reig_ref(rghtEig,DIR_X,+0,+0,+0)@,&
                       HY_NSPEC,@tf1(TransFluxXY,HY_SPEC_BEG)@,&
                       speciesScalar=.true.)

                 ! YX cross derivatives for Y states
                    SigmPtr => @spcsig_ref(hy_SpcSig,DIR_X,+0,-1,+0)@
                    SigcPtr => @spcsig_ref(hy_SpcSig,DIR_X,+0,+0,+0)@
                    SigpPtr => @spcsig_ref(hy_SpcSig,DIR_X,+0,+1,+0)@

                    call hy_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       @lam_ref(lambda,DIR_Y,+0,+0,+0)@,&
                       @leig_ref(leftEig,DIR_Y,+0,+0,+0)@,&
                       @reig_ref(rghtEig,DIR_Y,+0,+0,+0)@,&
                       HY_NSPEC,@tf1(TransFluxYX,HY_SPEC_BEG)@,&
                       speciesScalar=.true.)

                    @spc_rng_ref(hy_SpcR,1,HY_NSPEC,DIR_Z,+0,+0,+0)@ = @spc_rng_ref(hy_SpcR,1,HY_NSPEC,DIR_Z,+0,+0,+0)@&
                      +(@tf_rng(TransFluxXY,HY_SPEC_BEG,HY_SPEC_END)@+@tf_rng(TransFluxYX,HY_SPEC_BEG,HY_SPEC_END)@)*dt2dydz6

                    @spc_rng_ref(hy_SpcL,1,HY_NSPEC,DIR_Z,+0,+0,+0)@ = @spc_rng_ref(hy_SpcL,1,HY_NSPEC,DIR_Z,+0,+0,+0)@&
                      +(@tf_rng(TransFluxXY,HY_SPEC_BEG,HY_SPEC_END)@+@tf_rng(TransFluxYX,HY_SPEC_BEG,HY_SPEC_END)@)*dt2dydz6
#endif
                 ENDIF
#ifdef FLASH_UHD_HYDRO
                 Endif
#endif
                 Endif

                 if (hy_fallbackLowerCFL) then
#ifdef CFL_VAR
                    @fld(U,CFL_VAR,+0,+0,+0)@ = cellCfl
#else
                    minCfl = min(minCfl,cellCfl)
#endif
                 end if

              enddo ! i-loop
           enddo ! j-loop
        enddo ! k-loop

#ifndef CFL_VAR
        if (hy_fallbackLowerCFL) then
           !$omp critical (hy_crit_update_hy_cfl)
           hy_cfl = min(hy_cfl,minCfl)
           !$omp end critical (hy_crit_update_hy_cfl)
        end if
#endif

     end if !end of if (hy_use3dFullCTU) then

  endif ! End of if (.not. normalFieldUpdateOnly) then

#endif /* end of #if NDIM == 3 */




!#ifdef GLMGLM
#ifdef FLASH_UGLM_MHD
  do k=k0-2-k3+kHydro*k3,kmax+2+k3-kHydro*k3
     do j=j0-2-k3+kHydro*k2,jmax+2+k3-kHydro*k2
        do i=i0-2-k3+kHydro,imax+2+k3-kHydro
           ! Extra stencil is needed for 3D to correctly calculate transverse fluxes
           !(i.e., cross derivatives in x,y, & z)

           ! (1) costruct Godunov flux for GLM-Psi (GMLP_VAR)
           @spc3(BxStar,+0,+0,+0)@ = ( @fld(scrchFaceXPtr,HY_P06_FACEXPTR_VAR,-1,+0,+0)@&
                            +@fld(scrchFaceXPtr,HY_N06_FACEXPTR_VAR,+0,+0,+0)@)*0.5

           @spc3(BxStar,+0,+0,+0)@ = @spc3(BxStar,+0,+0,+0)@ - 0.5/hy_C_hyp*&
                           ( @fld(scrchFaceXPtr,HY_N09_FACEXPTR_VAR,+0,+0,+0)@ &
                            -@fld(scrchFaceXPtr,HY_P09_FACEXPTR_VAR,-1,+0,+0)@)
#if NDIM > 1
           @spc3(ByStar,+0,+0,+0)@ = ( @fld(scrchFaceYPtr,HY_P07_FACEYPTR_VAR,+0,-1,+0)@&
                            +@fld(scrchFaceYPtr,HY_N07_FACEYPTR_VAR,+0,+0,+0)@)*0.5

           @spc3(ByStar,+0,+0,+0)@ = @spc3(ByStar,+0,+0,+0)@ - 0.5/hy_C_hyp*&
                           ( @fld(scrchFaceYPtr,HY_N09_FACEYPTR_VAR,+0,+0,+0)@ &
                            -@fld(scrchFaceYPtr,HY_P09_FACEYPTR_VAR,+0,-1,+0)@)
#if NDIM == 3
           @spc3(BzStar,+0,+0,+0)@ = ( @fld(scrchFaceZPtr,HY_P08_FACEZPTR_VAR,+0,+0,-1)@&
                            +@fld(scrchFaceZPtr,HY_N08_FACEZPTR_VAR,+0,+0,+0)@)*0.5

           @spc3(BzStar,+0,+0,+0)@ = @spc3(BzStar,+0,+0,+0)@ - 0.5/hy_C_hyp*&
                           ( @fld(scrchFaceZPtr,HY_N09_FACEZPTR_VAR,+0,+0,+0)@ &
                            -@fld(scrchFaceZPtr,HY_P09_FACEZPTR_VAR,+0,+0,-1)@)
#endif
#endif


           ! (2) costruct Godunov flux for normal magnetic fields
           @spc3(GLMxStar,+0,+0,+0)@ = ( @fld(scrchFaceXPtr,HY_P09_FACEXPTR_VAR,-1,+0,+0)@&
                              +@fld(scrchFaceXPtr,HY_N09_FACEXPTR_VAR,+0,+0,+0)@)*0.5

           @spc3(GLMxStar,+0,+0,+0)@ = @spc3(GLMxStar,+0,+0,+0)@ - 0.5*hy_C_hyp*&
                             ( @fld(scrchFaceXPtr,HY_N06_FACEXPTR_VAR,+0,+0,+0)@ &
                              -@fld(scrchFaceXPtr,HY_P06_FACEXPTR_VAR,-1,+0,+0)@)
#if NDIM > 1
           @spc3(GLMyStar,+0,+0,+0)@ = ( @fld(scrchFaceYPtr,HY_P09_FACEYPTR_VAR,+0,-1,+0)@&
                              +@fld(scrchFaceYPtr,HY_N09_FACEYPTR_VAR,+0,+0,+0)@)*0.5

           @spc3(GLMyStar,+0,+0,+0)@ = @spc3(GLMyStar,+0,+0,+0)@ - 0.5*hy_C_hyp*&
                             ( @fld(scrchFaceYPtr,HY_N07_FACEYPTR_VAR,+0,+0,+0)@ &
                              -@fld(scrchFaceYPtr,HY_P07_FACEYPTR_VAR,+0,-1,+0)@)
#if NDIM == 3
           @spc3(GLMzStar,+0,+0,+0)@ = ( @fld(scrchFaceZPtr,HY_P09_FACEZPTR_VAR,+0,+0,-1)@&
                              +@fld(scrchFaceZPtr,HY_N09_FACEZPTR_VAR,+0,+0,+0)@)*0.5

           @spc3(GLMzStar,+0,+0,+0)@ = @spc3(GLMzStar,+0,+0,+0)@ - 0.5*hy_C_hyp*&
                             ( @fld(scrchFaceZPtr,HY_N08_FACEZPTR_VAR,+0,+0,+0)@ &
                              -@fld(scrchFaceZPtr,HY_P08_FACEZPTR_VAR,+0,+0,-1)@)
#endif
#endif

        enddo
     enddo
  enddo


  do k=k0-2-k3+kHydro*k3,kmax+2+k3-kHydro*k3
     do j=j0-2-k3+kHydro*k2,jmax+2+k3-kHydro*k2
        do i=i0-2-k3+kHydro,imax+2+k3-kHydro
           ! GLM psi x-Riemann state
           @fld(scrchFaceXPtr,HY_P09_FACEXPTR_VAR,+0,+0,+0)@ = @spc3(GLMxSTar,+1,+0,+0)@
           @fld(scrchFaceXPtr,HY_N09_FACEXPTR_VAR,+0,+0,+0)@ = @spc3(GLMxStar,+0,+0,+0)@

           ! magx
           @fld(scrchFaceXPtr,HY_P06_FACEXPTR_VAR,+0,+0,+0)@ = @spc3(BxStar,+1,+0,+0)@
           @fld(scrchFaceXPtr,HY_N06_FACEXPTR_VAR,+0,+0,+0)@ = @spc3(BxStar,+0,+0,+0)@

#if NDIM > 1
           ! GLM psi y-Riemann state
           @fld(scrchFaceYPtr,HY_P09_FACEYPTR_VAR,+0,+0,+0)@ = @spc3(GLMySTar,+0,+1,+0)@
           @fld(scrchFaceYPtr,HY_N09_FACEYPTR_VAR,+0,+0,+0)@ = @spc3(GLMyStar,+0,+0,+0)@

           @fld(scrchFaceYPtr,HY_P07_FACEYPTR_VAR,+0,+0,+0)@ = @spc3(ByStar,+0,+1,+0)@
           @fld(scrchFaceYPtr,HY_N07_FACEYPTR_VAR,+0,+0,+0)@ = @spc3(ByStar,+0,+0,+0)@

#if NDIM == 3
           ! GLM psi z-Riemann state
           @fld(scrchFaceZPtr,HY_P09_FACEZPTR_VAR,+0,+0,+0)@ = @spc3(GLMzSTar,+0,+0,+1)@
           @fld(scrchFaceZPtr,HY_N09_FACEZPTR_VAR,+0,+0,+0)@ = @spc3(GLMzStar,+0,+0,+0)@

           @fld(scrchFaceZPtr,HY_P08_FACEZPTR_VAR,+0,+0,+0)@ = @spc3(BzStar,+0,+0,+1)@
           @fld(scrchFaceZPtr,HY_N08_FACEZPTR_VAR,+0,+0,+0)@ = @spc3(BzStar,+0,+0,+0)@
#endif
#endif

        enddo
     enddo
  enddo

#endif /* ifdef FLASH_UGLM_MHD */
!#endif


  !! Release pointers
!!$  call Grid_releaseBlkPtr(blockID,U,CENTER)

!!$  ! MHD only-------------------------------------------------------------------------------
!!$#if defined(FLASH_USM_MHD) && NFACE_VARS > 0 && NDIM > 1
!!$  if (hy_order > 1) then
!!$     call Grid_releaseBlkPtr(blockID,Bx,FACEX)
!!$     call Grid_releaseBlkPtr(blockID,By,FACEY)
!!$     if (NDIM == 3) call Grid_releaseBlkPtr(blockID,Bz,FACEZ)
!!$  endif ! if (hy_order > 1) then
!!$#endif /* endif of if defined(FLASH_USM_MHD) && NFACE_VARS > 0 && NDIM > 1 */
!!$  ! MHD only-------------------------------------------------------------------------------


  !! Deallocate arrays
  !deallocate(DivU)
  if (.NOT. normalFieldUpdateOnly) then
     deallocate(sig)
     deallocate(lambda)
     deallocate(leftEig)
     deallocate(rghtEig)
  end if

contains
#include "FortranLangFeatures.fh"
  subroutine fallbackToFirstOrder(iDir,Wleft,Wright,Vc,spcL,spcR,U,i,j,k)
    integer,intent(IN)                        :: iDir
    real,   intent(OUT), dimension(:)         :: Wleft,Wright
    real,   intent(IN),  dimension(:)         :: Vc
    real,   POINTER_INTENT_IN, dimension(:,:,:,:,:),OPTIONAL :: spcL,spcR
    real,   POINTER_INTENT_IN, dimension(:,:,:,:)  ,OPTIONAL :: U
    integer,intent(IN),                       OPTIONAL :: i,j,k

    @wright_rng(Wleft,HY_DENS,HY_END_VARS-kGrav)@ = @vc_rng(HY_DENS,HY_END_VARS-kGrav)@
    @wright_rng(Wright,HY_DENS,HY_END_VARS-kGrav)@ = @vc_rng(HY_DENS,HY_END_VARS-kGrav)@
    if (hy_fullSpecMsFluxHandling .AND. hy_numXN > 0 &
         .AND. present(spcR)) then
       @spc_section(spcL,iDir,+0,+0,+0)@ = &
            @fld_rng(U,SPECIES_BEGIN,MASS_SCALARS_END,+0,+0,+0)@
       @spc_section(spcR,iDir,+0,+0,+0)@ = &
            @fld_rng(U,SPECIES_BEGIN,MASS_SCALARS_END,+0,+0,+0)@
    end if
    cellCfl = min(cellCfl, hy_cflFallbackFactor / real(NDIM))

  end subroutine fallbackToFirstOrder
End Subroutine hy_getRiemannState
