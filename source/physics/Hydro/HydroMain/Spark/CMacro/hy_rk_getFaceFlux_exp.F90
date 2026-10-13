!! Compute directional face fluxes for one Runge-Kutta stage using pencil work arrays.
!!Reorder(4):hy_starState,hy_fl[xyz]
subroutine hy_rk_getFaceFlux(stage, hy_starState, flat3d, hy_flx, hy_fly, hy_flz, &
                             limits, del, snake, flux, uPlusArray, uMinusArray, &
                             loGC)

   use Hydro_data, ONLY : hy_hybridRiemann, hy_cvisc, hy_C_hyp, hy_tiny, &
                          hy_smalldens, hy_smallpres, hy_smallX
   use Timers_interface, ONLY : Timers_start, Timers_stop
 
   implicit none
 
#include "Simulation.h"
#include "constants.h"
#include "Spark.h"
#define NRECON HY_NUM_VARS+NSPECIES+NMASS_SCALARS
 
   integer, intent(IN) :: stage, loGC(3)
   integer, intent(IN), dimension(LOW:HIGH,MDIM,MAXSTAGE) :: limits
   real, dimension(MDIM), intent(IN) :: del
   real, dimension(1:,loGC(1):,loGC(2):,loGC(3):), intent(IN) :: hy_starState
   real, dimension(1:,loGC(1):,loGC(2):,loGC(3):), intent(OUT) :: hy_flx, hy_fly, hy_flz
   real, dimension(:,:,:,:), intent(OUT), target :: uPlusArray, uMinusArray
   real, dimension(:,:,:,:), intent(OUT) :: snake
   real, dimension(:,:,:,:), intent(OUT), target :: flux
   real, dimension(loGC(1):,loGC(2):,loGC(3):), intent(IN) :: flat3d
   real, dimension(size(snake,2),size(snake,3),size(snake,4)) :: flat, grv, shck
   integer, dimension(LOW:HIGH,MDIM) :: activeLimits
 
   integer :: i1,i2,i3, n, g, dir, ierr
   integer, dimension(3) :: guardCells
   integer, dimension(LOW:HIGH,MDIM) :: dirLims
   integer :: i_global, j_global, k_global
   character(len = 2) :: dir_str

   interface
      subroutine reconstruct(i1,i2,i3,uPlusArray,uMinusArray,snake,flat)
         integer, intent(IN) :: i1,i2,i3
         real, dimension(:,:,:,:), intent(INOUT) :: uPlusArray, uMinusArray, snake
         real, dimension(:,:,:), intent(INOUT) :: flat
      end subroutine reconstruct

      subroutine riemann(i1,i2,i3,dir,uPlusArray,uMinusArray,shck,flux, &
                         hy_tiny,hy_hybridRiemann,hy_C_hyp)
         integer, intent(IN) :: i1,i2,i3,dir
         real, dimension(:,:,:,:), intent(INOUT) :: uPlusArray, uMinusArray
         real, dimension(:,:,:), intent(INOUT) :: shck
         real, dimension(:,:,:,:), intent(INOUT), target :: flux
         real, intent(IN) :: hy_tiny, hy_C_hyp
         logical, intent(IN) :: hy_hybridRiemann
      end subroutine riemann

      subroutine initializeSnake(dir,dirLims,guardCells,snake,grv,shck,flat, &
                                 hy_starState,flat3d,loGC)
         integer, intent(IN) :: dir
         integer, dimension(3), intent(IN) :: guardCells, loGC
         integer, dimension(LOW:HIGH,MDIM), intent(IN) :: dirLims
         real, dimension(:,:,:,:), intent(INOUT) :: snake
         real, dimension(:,:,:), intent(INOUT) :: grv, shck, flat
         real, dimension(1:,loGC(1):,loGC(2):,loGC(3):), intent(IN) :: hy_starState
         real, dimension(loGC(1):,loGC(2):,loGC(3):), intent(IN) :: flat3d
      end subroutine initializeSnake
   end interface

   activeLimits(:,:) = limits(:,:,stage)
   
   !$omp target data map(to: dir, dirLims, guardCells)



   ! !  Begin loop over zones
   do dir = 1,  NDIM

      select case(dir)
      case (IAXIS)
         dirLims(:,:) = activeLimits(:,:)
         guardCells(1) = activeLimits(LOW,IAXIS) - loGC(IAXIS)
         guardCells(2) = activeLimits(LOW,JAXIS) - loGC(JAXIS)
         guardCells(3) = activeLimits(LOW,KAXIS) - loGC(KAXIS)
         dir_str = "_x"
      case (JAXIS)
         dirLims(:,1) = activeLimits(:,JAXIS)
         dirLims(:,2) = activeLimits(:,IAXIS)
         dirLims(:,3) = activeLimits(:,KAXIS)
         guardCells(1) = activeLimits(LOW,JAXIS) - loGC(JAXIS)
         guardCells(2) = activeLimits(LOW,IAXIS) - loGC(IAXIS)
         guardCells(3) = activeLimits(LOW,KAXIS) - loGC(KAXIS)
         dir_str = "_y"
      case (KAXIS)
         dirLims(:,1) = activeLimits(:,KAXIS)
         dirLims(:,2) = activeLimits(:,IAXIS)
         dirLims(:,3) = activeLimits(:,JAXIS)
         guardCells(1) = activeLimits(LOW,KAXIS) - loGC(KAXIS)
         guardCells(2) = activeLimits(LOW,IAXIS) - loGC(IAXIS)
         guardCells(3) = activeLimits(LOW,JAXIS) - loGC(JAXIS)
         dir_str = "_z"
      end select 
      !$omp target update to(dir,guardCells,dirLims)
      ! call Timers_start("snake_init"//dir_str) 

      call initializeSnake(dir,dirLims,guardCells,snake,grv,shck,flat, &
                           hy_starState,flat3d,loGC) ! This is offloaded internally

      ! call Timers_stop("snake_init"//dir_str)

      ! call Timers_start("setHSE"//dir_str)
      !$omp target teams distribute parallel do collapse(2) & ! This collapse 2 is because there is a data dependency
      !$omp private(i1,i2,i3) shared(dirLims,dir,guardCells)
      do i3 = 1 + guardCells(3) , 1 + guardCells(3) + dirLims(HIGH,3) - dirLims(LOW,3)
         do i2 = 1 + guardCells(2) , 1 + guardCells(2) + dirLims(HIGH,2) - dirLims(LOW,2)
            do i1 = 1 + guardCells(1) - 1 , 1 + guardCells(1) + dirLims(HIGH,1) - dirLims(LOW,1) + 1
               call setHSE(i1,i2,i3,dir)
            enddo
         enddo
      enddo
      ! call Timers_stop("setHSE"//dir_str)

      ! call Timers_start("recon"//dir_str)

     !$omp target teams distribute parallel do collapse(3) &
     !$omp  private(i1,i2,i3) shared(dirLims,guardCells,uPlusArray, uMinusArray, snake, flat) default(none)
     do i3 =1 + guardCells(3) , 1 + guardCells(3) + dirLims(HIGH,3) - dirLims(LOW,3)
        do i2 =1 + guardCells(2) , 1 + guardCells(2) + dirLims(HIGH,2) - dirLims(LOW,2)
           do i1 =1 + guardCells(1) - 1 , 1 + guardCells(1) + dirLims(HIGH,1) - dirLims(LOW,1) + 1
              call reconstruct(i1,i2,i3,uPlusArray,uMinusArray,snake,flat)
           enddo
        enddo
     enddo
   !   call reconstructOffloaded_split_inlined(dirLims,guardCells)
      !$omp target teams distribute parallel do collapse(3) &
      !$omp private(i1,i2,i3) shared(dirLims,dir,guardCells)
      do i3 = 1 + guardCells(3) , 1 + guardCells(3) + dirLims(HIGH,3) - dirLims(LOW,3)
         do i2 = 1 + guardCells(2) , 1 + guardCells(2) + dirLims(HIGH,2) - dirLims(LOW,2)
            do i1 = 1 + guardCells(1) - 1 , 1 + guardCells(1) + dirLims(HIGH,1) - dirLims(LOW,1) + 1
               call ensurePhysicalState(uPlusArray(:,i1,i2,i3))
               call ensurePhysicalState(uMinusArray(:,i1,i2,i3))
            enddo 
         enddo
      enddo

      ! call Timers_stop("recon"//dir_str)

      ! call Timers_start("resetHSE"//dir_str)
      !$omp target teams distribute parallel do collapse(2) & ! This collapse 2 is because there is a data dependency
      !$omp private(i1,i2,i3) shared(dirLims,dir,guardCells)
      do i3 = 1 + guardCells(3) , 1 + guardCells(3) + dirLims(HIGH,3) - dirLims(LOW,3)
         do i2 = 1 + guardCells(2) , 1 + guardCells(2) + dirLims(HIGH,2) - dirLims(LOW,2)
            do i1 = 1 + guardCells(1) - 1 , 1 + guardCells(1) + dirLims(HIGH,1) - dirLims(LOW,1) + 1
               call resetHSE(i1,i2,i3,dir)
            enddo
         enddo
      enddo
      ! call Timers_stop("resetHSE"//dir_str)

      ! call Timers_start("riemann"//dir_str)
      !  Now call the Riemann solver to compute fluxes
      !$omp target teams distribute parallel do collapse(3) &
      !$omp private(i1,i2,i3) shared(dirLims,dir,guardCells,flux) default(none)
       do i3 = 1 + guardCells(3) ,1 + guardCells(3) + dirLims(HIGH,3) - dirLims(LOW,3)
          do i2 = 1 + guardCells(2) , 1 + guardCells(2) + dirLims(HIGH,2) - dirLims(LOW,2)
             do i1 = 1 + guardCells(1) , 1 + guardCells(1) + dirLims(HIGH,1) - dirLims(LOW,1) + 1
                 call riemann(i1,i2,i3,dir,uPlusArray,uMinusArray,shck,flux, &
                              hy_tiny,hy_hybridRiemann,hy_C_hyp)
            enddo
         enddo
      enddo
      !$omp target teams distribute parallel do collapse(3) &
      !$omp private(i1,i2,i3) shared(dirLims,dir,guardCells,flux) default(none)
      do i3 = 1 + guardCells(3) ,1 + guardCells(3) + dirLims(HIGH,3) - dirLims(LOW,3)
         do i2 = 1 + guardCells(2) , 1 + guardCells(2) + dirLims(HIGH,2) - dirLims(LOW,2)
            do i1 = 1 + guardCells(1) , 1 + guardCells(1) + dirLims(HIGH,1) - dirLims(LOW,1) + 1
               call avisc(i1,i2,i3,dir,NRECON)
               call mscalarFluxes(i1,i2,i3)
           enddo
        enddo
     enddo

      ! call Timers_stop("riemann"//dir_str)

       ! ***************
       ! Fluxes computed for one face of this zone
       ! Save the fluxes
       ! ***************
      !  call Timers_start("saveFluxes"//dir_str)

      !$omp target teams distribute parallel do collapse(3) &
      !$omp private(i1,i2,i3) shared(dirLims,dir,guardCells) default(none)
       do i3 = 1 + guardCells(3) ,1 + guardCells(3) + dirLims(HIGH,3) - dirLims(LOW,3)
         do i2 = 1 + guardCells(2) , 1 + guardCells(2) + dirLims(HIGH,2) - dirLims(LOW,2)
            do i1 = 1 + guardCells(1) , 1 + guardCells(1) + dirLims(HIGH,1) - dirLims(LOW,1) + 1
               call saveFluxes(dir,i1,i2,i3,guardCells,dirLims)
             enddo
          enddo
       enddo 
      !  call Timers_stop("saveFluxes"//dir_str)
      
    enddo
   !$omp end target data 
 contains
 
   !!Set loop dimensions based on the direction of the pencil set
   subroutine setLoop(dir)
     implicit none
     integer, intent(IN) :: dir
     !$omp declare target
     select case(dir)
     case (IAXIS)
        dirLims(:,:) = activeLimits(:,:)
     case (JAXIS)
        dirLims(:,1) = activeLimits(:,JAXIS)
        dirLims(:,2) = activeLimits(:,IAXIS)
        dirLims(:,3) = activeLimits(:,KAXIS)
     case (KAXIS)
        dirLims(:,1) = activeLimits(:,KAXIS)
        dirLims(:,2) = activeLimits(:,IAXIS)
        dirLims(:,3) = activeLimits(:,JAXIS)
     end select
 
   end subroutine setLoop
 
   subroutine saveFluxes(dir,i1,i2,i3,guardCells,dirLims)
     integer, intent(IN) :: dir,i1,i2,i3
     integer, dimension(3) :: guardCells
     integer, dimension(LOW:HIGH,MDIM) :: dirLims
     integer :: i_global, j_global, k_global
     !$omp declare target
     i_global = -1 + i1 + dirLims(LOW,1) - guardCells(1)
     j_global = -1 + i2 + dirLims(LOW,2) - guardCells(2)
     k_global = -1 + i3 + dirLims(LOW,3) - guardCells(3)
     select case(dir)
     case(IAXIS)
        hy_flx(:,i_global,j_global,k_global) = flux(:,i1,i2,i3)
     case (JAXIS)
        hy_fly(:,j_global,i_global,k_global) = flux(:,i1,i2,i3)
     case (KAXIS)
        hy_flz(:,j_global,k_global,i_global) = flux(:,i1,i2,i3)
     end select
   end subroutine saveFluxes
 
!!Account for fluxes that are proportional to mass flux 
subroutine mscalarFluxes(i1,i2,i3)
   implicit none
   real, pointer :: Xstar(:)
   integer, intent(IN) :: i1,i2,i3
   !$omp declare target
#if NSPECIES+NMASS_SCALARS==0
   return
#else
   if (flux(HY_MASS,i1,i2,i3) > 0.) then
      Xstar => uPlusArray(HY_NUM_VARS+1:NRECON,i1-1,i2+0,i3+0)
   else
      Xstar => uMinusArray(HY_NUM_VARS+1:NRECON,i1,i2,i3)
   endif
   flux(HY_NUM_FLUX+1:NFLUXES,i1,i2,i3) = Xstar*flux(HY_MASS,i1,i2,i3)
#endif
end subroutine mscalarFluxes
 
subroutine setHSE(i1,i2,i3,dir)
   ! Note there there are data dependencies in the i1 direction.
   implicit none
   integer, intent(IN) :: i1,i2,i3,dir
   real :: dx
   real :: accelM, accelP
   !$omp declare target
#ifndef GRAVITY
   return
#else
   dx = del(dir)
#ifdef GPOT_VAR
   accelM = grv(i1)   - grv(i1-1)
   accelP = grv(i1+1) - grv(i1)
#else
   accelM = 0.5*(grv(i1-1) + grv(i1)  )*dx
   accelP = 0.5*(grv(i1)   + grv(i1+1))*dx
#endif
   snake(HY_PRES,i1-1,i2+0,i3+0) = snake(HY_PRES,i1-1,i2+0,i3+0) - &
         0.5*(snake(HY_DENS,i1-1,i2+0,i3+0)+snake(HY_DENS,i1,i2,i3)) * accelM
   snake(HY_PRES,i1+1,i2+0,i3+0) = snake(HY_PRES,i1+1,i2+0,i3+0) + &
         0.5*(snake(HY_DENS,i1,i2,i3)+snake(HY_DENS,i1+1,i2+0,i3+0)) * accelP
#endif /* GRAVITY */
end subroutine setHSE

subroutine resetHSE(i1,i2,i3,dir)
   implicit none
   integer, intent(IN) :: i1,i2,i3,dir
   real :: dx
   real :: accelM, accelP
   !$omp declare target
#ifndef GRAVITY
   return
#else
   dx = del(dir)
#ifdef GPOT_VAR
   accelM = grv(i1)   - grv(i1-1)
   accelP = grv(i1+1) - grv(i1)
#else
   accelM = 0.5*(grv(i1-1) + grv(i1)  )*dx
   accelP = 0.5*(grv(i1)   + grv(i1+1))*dx
#endif
   snake(HY_PRES,i1-1,i2+0,i3+0) = snake(HY_PRES,i1-1,i2+0,i3+0) + &
         0.5*(snake(HY_DENS,i1-1,i2+0,i3+0)+snake(HY_DENS,i1,i2,i3)) * accelM
   snake(HY_PRES,i1+1,i2+0,i3+0) = snake(HY_PRES,i1+1,i2+0,i3+0) - &
         0.5*(snake(HY_DENS,i1,i2,i3)+snake(HY_DENS,i1+1,i2+0,i3+0)) * accelP
#endif /* GRAVITY */
end subroutine resetHSE
 
!! Artificial viscosity as in Colella and Woodward
!! Calculation of cvisc here is different from the UHD solver
!! I am not including the transverse velocity terms to keep everything
!! neat and 1D for the pencil data
subroutine avisc(i1,i2,i3,dir,nvars) !flux,V,i1,nvars,dir)
   implicit none
   integer, intent(IN) :: nvars
   integer, intent(IN) :: i1,i2,i3,dir
   real :: cvisc, VenerLo, VenerHi
   !$omp declare target
   cvisc = hy_cvisc*max(-(snake(HY_VELX+dir-1,i1,i2,i3) - snake(HY_VELX+dir-1,i1-1,i2+0,i3+0)),0.)

   ! Construct minus and plus TOTAL energy densities
   VenerLo = snake(HY_DENS,i1-1,i2+0,i3+0)*0.5* &
            dot_product(snake(HY_VELX:HY_VELZ,i1-1,i2+0,i3+0), &
                        snake(HY_VELX:HY_VELZ,i1-1,i2+0,i3+0)) &
            + snake(HY_RHOE,i1-1,i2+0,i3+0)
   VenerHi = snake(HY_DENS,i1,i2,i3)*0.5* &
            dot_product(snake(HY_VELX:HY_VELZ,i1,i2,i3), &
                        snake(HY_VELX:HY_VELZ,i1,i2,i3)) &
            + snake(HY_RHOE,i1,i2,i3)

   flux(HY_MASS:HY_ENER,i1,i2,i3) = &
         flux(HY_MASS:HY_ENER,i1,i2,i3) &
         +cvisc*(/snake(HY_DENS,i1-1,i2+0,i3+0)                 - snake(HY_DENS,i1,i2,i3)&
         ,        snake(HY_DENS,i1-1,i2+0,i3+0) &
                * snake(HY_VELX,i1-1,i2+0,i3+0) &
                - snake(HY_DENS,i1,i2,i3)*snake(HY_VELX,i1,i2,i3)&
         ,        snake(HY_DENS,i1-1,i2+0,i3+0) &
                * snake(HY_VELY,i1-1,i2+0,i3+0) &
                - snake(HY_DENS,i1,i2,i3)*snake(HY_VELY,i1,i2,i3)&
         ,        snake(HY_DENS,i1-1,i2+0,i3+0) &
                * snake(HY_VELZ,i1-1,i2+0,i3+0) &
                - snake(HY_DENS,i1,i2,i3)*snake(HY_VELZ,i1,i2,i3)&
         ,        VenerLo                         - VenerHi/)
#ifdef SPARK_GLM
   flux(HY_FMGX:HY_FPSI,i1,i2,i3) = &
         flux(HY_FMGX:HY_FPSI,i1,i2,i3) &
         +cvisc*(/snake(HY_MAGX,i1-1,i2+0,i3+0)                 - snake(HY_MAGX,i1,i2,i3)&
         ,        snake(HY_MAGY,i1-1,i2+0,i3+0)                 - snake(HY_MAGY,i1,i2,i3)&
         ,        snake(HY_MAGZ,i1-1,i2+0,i3+0)                 - snake(HY_MAGZ,i1,i2,i3)&
         ,        snake(HY_PSIB,i1-1,i2+0,i3+0)                 - snake(HY_PSIB,i1,i2,i3)/)
#endif
end subroutine avisc
 
! !~ Flattening has not been tested yet in FLASH5, only 1D & 2D runs so far.
! subroutine flattening(flat3d,limits)
!    !! This follows Miller & Colella 2002
!    use Hydro_data, ONLY : hy_starState
!    implicit none
!    integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits
!    !real, intent(OUT) :: flat3d(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
!    !real :: flatTilde(NDIM,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
!    real, intent(OUT) :: flat3d(xLoGC:xHiGC,yLoGC:yHiGC,zLoGC:zHiGC)
!    real :: flatTilde(NDIM,xLoGC:xHiGC,yLoGC:yHiGC,zLoGC:zHiGC)
!    real :: beta, Z
!    real, parameter :: betaMin = 0.75, betaMax = 0.85
!    real, parameter :: Zmin = 0.25, Zmax = 0.75
!    integer :: i,j,k, kx, ky, kz
!       pLo = blkLimitsGC(LOW,maxloc(blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM),1))
!    pHi = blkLimitsGC(HIGH,maxloc(blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM),1))
!    kx = 1
!    ky = 0
!    kz = 0
! #if NDIM>1
!    ky = 1
! #if NDIM==3
!    kz = 1
! #endif
! #endif
!    do k=limits(LOW,KAXIS)-kz, limits(HIGH,KAXIS)+kz
!       do j=limits(LOW,JAXIS)-ky, limits(HIGH,JAXIS)+ky
!          do i=limits(LOW,IAXIS)-kx, limits(HIGH,IAXIS)+kx
!             !1e-99 used to be TINY(1.0) but that produced Seg Faults
!             beta = abs(hy_starState(PRES_VAR,i+1,j,k)-hy_starState(PRES_VAR,i-1,j,k)) &
!                   / max(1e-99, abs(hy_starState(PRES_VAR,i+2,j,k)-hy_starState(PRES_VAR,i-2,j,k)))
!             Z    = abs(hy_starState(PRES_VAR,i+1,j,k)-hy_starState(PRES_VAR,i-1,j,k)) &
!                   / (hy_starState(GAMC_VAR,i,j,k)*hy_starState(PRES_VAR,i,j,k))
!             flatTilde(IAXIS,i,j,k) = max(0.,min(1.,(betaMax-beta)/(betaMax-betaMin)))
!             if (hy_starState(VELX_VAR,i+1,j,k)<hy_starState(VELX_VAR,i,j,k)) then
!                flatTilde(IAXIS,i,j,k) = max(flatTilde(IAXIS,i,j,k), &
!                      &                   min(1., (Zmax-Z)/(Zmax-Zmin)))
!             else
!                flatTilde(IAXIS,i,j,k) = 1.0
!             end if
! #if NDIM>1
!             beta = abs(hy_starState(PRES_VAR,i,j+1,k)-hy_starState(PRES_VAR,i,j-1,k)) &
!                   / max(1.e-99,abs(hy_starState(PRES_VAR,i,j+2,k)-hy_starState(PRES_VAR,i,j-2,k)))
!             Z    = abs(hy_starState(PRES_VAR,i,j+1,k)-hy_starState(PRES_VAR,i,j-1,k)) &
!                   / (hy_starState(GAMC_VAR,i,j,k)*hy_starState(PRES_VAR,i,j,k))
!             flatTilde(JAXIS,i,j,k) = max(0.,min(1.,(betaMax-beta)/(betaMax-betaMin)))
!             if (hy_starState(VELY_VAR,i,j+1,k)<hy_starState(VELY_VAR,i,j,k)) then
!                flatTilde(JAXIS,i,j,k) = max(flatTilde(JAXIS,i,j,k), &
!                      &                   min(1., (Zmax-Z)/(Zmax-Zmin)))
!             else
!                flatTilde(JAXIS,i,j,k) = 1.0
!             end if
! #if NDIM==3
!             beta = abs(hy_starState(PRES_VAR,i,j,k+1)-hy_starState(PRES_VAR,i,j,k-1)) &
!                   / max(1.e-99,abs(hy_starState(PRES_VAR,i,j,k+2)-hy_starState(PRES_VAR,i,j,k-2)))
!             Z    = abs(hy_starState(PRES_VAR,i,j,k+1)-hy_starState(PRES_VAR,i,j,k-1)) &
!                   / (hy_starState(GAMC_VAR,i,j,k)*hy_starState(PRES_VAR,i,j,k))
!             flatTilde(KAXIS,i,j,k) = max(0.,min(1.,(betaMax-beta)/(betaMax-betaMin)))
!             if (hy_starState(VELZ_VAR,i,j,k+1)<hy_starState(VELZ_VAR,i,j,k)) then
!                flatTilde(KAXIS,i,j,k) = max(flatTilde(KAXIS,i,j,k), &
!                      &                   min(1., (Zmax-Z)/(Zmax-Zmin)))
!             else
!                flatTilde(KAXIS,i,j,k) = 1.0
!             end if
! #endif
! #endif
!          end do
!       end do
!    end do
!    do k=limits(LOW,KAXIS)-kz, limits(HIGH,KAXIS)+kz
!       do j=limits(LOW,JAXIS)-ky, limits(HIGH,JAXIS)+ky
!          do i=limits(LOW,IAXIS)-kx, limits(HIGH,IAXIS)+kx
!             flat3d(i,j,k) = minval(flatTilde(1:NDIM,i,j,k))
! #ifdef FLAT_VAR
!             hy_starState(FLAT_VAR,i,j,k) = flat3d(i,j,k)
! #endif
!          end do
!       end do
!    end do
! end subroutine flattening
 
subroutine ensurePhysicalState(state)
   implicit none
   real, dimension(NRECON), intent(INOUT), target  :: state
   integer :: s
   real :: spcSumInv
   real, pointer :: spc(:)
   !$omp declare target
   state(HY_DENS) = max(hy_smalldens, state(HY_DENS))
   state(HY_PRES) = max(hy_smallpres, state(HY_PRES))
#if NSPECIES>0
   ! Limit and renormalize the species.
   spc => state(HY_NUM_VARS+1:HY_NUM_VARS+NSPECIES)
   do s = 1, NSPECIES
      spc(s) = max(hy_smallX,min(1.0,spc(s)))
   end do
   spcSumInv = 1./sum(spc(1:NSPECIES))
   spc = spc*spcSumInv
#endif
end subroutine ensurePhysicalState
 
end subroutine hy_rk_getFaceFlux
 
 subroutine get_scratch_indices(dirLims, i_new, j_new, k_new, i, j, k, guardCells)
   implicit none
   integer ,intent(in) :: i, j ,k
   integer, dimension(3) :: guardCells
   integer ,intent(in), dimension(LOW:HIGH,MDIM)  ::  dirLims
   integer ,intent(OUT) :: i_new, j_new, k_new
   !$omp declare target
   i_new = 1 + i - dirLims(LOW,1) + guardCells(1)
   j_new = 1 + j - dirLims(LOW,2) + guardCells(2)
   k_new = 1 + k - dirLims(LOW,3) + guardCells(3)
 
 end subroutine get_scratch_indices

 subroutine initializeSnake(dir,dirLims,guardCells,snake,grv,shck,flat, &
                            hy_starState,flat3d,loGC)
implicit none
integer, intent(IN) :: dir
integer :: l,i,j,k,n
integer, dimension(3), intent(IN) :: guardCells, loGC
integer :: i_s, j_s, k_s   ! Scratch indices
integer, dimension(LOW:HIGH,MDIM), intent(IN) :: dirLims
real, dimension(:,:,:,:), intent(INOUT) :: snake
real, dimension(:,:,:), intent(INOUT) :: grv, shck, flat
real, dimension(1:,loGC(1):,loGC(2):,loGC(3):), intent(IN) :: hy_starState
real, dimension(loGC(1):,loGC(2):,loGC(3):), intent(IN) :: flat3d
!!$omp declare target

! Give snake proper values
! Maps the hk_solnData to a snake, such that it is contiguous in memory in the direction of iteration.
 !$omp target teams distribute parallel do collapse(3) default(none) shared(dir,hy_starState,snake,grv,shck,flat,flat3d,dirLims,guardCells) private(i,j,k,n,i_s, j_s, k_s)
   do k = dirLims(LOW,3) - guardCells(3), dirLims(HIGH,3) + guardCells(3)
      do j = dirLims(LOW,2) - guardCells(2), dirLims(HIGH,2) + guardCells(2)
         do i = dirLims(LOW,1) - guardCells(1), dirLims(HIGH,1) + guardCells(1)
            call get_scratch_indices(dirLims,i_s,j_s,k_s,i,j,k, guardCells)
            if (dir == IAXIS) then
            snake(HY_DENS,i_s,j_s,k_s) = hy_starState(DENS_VAR,i,j,k)
            snake(HY_VELX,i_s,j_s,k_s) = hy_starState(VELX_VAR,i,j,k)
            snake(HY_VELY,i_s,j_s,k_s) = hy_starState(VELY_VAR,i,j,k)
            snake(HY_VELZ,i_s,j_s,k_s) = hy_starState(VELZ_VAR,i,j,k)
            snake(HY_PRES,i_s,j_s,k_s) = hy_starState(PRES_VAR,i,j,k)
            snake(HY_GAMC,i_s,j_s,k_s) = hy_starState(GAMC_VAR,i,j,k)
            snake(HY_RHOE,i_s,j_s,k_s) = hy_starState(DENS_VAR,i,j,k)*hy_starState(EINT_VAR,i,j,k)
#ifdef SPARK_GLM
            snake(HY_MAGX,i_s,j_s,k_s) = hy_starState(MAGX_VAR,i,j,k)
            snake(HY_MAGY,i_s,j_s,k_s) = hy_starState(MAGY_VAR,i,j,k)
            snake(HY_MAGZ,i_s,j_s,k_s) = hy_starState(MAGZ_VAR,i,j,k)
            snake(HY_PSIB,i_s,j_s,k_s) = hy_starState(PSIB_VAR,i,j,k)
#endif
#if NSPECIES+NMASS_SCALARS>0
            do n=SPECIES_BEGIN, MASS_SCALARS_END
               snake(HY_NUM_VARS+1+n-SPECIES_BEGIN,i_s,j_s,k_s)    = hy_starState(n,i,j,k)
            enddo
#endif
#ifdef GRAVITY
#ifdef GPOT_VAR
            grv(i_s,j_s,k_s) = hy_starState(GPOT_VAR,i,j,k)
#else
            grv(i_s,j_s,k_s) = 0.0
#endif
#endif
#ifdef SHOK_VAR
            shck(i_s,j_s,k_s) = hy_starState(SHOK_VAR,i,j,k)
#else
            shck(i_s,j_s,k_s) = 0.0
#endif
            flat(i_s,j_s,k_s) = flat3d(i,j,k)
            
else if (dir == JAXIS) then
            snake(HY_DENS,i_s,j_s,k_s) = hy_starState(DENS_VAR,j,i,k)
            snake(HY_VELX,i_s,j_s,k_s) = hy_starState(VELX_VAR,j,i,k)
            snake(HY_VELY,i_s,j_s,k_s) = hy_starState(VELY_VAR,j,i,k)
            snake(HY_VELZ,i_s,j_s,k_s) = hy_starState(VELZ_VAR,j,i,k)
            snake(HY_PRES,i_s,j_s,k_s) = hy_starState(PRES_VAR,j,i,k)
            snake(HY_GAMC,i_s,j_s,k_s) = hy_starState(GAMC_VAR,j,i,k)
            snake(HY_RHOE,i_s,j_s,k_s) = hy_starState(DENS_VAR,j,i,k)*hy_starState(EINT_VAR,j,i,k)
#ifdef SPARK_GLM
            snake(HY_MAGX,i_s,j_s,k_s) = hy_starState(MAGX_VAR,j,i,k)
            snake(HY_MAGY,i_s,j_s,k_s) = hy_starState(MAGY_VAR,j,i,k)
            snake(HY_MAGZ,i_s,j_s,k_s) = hy_starState(MAGZ_VAR,j,i,k)
            snake(HY_PSIB,i_s,j_s,k_s) = hy_starState(PSIB_VAR,j,i,k)
#endif
#if NSPECIES+NMASS_SCALARS>0
            do n=SPECIES_BEGIN, MASS_SCALARS_END
               snake(HY_NUM_VARS+1+n-SPECIES_BEGIN,i_s,j_s,k_s)    = hy_starState(n,j,i,k)
            enddo
#endif
#ifdef GRAVITY
#ifdef GPOT_VAR
            grv(i_s,j_s,k_s) = hy_starState(GPOT_VAR,j,i,k)
#else
            grv(i_s,j_s,k_s) = 0.0
#endif
#endif
#ifdef SHOK_VAR
            shck(i_s,j_s,k_s) = hy_starState(SHOK_VAR,j,i,k)
#else
            shck(i_s,j_s,k_s) = 0.0
#endif
            flat(i_s,j_s,k_s) = flat3d(j,i,k)

else if (dir == KAXIS) then
            snake(HY_DENS,i_s,j_s,k_s) = hy_starState(DENS_VAR,j,k,i)
            snake(HY_VELX,i_s,j_s,k_s) = hy_starState(VELX_VAR,j,k,i)
            snake(HY_VELY,i_s,j_s,k_s) = hy_starState(VELY_VAR,j,k,i)
            snake(HY_VELZ,i_s,j_s,k_s) = hy_starState(VELZ_VAR,j,k,i)
            snake(HY_PRES,i_s,j_s,k_s) = hy_starState(PRES_VAR,j,k,i)
            snake(HY_GAMC,i_s,j_s,k_s) = hy_starState(GAMC_VAR,j,k,i)
            snake(HY_RHOE,i_s,j_s,k_s) = hy_starState(DENS_VAR,j,k,i)*hy_starState(EINT_VAR,j,k,i)
#ifdef SPARK_GLM
            snake(HY_MAGX,i_s,j_s,k_s) = hy_starState(MAGX_VAR,j,k,i)
            snake(HY_MAGY,i_s,j_s,k_s) = hy_starState(MAGY_VAR,j,k,i)
            snake(HY_MAGZ,i_s,j_s,k_s) = hy_starState(MAGZ_VAR,j,k,i)
            snake(HY_PSIB,i_s,j_s,k_s) = hy_starState(PSIB_VAR,j,k,i)
#endif
#if NSPECIES+NMASS_SCALARS>0
            do n=SPECIES_BEGIN, MASS_SCALARS_END
               snake(HY_NUM_VARS+1+n-SPECIES_BEGIN,i_s,j_s,k_s)    = hy_starState(n,j,k,i)
            enddo
#endif
#ifdef GRAVITY
#ifdef GPOT_VAR
            grv(i_s,j_s,k_s) = hy_starState(GPOT_VAR,j,k,i)
#else
            grv(i_s,j_s,k_s) = 0.0
#endif
#endif
#ifdef SHOK_VAR
            shck(i_s,j_s,k_s) = hy_starState(SHOK_VAR,j,k,i)
#else
            shck(i_s,j_s,k_s) = 0.0
#endif
            flat(i_s,j_s,k_s) = flat3d(j,k,i)
         endif

         enddo
      enddo
   enddo



end subroutine initializeSnake


#include "reconstruct.F90"
#include "riemann.F90"
