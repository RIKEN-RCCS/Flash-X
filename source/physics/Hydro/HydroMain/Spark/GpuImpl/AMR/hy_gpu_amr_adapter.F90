#include "Simulation.h"
#include "constants.h"
#include "Spark.h"
#include "Eos.h"
#if NDIM != 2 && NDIM != 3
#error HY_GPU_AMR requires NDIM=2 or NDIM=3
#endif
#if defined(HY_RK3) || defined(SPARK_GLM) || defined(GRAVITY)
#error HY_GPU_AMR supports RK2 hydrodynamics without gravity or MHD
#endif
#if NSPECIES+NMASS_SCALARS > 0 || NSTENCIL < 3
#error HY_GPU_AMR requires WENO5 with no species or mass scalars
#endif
#if !defined(EOS_GAMMA_INCLUDED) || defined(EOS_HYBRID_INCLUDED)
#error HY_GPU_AMR requires the non-hybrid Gamma EOS
#endif
#ifdef INDEXREORDER
#error HY_GPU_AMR adapter currently requires variable-first Flash-X host arrays
#endif
#if defined(OMP_OL) || defined(OACC)
#error HY_GPU_AMR must not be combined with the legacy per-block offload path
#endif

module hy_gpu_amr_adapter
  use hy_gpu_amr_advance_module, only: hy_gpu_amr_advance,rk,hydro_parameters,HY_OK
  use hy_gpu_amr_timing
  use Grid_interface, only: Grid_getTileIterator,Grid_releaseTileIterator,Grid_putFluxData
  use Grid_tile, only: Grid_tile_t
  use Grid_iterator, only: Grid_iterator_t
  use Driver_interface, only: Driver_abort
  use Eos_interface, only: Eos_multiDim
  use hy_rk_interface, only: hy_rk_shockDetect
  use RuntimeParameters_interface, only: RuntimeParameters_get
  use omp_lib
  implicit none
  private
  public :: hy_gpu_amr_hydro,hy_gpu_amr_check
contains
  subroutine hy_gpu_amr_check()
    use Hydro_data, only: hy_telescoping,hy_fluxCorrectPerLevel,hy_geometry,hy_threadWithinBlock
    if(.not.hy_telescoping.or.hy_fluxCorrectPerLevel.or.hy_geometry/=CARTESIAN.or.hy_threadWithinBlock) &
      call Driver_abort('HY_GPU_AMR requires Cartesian telescoping and ordinary flux correction')
  end subroutine
  subroutine hy_gpu_amr_hydro(dt)
    ! Configuration is read here; numerical kernels receive all values explicitly.
    use Hydro_data, only: hy_telescoping,hy_fluxCorrectPerLevel,hy_geometry,hy_fluxCorrect, &
      hy_threadWithinBlock,hy_flattening,hy_hybridRiemann,hy_cvisc,hy_smalldens,hy_smallpres, &
      hy_smallE,hy_tiny,hy_coeffArray,hy_weights
    use Eos_data, only: eos_gamma,eos_smalle,eos_eintSwitch
    real,intent(in)::dt
    type(hy_amr_profile)::profile
    real(8)::adapterStamp,phaseStamp
    type(Grid_iterator_t)::iterator
    type(Grid_tile_t)::tile
    type(Grid_tile_t),allocatable::tiles(:)
    type(hydro_parameters)::params
    real,pointer::u(:,:,:,:)
    real,allocatable::vc(:,:,:),fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
    real(rk),allocatable::state(:,:,:,:,:),reference(:,:,:,:,:),minus(:,:,:,:,:),plus(:,:,:,:,:),flat(:,:,:,:)
    real(rk),allocatable::flux(:,:,:,:,:,:),buffers(:,:,:,:,:,:)
    real::spacing(3),batchSpacing(3)
    real(rk)::coefficients(3,2),weights(2)
    integer::capacity,count,nx,ny,nz,ng,lo(3),hi(3),lc(3),hc(3),i,j,k,v,b,ii,jj,kk,status,failedStage
    integer::stateMap(8),fluxMap(5),hostFlag,batchLevel
    logical::allocatedBatch
    character(80)::failedOperation
    character(256)::abortMessage
    call profile_start(P_TOTAL,adapterStamp)
    if(.not.hy_telescoping.or.hy_fluxCorrectPerLevel.or.hy_geometry/=CARTESIAN.or.hy_threadWithinBlock) &
      call Driver_abort('HY_GPU_AMR requires Cartesian telescoping and ordinary flux correction')
    if(NGUARD<NSTENCIL+3) call Driver_abort('HY_GPU_AMR: insufficient telescoping guards')
    if(eos_gamma<=1) call Driver_abort('HY_GPU_AMR: invalid gamma')
    if(NFLUXES/=5) call Driver_abort('HY_GPU_AMR: expected five hydrodynamic fluxes')
    call RuntimeParameters_get('hy_gpuBatchSize',capacity)
    if(capacity<1) call Driver_abort('HY_GPU_AMR: hy_gpuBatchSize must be positive')
    call profile_start(P_PROBE,phaseStamp)
    hostFlag=1
    !$omp target map(from:hostFlag)
    if(.not.omp_is_initial_device()) hostFlag=0
    !$omp end target
    call profile_stop(P_PROBE,phaseStamp,profile)
#ifdef HY_GPU_AMR_CPU_TEST
    if(hostFlag/=1) call Driver_abort('HY_GPU_AMR CPU test: set OMP_TARGET_OFFLOAD=DISABLED')
    print *, 'HY_GPU_AMR: CPU testing mode; OpenMP target regions execute on host'
#else
    if(hostFlag/=0) call Driver_abort('HY_GPU_AMR: target device unavailable')
#endif
    params%flattening=hy_flattening; params%hybrid=hy_hybridRiemann
    params%viscosity=hy_cvisc; params%density_floor=hy_smalldens; params%pressure_floor=hy_smallpres
    params%tiny_state=hy_tiny; params%shock_marker_present=.false.
#ifdef SHOK_VAR
    params%shock_marker_present=.true.
#endif
    coefficients=transpose(hy_coeffArray(1:2,:)); weights=hy_weights(1:2)
    stateMap=[DENS_VAR,VELX_VAR,VELY_VAR,VELZ_VAR,PRES_VAR,GAMC_VAR,EINT_VAR,ENER_VAR]
    fluxMap=[HY_MASS,HY_XMOM,HY_YMOM,HY_ZMOM,HY_ENER]
    allocate(tiles(capacity)); count=0; allocatedBatch=.false.; ng=NGUARD
    call Grid_getTileIterator(iterator,LEAF,tiling=.false.)
    do while(iterator%isValid())
      call iterator%currentTile(tile)
      lo=tile%limits(LOW,:); hi=tile%limits(HIGH,:)
      lc=tile%blkLimitsGC(LOW,:); hc=tile%blkLimitsGC(HIGH,:)
      call tile%deltas(spacing)
#if NDIM == 2
      spacing(3)=1.0
#endif
#if NDIM == 2
      if(hi(3)/=lo(3).or.lc(3)/=lo(3).or.hc(3)/=hi(3)) &
        call Driver_abort('HY_GPU_AMR: inactive z must contain exactly one unguarded plane')
#endif
      if(.not.allocatedBatch) then
        nx=hi(1)-lo(1)+1; ny=hi(2)-lo(2)+1; nz=hi(3)-lo(3)+1
        call allocate_batch()
      end if
      if(hi(1)-lo(1)+1/=nx.or.hi(2)-lo(2)+1/=ny.or.hi(3)-lo(3)+1/=nz) &
        call Driver_abort('HY_GPU_AMR: only equal-sized full blocks are supported')
      if(any(lc(1:NDIM)/=lo(1:NDIM)-ng).or.any(hc(1:NDIM)/=hi(1:NDIM)+ng)) &
        call Driver_abort('HY_GPU_AMR: unexpected prepared guard extent')
      if(count>0) then
        if(tile%level/=batchLevel.or.any(spacing/=batchSpacing)) call flush_batch(2)
      end if
      if(count==0) then
        batchSpacing=spacing; batchLevel=tile%level
      end if
      count=count+1; tiles(count)=tile
      call tile%getDataPtr(u,CENTER)
      ! Same CPU shock preparation as Hydro_prepBlock; device copies establish U0/U*.
      call profile_start(P_PREP,phaseStamp)
      allocate(vc(lc(1):hc(1),lc(2):hc(2),lc(3):hc(3)))
      call hy_rk_shockDetect(u,vc,tile%blkLimitsGC,lc)
      deallocate(vc)
      call profile_stop(P_PREP,phaseStamp,profile)
      call profile_start(P_PACK,phaseStamp)
      do k=lc(3),hc(3)
      do j=lc(2),hc(2)
        do i=lc(1),hc(1)
          ii=i-lo(1)+1; jj=j-lo(2)+1; kk=k-lo(3)+1
          do v=1,8
            state(ii,jj,kk,v,count)=u(stateMap(v),i,j,k)
          end do
          state(ii,jj,kk,9,count)=0
#ifdef SHOK_VAR
          state(ii,jj,kk,9,count)=u(SHOK_VAR,i,j,k)
#endif
        end do
      end do
      end do
      call tile%releaseDataPtr(u,CENTER)
      call profile_stop(P_PACK,phaseStamp,profile)
      if(count==capacity) call flush_batch(1)
      call iterator%next()
    end do
    if(count>0) call flush_batch(3)
    call Grid_releaseTileIterator(iterator)
    call profile_start(P_RELEASE,phaseStamp)
    if(allocatedBatch) then
      !$omp target exit data map(delete:state,reference,minus,plus,flat,flux,buffers)
      deallocate(state,reference,minus,plus,flat,flux,buffers)
    end if
    deallocate(tiles)
    call profile_stop(P_RELEASE,phaseStamp,profile)
    call profile_stop(P_TOTAL,adapterStamp,profile)
    call profile_write(profile,capacity)
  contains
    subroutine allocate_batch()
      real(8)::stamp
      call profile_start(P_ALLOC,stamp)
      allocate(state(1-ng:nx+ng,1-ng:ny+ng,1-(NDIM-2)*ng:nz+(NDIM-2)*ng,9,capacity),reference(1-ng:nx+ng,1-ng:ny+ng, &
  & 1-(NDIM-2)*ng:nz+(NDIM-2)*ng,9,capacity))
      allocate(minus(1-ng:nx+ng,1-ng:ny+ng,1-(NDIM-2)*ng:nz+(NDIM-2)*ng,7,capacity),plus(1-ng:nx+ng,1-ng:ny+ng,1-( &
  & NDIM-2)*ng:nz+(NDIM-2)*ng,7,capacity))
      allocate(flat(1-ng:nx+ng,1-ng:ny+ng,1-(NDIM-2)*ng:nz+(NDIM-2)*ng,capacity))
      allocate(flux(1-ng:nx+ng+1,1-ng:ny+ng+1,1-(NDIM-2)*ng:nz+(NDIM-2)*ng+(NDIM-2),5,NDIM,capacity),buffers(nx+1, &
  & ny+1,nz+(NDIM-2),5,NDIM,capacity))
      !$omp target enter data map(alloc:state,reference,minus,plus,flat,flux,buffers)
      allocatedBatch=.true.
      call profile_stop(P_ALLOC,stamp,profile)
    end subroutine
    subroutine flush_batch(reason)
      integer,intent(in)::reason
      real(8)::stamp
      integer::q,ix,jy,kz,ll(3),hh(3),lg(3)
      profile%batches=profile%batches+1
      profile%blocks=profile%blocks+count
      profile%minBatch=min(profile%minBatch,count)
      profile%maxBatch=max(profile%maxBatch,count)
      select case(reason)
      case(1)
        profile%capacityFlush=profile%capacityFlush+1
      case(2)
        profile%levelFlush=profile%levelFlush+1
      case(3)
        profile%tailFlush=profile%tailFlush+1
      end select
      ! Prefix sections are contiguous; the final incomplete batch uses count, not capacity.
      call hy_gpu_amr_advance(state(:,:,:,:,1:count),reference(:,:,:,:,1:count),minus(:,:,:,:,1:count), &
        plus(:,:,:,:,1:count),flat(:,:,:,1:count),flux(:,:,:,:,:,1:count),buffers(:,:,:,:,:,1:count), &
        nx,ny,nz,count,ng,NSTENCIL,batchSpacing(1),batchSpacing(2),batchSpacing(3),dt,coefficients,weights,params, &
        hy_smallE,eos_gamma,eos_smalle,eos_eintSwitch,status,failedStage,failedOperation,profile)
      if(status/=HY_OK) then
        print *, 'HY_GPU_AMR failure: status,stage,level,batch count=',status,failedStage,batchLevel,count
        print *, 'HY_GPU_AMR input: nx,ny,nz,guards,stencil=',nx,ny,nz,ng,NSTENCIL
        print *, 'HY_GPU_AMR input: spacing,dt,gamma,floors,switch=',batchSpacing,dt,eos_gamma, &
          hy_smallE,eos_smalle,eos_eintSwitch
        write(abortMessage,'(A,I0,A,I0,A,A,A,I0,A,I0)') &
          'HY_GPU_AMR: status=',status,', stage=',failedStage,', operation=',trim(failedOperation), &
          ', level=',batchLevel,', batch count=',count
        call Driver_abort(trim(abortMessage))
      end if
      do q=1,count
        call profile_start(P_UNPACK,stamp)
        ll=tiles(q)%limits(LOW,:); hh=tiles(q)%limits(HIGH,:); lg=tiles(q)%blkLimitsGC(LOW,:)
        call tiles(q)%getDataPtr(u,CENTER)
        do kz=1,nz
        do jy=1,ny
          do ix=1,nx
            do v=1,8
              u(stateMap(v),ll(1)+ix-1,ll(2)+jy-1,ll(3)+kz-1)=state(ix,jy,kz,v,q)
            end do
          end do
        end do
        end do
        call profile_stop(P_UNPACK,stamp,profile)
        ! Populate all mapped EOS outputs (TEMP, GAME, entropy if present) on CPU.
        call profile_start(P_HOST_EOS,stamp)
        call Eos_multiDim(MODE_DENS_EI,tiles(q)%limits,lg,u)
        call profile_stop(P_HOST_EOS,stamp,profile)
        call tiles(q)%releaseDataPtr(u,CENTER)
        if(hy_fluxCorrect) then
          call profile_start(P_PUBLISH,stamp)
          allocate(fx(NFLUXES,ll(1):hh(1)+1,ll(2):hh(2),ll(3):hh(3)))
          allocate(fy(NFLUXES,ll(1):hh(1),ll(2):hh(2)+1,ll(3):hh(3)))
          allocate(fz(NFLUXES,ll(1):hh(1),ll(2):hh(2),ll(3):hh(3)+1))
          fx=0; fy=0; fz=0
          do v=1,5
            do kz=1,nz
              do jy=1,ny
                do ix=1,nx+1
                  fx(fluxMap(v),ll(1)+ix-1,ll(2)+jy-1,ll(3)+kz-1)=buffers(ix,jy,kz,v,1,q)
                end do
              end do
              do jy=1,ny+1
                do ix=1,nx
                  fy(fluxMap(v),ll(1)+ix-1,ll(2)+jy-1,ll(3)+kz-1)=buffers(ix,jy,kz,v,2,q)
                end do
              end do
            end do
#if NDIM == 3
            do kz=1,nz+1
              do jy=1,ny
                do ix=1,nx
                  fz(fluxMap(v),ll(1)+ix-1,ll(2)+jy-1,ll(3)+kz-1)=buffers(ix,jy,kz,v,3,q)
                end do
              end do
            end do
#endif
          end do
          ! Match Hydro_math.ini: default telescoping Grid flux-buffer convention.
          call Grid_putFluxData(tiles(q),fx,fy,fz,ll)
          deallocate(fx,fy,fz)
          call profile_stop(P_PUBLISH,stamp,profile)
        end if
      end do
      count=0
    end subroutine
  end subroutine
end module
