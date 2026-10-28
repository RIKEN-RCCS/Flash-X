#include "Simulation.h"
! No Grid operations or mutable module configuration in this numerical driver.
module hy_gpu_amr_advance_module
  use hy_gpu_amr_kernels
  use hy_gpu_amr_timing
  implicit none
contains
  subroutine hy_gpu_amr_advance(state,reference,minus,plus,flat,flux,buffers,nx,ny,nz,nb,ng,stencil, &
      dx,dy,dz,dt,coefficients,weights,params,energyFloor,gamma,eosFloor,eintSwitch,status,failedStage,failedOperation,profile)
    type(hy_amr_profile),intent(inout)::profile
    real(8)::phaseStamp,advanceStamp
    integer,intent(in)::nx,ny,nz,nb,ng,stencil
    real(rk),intent(inout)::state(1-ng:nx+ng,1-ng:ny+ng,1-(NDIM-2)*ng:nz+(NDIM-2)*ng,9,nb)
    real(rk),intent(inout)::reference(1-ng:nx+ng,1-ng:ny+ng,1-(NDIM-2)*ng:nz+(NDIM-2)*ng,9,nb)
    real(rk),intent(inout)::minus(1-ng:nx+ng,1-ng:ny+ng,1-(NDIM-2)*ng:nz+(NDIM-2)*ng,7,nb),plus(1-ng:nx+ng,1-ng:ny+ng, &
  & 1-(NDIM-2)*ng:nz+(NDIM-2)*ng,7,nb)
    real(rk),intent(inout)::flat(1-ng:nx+ng,1-ng:ny+ng,1-(NDIM-2)*ng:nz+(NDIM-2)*ng,nb)
    real(rk),intent(inout)::flux(1-ng:nx+ng+1,1-ng:ny+ng+1,1-(NDIM-2)*ng:nz+(NDIM-2)*ng+(NDIM-2),5,NDIM,nb),buffers( &
  & 1:nx+1,1:ny+1,1:nz+(NDIM-2),5,NDIM,nb)
    real(rk),intent(in)::dx,dy,dz,dt,coefficients(3,2),weights(2),energyFloor,gamma,eosFloor,eintSwitch
    type(hydro_parameters),intent(in)::params
    integer,intent(out)::status,failedStage
    character(*),intent(out)::failedOperation
    integer::stage,lower(3),upper(3),i,j,k,b,v,d,ex,ey,ez,error,diagnosticError
    real(rk)::ekin,internal,candidate,wt
    status=HY_BAD_INPUT; failedStage=0; failedOperation='input dimensions/guard width'
    if(min(nx,ny,nz,nb)<1.or.stencil<3.or.ng<stencil+3) return
    failedOperation='nonfinite input parameters'
    if(.not.all(ieee_is_finite([dx,dy,dz,dt,energyFloor,gamma,eosFloor,eintSwitch]))) return
    failedOperation='input spacing/timestep/EOS parameter range'
    if(min(dx,dy,dz,dt)<=0.or.gamma<=1.or.min(energyFloor,eosFloor,eintSwitch)<0) return
    failedOperation='nonfinite RK coefficients/weights'
    if(.not.all(ieee_is_finite(coefficients)).or..not.all(ieee_is_finite(weights))) return
    call profile_start(P_ADVANCE,advanceStamp)
    ! Owner has already allocated persistent device buffers for the maximum batch.
    !$omp target data map(alloc:state,reference,minus,plus,flat,flux,buffers)
    call profile_start(P_TO,phaseStamp)
    !$omp target update to(state)
    call profile_stop(P_TO,phaseStamp,profile)
    call profile_start(P_SAVE,phaseStamp)
    !$omp target teams distribute parallel do collapse(4) private(v)
    do b=1,nb
      do k=1-(NDIM-2)*ng,nz+(NDIM-2)*ng
      do j=1-ng,ny+ng
        do i=1-ng,nx+ng
          do v=1,9
            reference(i,j,k,v,b)=state(i,j,k,v,b)
          end do
        end do
      end do
    end do
      end do
    !$omp end target teams distribute parallel do
    call profile_stop(P_SAVE,phaseStamp,profile)
    do stage=1,2
      failedStage=stage
      lower=1; upper=[nx,ny,nz]
      lower(1:NDIM)=1-(2-stage)*stencil
      upper(1:NDIM)=upper(1:NDIM)+(2-stage)*stencil
      failedOperation='flattening'
      call profile_start(P_FLAT,phaseStamp)
      call flatten(state,flat,nx,ny,nz,nb,ng,lower,upper,params,error)
      call profile_stop(P_FLAT,phaseStamp,profile)
      if(error/=HY_OK) then
        status=error; exit
      end if
      do d=1,NDIM
        write(failedOperation,'(A,I0)') 'reconstruction direction ',d
        call profile_start(P_RECON,phaseStamp)
        call reconstruct(state,flat,minus,plus,nx,ny,nz,nb,ng,lower,upper,params,d,error)
        call profile_stop(P_RECON,phaseStamp,profile)
        if(error/=HY_OK) exit
        write(failedOperation,'(A,I0)') 'Riemann flux direction ',d
        call profile_start(P_RIEMANN,phaseStamp)
        call faces(state,minus,plus,flux,nx,ny,nz,nb,ng,lower,upper,params,d,error)
        call profile_stop(P_RIEMANN,phaseStamp,profile)
        if(error/=HY_OK) then
          !$omp target update from(state,minus,plus)
          call diagnose_faces(state,minus,plus,flux,nx,ny,nz,nb,ng,lower,upper,params,d,diagnosticError)
          if(diagnosticError==HY_OK) print *, 'HY_GPU_AMR: face failure did not reproduce in host replay'
          exit
        end if
      end do
      if(error/=HY_OK) then
        status=error; exit
      end if
      failedOperation='solution update'
      call profile_start(P_UPDATE,phaseStamp)
      call update(reference,state,flux,nx,ny,nz,nb,ng,lower,upper,dx,dy,dz,dt, &
                  coefficients(:,stage),params%density_floor,energyFloor,error)
      call profile_stop(P_UPDATE,phaseStamp,profile)
      if(error/=HY_OK) then
        status=error; exit
      end if
      ! Flash-X telescoping SPFS requires weighted flux DENSITIES: no dt or area.
      call profile_start(P_FLUX,phaseStamp)
      wt=weights(stage)
      do d=1,NDIM
        ex=0; ey=0; ez=0
        if(d==1) ex=1
        if(d==2) ey=1
        if(d==3) ez=1
        !$omp target teams distribute parallel do collapse(4) private(v)
        do b=1,nb
          do k=1,nz+ez
            do j=1,ny+ey
              do i=1,nx+ex
                do v=1,5
                  if(stage==1) buffers(i,j,k,v,d,b)=0
                  buffers(i,j,k,v,d,b)=buffers(i,j,k,v,d,b)+wt*flux(i,j,k,v,d,b)
                end do
              end do
            end do
          end do
        end do
        !$omp end target teams distribute parallel do
      end do
      call profile_stop(P_FLUX,phaseStamp,profile)
      ! MODE_DENS_EI preprocessing/postprocessing from Eos.ini, then gamma-law p.
      ! Temperature/GAME/entropy are restored by the host EOS after unpacking.
      call profile_start(P_EOS,phaseStamp)
      failedOperation='gamma-law EOS recovery'
      error=HY_OK
      !$omp target teams distribute parallel do collapse(4) private(ekin,internal,candidate) reduction(max:error)
      do b=1,nb
        do k=lower(3),upper(3)
        do j=lower(2),upper(2)
          do i=lower(1),upper(1)
            ekin=0.5_rk*sum(state(i,j,k,2:4,b)**2)
            internal=state(i,j,k,7,b); candidate=state(i,j,k,8,b)-ekin
            if(candidate>max(eosFloor,eintSwitch*ekin)) internal=candidate
            internal=max(internal,eosFloor)
            state(i,j,k,7,b)=internal; state(i,j,k,8,b)=internal+ekin
            state(i,j,k,5,b)=state(i,j,k,1,b)*internal*(gamma-1._rk)
            state(i,j,k,6,b)=gamma
            if(.not.all(ieee_is_finite(state(i,j,k,1:8,b)))) error=HY_BAD_STATE
          end do
        end do
      end do
        end do
      !$omp end target teams distribute parallel do
      call profile_stop(P_EOS,phaseStamp,profile)
      if(error/=HY_OK) then
        status=error; exit
      end if
      status=HY_OK
    end do
    if(status==HY_OK) then
      call profile_start(P_FROM,phaseStamp)
      !$omp target update from(state,buffers)
      call profile_stop(P_FROM,phaseStamp,profile)
      failedStage=0; failedOperation='none'
    end if
    !$omp end target data
    call profile_stop(P_ADVANCE,advanceStamp,profile)
  end subroutine
end module
