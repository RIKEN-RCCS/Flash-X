! No Grid operations or mutable module configuration in this numerical driver.
module hy_gpu_amr_advance_module
  use hy_gpu_amr_kernels
  implicit none
contains
  subroutine hy_gpu_amr_advance(state,reference,minus,plus,flat,flux,buffers,nx,ny,nb,ng,stencil, &
      dx,dy,dt,coefficients,weights,params,energyFloor,gamma,eosFloor,eintSwitch,status,failedStage)
    integer,intent(in)::nx,ny,nb,ng,stencil
    real(rk),intent(inout)::state(1-ng:nx+ng,1-ng:ny+ng,9,nb)
    real(rk),intent(inout)::reference(1-ng:nx+ng,1-ng:ny+ng,9,nb)
    real(rk),intent(inout)::minus(1-ng:nx+ng,1-ng:ny+ng,7,nb),plus(1-ng:nx+ng,1-ng:ny+ng,7,nb)
    real(rk),intent(inout)::flat(1-ng:nx+ng,1-ng:ny+ng,nb)
    real(rk),intent(inout)::flux(1-ng:nx+ng+1,1-ng:ny+ng+1,5,2,nb),buffers(1:nx+1,1:ny+1,5,2,nb)
    real(rk),intent(in)::dx,dy,dt,coefficients(3,2),weights(2),energyFloor,gamma,eosFloor,eintSwitch
    type(hydro_parameters),intent(in)::params
    integer,intent(out)::status,failedStage
    integer::stage,lower(2),upper(2),i,j,b,v,d,ex,ey,error
    real(rk)::ekin,internal,candidate,wt
    status=HY_BAD_INPUT; failedStage=0
    if(min(nx,ny,nb)<1.or.stencil<3.or.ng<stencil+3) return
    if(.not.all(ieee_is_finite([dx,dy,dt,energyFloor,gamma,eosFloor,eintSwitch]))) return
    if(min(dx,dy,dt)<=0.or.gamma<=1.or.min(energyFloor,eosFloor,eintSwitch)<0) return
    if(.not.all(ieee_is_finite(coefficients)).or..not.all(ieee_is_finite(weights))) return
    ! Owner has already allocated persistent device buffers for the maximum batch.
    !$omp target data map(alloc:state,reference,minus,plus,flat,flux,buffers)
    !$omp target update to(state)
    !$omp target teams distribute parallel do collapse(3) private(v)
    do b=1,nb
      do j=1-ng,ny+ng
        do i=1-ng,nx+ng
          do v=1,9
            reference(i,j,v,b)=state(i,j,v,b)
          end do
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    do stage=1,2
      failedStage=stage
      lower=1-(2-stage)*stencil; upper=[nx,ny]+(2-stage)*stencil
      call flatten(state,flat,nx,ny,nb,ng,lower,upper,params,error)
      if(error/=HY_OK) then
        status=error; exit
      end if
      do d=1,2
        call reconstruct(state,flat,minus,plus,nx,ny,nb,ng,lower,upper,params,d,error)
        call faces(state,minus,plus,flux,nx,ny,nb,ng,lower,upper,params,d,error)
        if(error/=HY_OK) exit
      end do
      if(error/=HY_OK) then
        status=error; exit
      end if
      call update(reference,state,flux,nx,ny,nb,ng,lower,upper,dx,dy,dt, &
                  coefficients(:,stage),params%density_floor,energyFloor,error)
      if(error/=HY_OK) then
        status=error; exit
      end if
      ! Flash-X telescoping SPFS requires weighted flux DENSITIES: no dt or area.
      wt=weights(stage)
      do d=1,2
        ex=0; ey=0
        if(d==1) ex=1
        if(d==2) ey=1
        !$omp target teams distribute parallel do collapse(3) private(v)
        do b=1,nb
          do j=1,ny+ey
            do i=1,nx+ex
              do v=1,5
                if(stage==1) buffers(i,j,v,d,b)=0
                buffers(i,j,v,d,b)=buffers(i,j,v,d,b)+wt*flux(i,j,v,d,b)
              end do
            end do
          end do
        end do
        !$omp end target teams distribute parallel do
      end do
      ! MODE_DENS_EI preprocessing/postprocessing from Eos.ini, then gamma-law p.
      ! Temperature/GAME/entropy are restored by the host EOS after unpacking.
      error=HY_OK
      !$omp target teams distribute parallel do collapse(3) private(ekin,internal,candidate) reduction(max:error)
      do b=1,nb
        do j=lower(2),upper(2)
          do i=lower(1),upper(1)
            ekin=0.5_rk*sum(state(i,j,2:4,b)**2)
            internal=state(i,j,7,b); candidate=state(i,j,8,b)-ekin
            if(candidate>max(eosFloor,eintSwitch*ekin)) internal=candidate
            internal=max(internal,eosFloor)
            state(i,j,7,b)=internal; state(i,j,8,b)=internal+ekin
            state(i,j,5,b)=state(i,j,1,b)*internal*(gamma-1._rk)
            state(i,j,6,b)=gamma
            if(.not.all(ieee_is_finite(state(i,j,1:8,b)))) error=HY_BAD_STATE
          end do
        end do
      end do
      !$omp end target teams distribute parallel do
      if(error/=HY_OK) then
        status=error; exit
      end if
      status=HY_OK
    end do
    if(status==HY_OK) then
      !$omp target update from(state,buffers)
      failedStage=0
    end if
    !$omp end target data
  end subroutine
end module
