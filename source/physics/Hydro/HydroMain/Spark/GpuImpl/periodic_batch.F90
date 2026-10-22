! Uniform 1D periodic batched baseline. No host accesses to resident arrays inside timesteps.
program periodic_batch
  use gpu_batch
  use omp_lib
  implicit none
  integer::n,nb,ncells,b,i,v,g,donor,db,di,stage,step,error,hostFlag,ios,unit
  real(rk),allocatable::state(:,:,:),reference(:,:,:),minus(:,:,:),plus(:,:,:),flux(:,:,:),integrated(:,:,:),flat(:,:),sound(:,:)
  real(rk)::dx,dt,time,finish,rate,speed,a,mp,ms,dv,pg,marker,rho,p,ei,x
  real(rk)::coeff(3,3),weights(3),mass,energy,momentum,t0,t1,imbalance
  integer::status
  type(hydro_parameters)::params
  character(256)::arg,path
  ncells=256; nb=4; finish=0.05_rk; path='gpu_sod.csv'
  call get_command_argument(1,arg)
  if(len_trim(arg)>0) read(arg,*) ncells
  call get_command_argument(2,arg)
  if(len_trim(arg)>0) read(arg,*) nb
  call get_command_argument(3,arg)
  if(len_trim(arg)>0) read(arg,*) finish
  call get_command_argument(4,arg)
  if(len_trim(arg)>0) path=arg
  if(nb<1.or.ncells<8.or.finish<=0) stop 1
  if(mod(ncells,nb)/=0) stop 1
  n=ncells/nb; dx=1._rk/ncells
  ! Refuse silent host fallback. Set GPU_ALLOW_HOST=1 solely for correctness development.
  call get_environment_variable('GPU_ALLOW_HOST',arg,status=ios)
  hostFlag=1
  !$omp target map(from:hostFlag)
  if(.not.omp_is_initial_device()) hostFlag=0
  !$omp end target
  if(hostFlag==1.and.trim(arg)/='1') then
    print *, 'No target device: use OMP_TARGET_OFFLOAD=MANDATORY on RIKYU'
    stop 2
  end if
  print *, 'initial_device=',hostFlag,' target_devices=',omp_get_num_devices()
  allocate(state(-2:n+3,9,nb),reference(-2:n+3,9,nb),minus(-2:n+3,7,nb),plus(-2:n+3,7,nb))
  allocate(flux(n+1,5,nb),integrated(n+1,5,nb),flat(-2:n+3,nb),sound(-2:n+3,nb))
  state=0
  do b=1,nb
    do i=1,n
      x=((b-1)*n+i-0.5_rk)*dx
      rho=0.125_rk; p=0.1_rk
      if(x<0.5_rk) then
        rho=1; p=1
      end if
      ei=p/(0.4_rk*rho)
      state(i,1,b)=rho; state(i,5,b)=p; state(i,6,b)=1.4_rk
      state(i,7,b)=ei; state(i,8,b)=ei
    end do
  end do
  coeff(:,1)=[1._rk,0._rk,1._rk]
  coeff(:,2)=[0.75_rk,0.25_rk,0.25_rk]
  coeff(:,3)=[1._rk/3,2._rk/3,2._rk/3]
  weights=[1._rk/6,1._rk/6,2._rk/3]
  params%hybrid=.true.; params%shock_marker_present=.true.; params%viscosity=0.1_rk
  time=0; step=0
  t0=omp_get_wtime()
  !$omp target data map(tofrom:state) map(alloc:reference,minus,plus,flux,integrated,flat,sound)
  t1=omp_get_wtime()
  do while(time<finish)
    call guards()
    call shocks()
    call timestep_rate(state,n,nb,dx,rate,error)
    if(error/=0.or.rate<=0) stop 3
    dt=min(0.4_rk/rate,finish-time)
    !$omp target teams distribute parallel do collapse(2) private(v)
    do b=1,nb
      do i=-2,n+3
        do v=1,9
          reference(i,v,b)=state(i,v,b)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    do stage=1,3
      call guards()
      call flatten(state,flat,n,nb,params,error)
      if(error/=0) stop 4
      call reconstruct(state,flat,minus,plus,n,nb,params,error)
      call faces(state,minus,plus,flux,n,nb,params,error)
      if(error/=0) stop 5
      ! Negative boundary face owns each shared periodic transfer.
      !$omp target teams distribute parallel do collapse(2) private(db)
      do b=1,nb
        do v=1,5
          db=mod(b,nb)+1
          flux(n+1,v,b)=flux(1,v,db)
        end do
      end do
      !$omp end target teams distribute parallel do
      call update(reference,state,flux,integrated,n,nb,dx,dt,coeff(:,stage),weights(stage),stage==1,error)
      if(error/=0) stop 6
      !$omp target teams distribute parallel do collapse(2)
      do b=1,nb
        do i=1,n
          state(i,5,b)=(state(i,6,b)-1)*state(i,1,b)*state(i,7,b)
        end do
      end do
      !$omp end target teams distribute parallel do
    end do
    if(time+dt<=time) stop 7
    time=time+dt; step=step+1
    if(step>10000000) stop 7
  end do
  print *, 'resident stepping seconds=',omp_get_wtime()-t1
  !$omp end target data
  print *, 'including transfer seconds=',omp_get_wtime()-t0
  mass=0; momentum=0; energy=0
  open(newunit=unit,file=trim(path),status='replace')
  write(unit,'(A)') 'x,density,velocity,pressure,internal_energy,total_energy'
  do b=1,nb
    do i=1,n
      x=((b-1)*n+i-0.5_rk)*dx
      mass=mass+dx*state(i,1,b)
      momentum=momentum+dx*state(i,1,b)*state(i,2,b)
      energy=energy+dx*state(i,1,b)*state(i,8,b)
      write(unit,'(ES24.16,5(",",ES24.16))') x,state(i,1,b),state(i,2,b),state(i,5,b),state(i,7,b),state(i,8,b)
    end do
  end do
  close(unit)
  print *, 'steps,time=',step,time
  print *, 'mass,momentum,energy errors=',mass-0.5625_rk,momentum,energy-1.375_rk
  if(maxval(abs([mass-0.5625_rk,momentum,energy-1.375_rk]))>1.e-9_rk) stop 8
contains
  subroutine guards()
    integer::bb,ii,vv,gg,dd,donorb
    !$omp target teams distribute parallel do collapse(2) private(vv,gg,dd,donorb)
    do bb=1,nb
      do ii=-2,n+3
        if(ii>=1.and.ii<=n) cycle
        gg=modulo((bb-1)*n+ii-1,ncells)
        donorb=gg/n+1; dd=mod(gg,n)+1
        do vv=1,9
          state(ii,vv,bb)=state(dd,vv,donorb)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine

  subroutine shocks()
    ! Expanded unchanged acoustic/shock helpers are provided by the generated kernel module.
    call batch_shocks(state,sound,n,nb,params%tiny_state,error)
    if(error/=0) stop 9
    call guards()
  end subroutine
end program
