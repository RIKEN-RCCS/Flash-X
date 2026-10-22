! 2D periodic Cartesian blocks. Reuse unchanged arithmetic; directional scratch is sequential.
program periodic_batch_2d
  use gpu_batch_2d
  implicit none
  integer::gx,gy,bx,by,nx,ny,nb,b,i,j,v,d,stage,step,error,hostFlag,ios,unit,stepLimit,test
  integer::gix,giy,ib,jb,neighbor
  real(rk),allocatable::state(:,:,:,:),reference(:,:,:,:),minus(:,:,:,:),plus(:,:,:,:)
  real(rk),allocatable::flux(:,:,:,:,:),integrated(:,:,:,:,:),flat(:,:,:),sound(:,:,:)
  real(rk)::dx,dy,x,y,rho,p,ei,coord,time,finish,dt,rate,coeff(3,3),weights(3),t0,t1
  real(rk)::mass0,energy0,mass,energy,mx,my,volume,elapsed,resident
  type(hydro_parameters)::params
  character(256)::arg,path
  gx=128; gy=128; bx=2; by=2; finish=0.05_rk; path='gpu2d.csv'; stepLimit=0; test=1
  ! Arguments: global nx, global ny, blocks-x, blocks-y, time, CSV, case, max steps.
  call get_command_argument(1,arg)
  if(len_trim(arg)>0) read(arg,*) gx
  call get_command_argument(2,arg)
  if(len_trim(arg)>0) read(arg,*) gy
  call get_command_argument(3,arg)
  if(len_trim(arg)>0) read(arg,*) bx
  call get_command_argument(4,arg)
  if(len_trim(arg)>0) read(arg,*) by
  call get_command_argument(5,arg)
  if(len_trim(arg)>0) read(arg,*) finish
  call get_command_argument(6,arg)
  if(len_trim(arg)>0) path=arg
  call get_command_argument(7,arg)
  if(len_trim(arg)>0) read(arg,*) test
  call get_command_argument(8,arg)
  if(len_trim(arg)>0) read(arg,*) stepLimit
  if(min(gx,gy)<8.or.min(bx,by)<1.or.finish<=0.or.stepLimit<0) stop 1
  if(mod(gx,bx)/=0.or.mod(gy,by)/=0.or.test<1.or.test>3) stop 1
  nx=gx/bx; ny=gy/by; nb=bx*by
  dx=1._rk/gx; dy=1._rk/gy; volume=dx*dy
  hostFlag=1
  !$omp target map(from:hostFlag)
  if(.not.omp_is_initial_device()) hostFlag=0
  !$omp end target
  call get_environment_variable('GPU_ALLOW_HOST',arg,status=ios)
  if(hostFlag==1.and.trim(arg)/='1') stop 2
  print *, 'initial_device=',hostFlag,' target_devices=',omp_get_num_devices()
  allocate(state(-2:nx+3,-2:ny+3,9,nb),reference(-2:nx+3,-2:ny+3,9,nb))
  allocate(minus(-2:nx+3,-2:ny+3,7,nb),plus(-2:nx+3,-2:ny+3,7,nb))
  allocate(flux(nx+1,ny+1,5,2,nb),integrated(nx+1,ny+1,5,2,nb))
  allocate(flat(-2:nx+3,-2:ny+3,nb),sound(-2:nx+3,-2:ny+3,nb))
  state=0; mass0=0; energy0=0
  do b=1,nb
    ib=mod(b-1,bx); jb=(b-1)/bx
    do j=1,ny
      do i=1,nx
        x=(ib*nx+i-0.5_rk)*dx; y=(jb*ny+j-0.5_rk)*dy
        coord=x
        if(test==2) coord=y
        if(test==3) coord=modulo(x+y,1._rk)
        rho=0.125_rk; p=0.1_rk
        if(coord<0.5_rk) then
          rho=1; p=1
        end if
        ei=p/(0.4_rk*rho)
        state(i,j,1,b)=rho; state(i,j,5,b)=p; state(i,j,6,b)=1.4_rk
        state(i,j,7,b)=ei; state(i,j,8,b)=ei
        mass0=mass0+volume*rho; energy0=energy0+volume*rho*ei
      end do
    end do
  end do
  coeff(:,1)=[1._rk,0._rk,1._rk]
  coeff(:,2)=[0.75_rk,0.25_rk,0.25_rk]
  coeff(:,3)=[1._rk/3,2._rk/3,2._rk/3]
  weights=[1._rk/6,1._rk/6,2._rk/3]
  params%hybrid=.true.; params%shock_marker_present=.true.; params%viscosity=0.1_rk
  time=0; step=0; t0=omp_get_wtime()
  !$omp target data map(tofrom:state) map(alloc:reference,minus,plus,flux,integrated,flat,sound)
  t1=omp_get_wtime()
  do while(time<finish)
    if(stepLimit>0.and.step>=stepLimit) exit
    call guards()
    call batch_shocks(state,sound,nx,ny,nb,params%tiny_state,error)
    if(error/=0) stop 3
    call guards()
    call timestep_rate(state,nx,ny,nb,dx,dy,rate,error)
    if(error/=0.or.rate<=0) stop 3
    dt=min(0.4_rk/rate,finish-time)
    !$omp target teams distribute parallel do collapse(3) private(v)
    do b=1,nb
      do j=-2,ny+3
        do i=-2,nx+3
          do v=1,9
            reference(i,j,v,b)=state(i,j,v,b)
          end do
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    do stage=1,3
      call guards()
      call flatten(state,flat,nx,ny,nb,params,error)
      if(error/=0) stop 4
      do d=1,2
        call reconstruct(state,flat,minus,plus,nx,ny,nb,params,d,error)
        call faces(state,minus,plus,flux,nx,ny,nb,params,d,error)
        if(error/=0) stop 5
      end do
      call synchronize_faces()
      call update(reference,state,flux,integrated,nx,ny,nb,dx,dy,dt,coeff(:,stage),weights(stage),stage==1,error)
      if(error/=0) stop 6
      !$omp target teams distribute parallel do collapse(3)
      do b=1,nb
        do j=1,ny
          do i=1,nx
            state(i,j,5,b)=(state(i,j,6,b)-1)*state(i,j,1,b)*state(i,j,7,b)
          end do
        end do
      end do
      !$omp end target teams distribute parallel do
    end do
    if(time+dt<=time) stop 7
    time=time+dt; step=step+1
    if(step>10000000) stop 7
  end do
  resident=omp_get_wtime()-t1
  !$omp end target data
  elapsed=omp_get_wtime()-t0
  print *, 'resident stepping seconds=',resident
  print *, 'including transfer seconds=',elapsed
  print *, 'steps,time=',step,time
  mass=0; energy=0; mx=0; my=0
  open(newunit=unit,file=trim(path),status='replace')  
  write(unit,'(A)') 'x,y,density,velocity_x,velocity_y,pressure,internal_energy,total_energy'
  ! Global row-major order is independent of decomposition.
  do giy=1,gy
    do gix=1,gx
      ib=(gix-1)/nx; jb=(giy-1)/ny; b=jb*bx+ib+1
      i=mod(gix-1,nx)+1; j=mod(giy-1,ny)+1
      rho=state(i,j,1,b)
      mass=mass+volume*rho; energy=energy+volume*rho*state(i,j,8,b)
      mx=mx+volume*rho*state(i,j,2,b); my=my+volume*rho*state(i,j,3,b)
      write(unit,'(ES24.16,7(",",ES24.16))') (gix-0.5_rk)*dx,(giy-0.5_rk)*dy, &
        rho,state(i,j,2,b),state(i,j,3,b),state(i,j,5,b),state(i,j,7,b),state(i,j,8,b)
    end do
  end do
  close(unit)
  print *, 'mass,momentum_x,momentum_y,energy errors=',mass-mass0,mx,my,energy-energy0
  if(maxval(abs([mass-mass0,mx,my,energy-energy0]))>1.e-9_rk) stop 8
contains
  subroutine guards()
    integer::bb,ii,jj,vv,globalx,globaly,donorb,donori,donorj,blockx,blocky
    !$omp target teams distribute parallel do collapse(3) &
    !$omp& private(vv,globalx,globaly,donorb,donori,donorj,blockx,blocky)
    do bb=1,nb
      do jj=-2,ny+3
        do ii=-2,nx+3
          if(ii>=1.and.ii<=nx.and.jj>=1.and.jj<=ny) cycle
          blockx=mod(bb-1,bx); blocky=(bb-1)/bx
          globalx=modulo(blockx*nx+ii-1,gx); globaly=modulo(blocky*ny+jj-1,gy)
          donorb=(globaly/ny)*bx+globalx/nx+1
          donori=mod(globalx,nx)+1; donorj=mod(globaly,ny)+1
          do vv=1,9
            state(ii,jj,vv,bb)=state(donori,donorj,vv,donorb)
          end do
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine
  subroutine synchronize_faces()
    integer::bb,jj,ii,vv,neighborx,neighbory,blockx,blocky
    !$omp target teams distribute parallel do collapse(2) private(vv,neighborx,blockx,blocky)
    do bb=1,nb
      do jj=1,ny
        blockx=mod(bb-1,bx); blocky=(bb-1)/bx
        neighborx=blocky*bx+mod(blockx+1,bx)+1
        do vv=1,5
          flux(nx+1,jj,vv,1,bb)=flux(1,jj,vv,1,neighborx)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    !$omp target teams distribute parallel do collapse(2) private(vv,neighbory,blockx,blocky)
    do bb=1,nb
      do ii=1,nx
        blockx=mod(bb-1,bx); blocky=(bb-1)/bx
        neighbory=mod(blocky+1,by)*bx+blockx+1
        do vv=1,5
          flux(ii,ny+1,vv,2,bb)=flux(ii,1,vv,2,neighbory)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine
end program
