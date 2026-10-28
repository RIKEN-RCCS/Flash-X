#include "Simulation.h"
#include "constants.h"
! Optional verification output. All access is on host, after the Hydro reflux pass.
module hy_gpu_amr_dump_module
  use Grid_interface, only: Grid_getTileIterator,Grid_releaseTileIterator
  use Grid_iterator, only: Grid_iterator_t
  use Grid_tile, only: Grid_tile_t
  use Driver_interface, only: Driver_getNStep,Driver_getMype,Driver_abort
  use RuntimeParameters_interface, only: RuntimeParameters_get
  implicit none
  private
  public :: hy_gpu_amr_dump
contains
  subroutine hy_gpu_amr_dump(time,dt)
    real,intent(in)::time,dt
    type(Grid_iterator_t)::iterator
    type(Grid_tile_t)::tile
    real,pointer::u(:,:,:,:)=>null()
    real::box(2,3),spacing(3),position(3),volume,totals(5),rho
    integer::interval,step,rank,unit,ios,i,j,k,lo(3),hi(3),cells
    character(256)::prefix
    character(512)::filename
    call RuntimeParameters_get('hy_gpuDumpInterval',interval)
    if(interval<=0) return
    call Driver_getNStep(step)
    if(mod(step,interval)/=0) return
    call RuntimeParameters_get('hy_gpuDumpPrefix',prefix)
    call Driver_getMype(MESH_COMM,rank)
    write(filename,'(A,"_step",I10.10,"_rank",I6.6,".csv")') trim(prefix),step,rank
    open(newunit=unit,file=trim(filename),status='replace',action='write',iostat=ios)
    if(ios/=0) call Driver_abort('HY_GPU_AMR: cannot open verification dump; check prefix directory')
    write(unit,'(A)') 'step,time,dt,level,i,j,k,x,y,z,dx,dy,dz,density,velocity_x,velocity_y,' // &
      'velocity_z,pressure,internal_energy,total_energy'
    totals=0; cells=0
    call Grid_getTileIterator(iterator,LEAF,tiling=.false.)
    do while(iterator%isValid())
      call iterator%currentTile(tile)
      lo=tile%limits(LOW,:); hi=tile%limits(HIGH,:)
      call tile%boundBox(box)
      call tile%deltas(spacing)
#if NDIM == 2
      box(:,3)=0; spacing(3)=1
#endif
      volume=product(spacing(1:NDIM))
      call tile%getDataPtr(u,CENTER)
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            position=box(LOW,:)+([i,j,k]-lo+0.5)*spacing
#if NDIM == 2
            position(3)=0
#endif
            rho=u(DENS_VAR,i,j,k)
            totals(1)=totals(1)+volume*rho
            totals(2:4)=totals(2:4)+volume*rho*[u(VELX_VAR,i,j,k),u(VELY_VAR,i,j,k),u(VELZ_VAR,i,j,k)]
            totals(5)=totals(5)+volume*rho*u(ENER_VAR,i,j,k)
            cells=cells+1
            write(unit,'(I0,2(",",ES26.17E3),4(",",I0),13(",",ES26.17E3))',iostat=ios) &
              step,time,dt,tile%level,i,j,k,position,spacing,rho,u(VELX_VAR,i,j,k),u(VELY_VAR,i,j,k), &
              u(VELZ_VAR,i,j,k),u(PRES_VAR,i,j,k),u(EINT_VAR,i,j,k),u(ENER_VAR,i,j,k)
            if(ios/=0) call Driver_abort('HY_GPU_AMR: error writing verification dump')
          end do
        end do
      end do
      call tile%releaseDataPtr(u,CENTER)
      call iterator%next()
    end do
    call Grid_releaseTileIterator(iterator)
    close(unit,iostat=ios)
    if(ios/=0) call Driver_abort('HY_GPU_AMR: error closing verification dump')
    print *, 'HY_GPU_AMR dump: ',trim(filename),' owned cells=',cells
    print *, 'HY_GPU_AMR rank-local mass,Px,Py,Pz,energy=',totals
  end subroutine
end module
