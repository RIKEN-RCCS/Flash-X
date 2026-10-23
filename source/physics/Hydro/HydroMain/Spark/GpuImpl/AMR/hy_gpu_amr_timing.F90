#include "constants.h"
! Host-only instrumentation. Timings/counters are owned by each Hydro invocation.
module hy_gpu_amr_timing
  use Timers_interface, only: Timers_start,Timers_stop
  use omp_lib, only: omp_get_wtime
  implicit none
  integer,parameter::P_TOTAL=1,P_PROBE=2,P_ALLOC=3,P_PREP=4,P_PACK=5,P_ADVANCE=6,P_TO=7,P_SAVE=8
  integer,parameter::P_FLAT=9,P_RECON=10,P_RIEMANN=11,P_UPDATE=12,P_FLUX=13,P_EOS=14,P_FROM=15
  integer,parameter::P_UNPACK=16,P_HOST_EOS=17,P_PUBLISH=18,P_RELEASE=19,NPHASE=19
  character(28),parameter::phaseNames(NPHASE)=[character(28):: &
    'GPU AMR adapter','GPU AMR device probe','GPU AMR allocation','GPU AMR CPU preparation', &
    'GPU AMR packing','GPU AMR batch advance','GPU AMR transfer to','GPU AMR reference copy', &
    'GPU AMR flattening','GPU AMR reconstruction','GPU AMR Riemann','GPU AMR update', &
    'GPU AMR flux accumulation','GPU AMR device EOS','GPU AMR transfer from', &
    'GPU AMR unpacking','GPU AMR host EOS','GPU AMR flux publication','GPU AMR release']
  type :: hy_amr_profile
    real(8)::seconds(NPHASE)=0
    integer::batches=0,blocks=0,minBatch=huge(0),maxBatch=0,capacityFlush=0,levelFlush=0,tailFlush=0
  end type
contains
  subroutine profile_start(phase,stamp)
    integer,intent(in)::phase
    real(8),intent(out)::stamp
    call Timers_start(trim(phaseNames(phase)))
    stamp=omp_get_wtime()
  end subroutine
  subroutine profile_stop(phase,stamp,profile)
    integer,intent(in)::phase
    real(8),intent(in)::stamp
    type(hy_amr_profile),intent(inout)::profile
    profile%seconds(phase)=profile%seconds(phase)+omp_get_wtime()-stamp
    call Timers_stop(trim(phaseNames(phase)))
  end subroutine
  subroutine profile_write(profile,capacity)
    use Driver_interface, only: Driver_getNStep,Driver_getMype,Driver_abort
    use RuntimeParameters_interface, only: RuntimeParameters_get
    type(hy_amr_profile),intent(in)::profile
    integer,intent(in)::capacity
    integer::interval,step,rank,unit,ios,n,smallest
    character(256)::prefix
    character(512)::filename
    real(8)::average,occupancy
    call RuntimeParameters_get('hy_gpuTimingInterval',interval)
    if(interval<=0) return
    call Driver_getNStep(step)
    if(mod(step,interval)/=0) return
    call Driver_getMype(MESH_COMM,rank)
    call RuntimeParameters_get('hy_gpuTimingPrefix',prefix)
    write(filename,'(A,"_step",I10.10,"_rank",I6.6,".csv")') trim(prefix),step,rank
    open(newunit=unit,file=trim(filename),status='replace',action='write',iostat=ios)
    if(ios/=0) call Driver_abort('HY_GPU_AMR: cannot open timing CSV')
    write(unit,'(A)',advance='no') 'step,rank,capacity,batches,blocks,min_batch,max_batch,' // &
      'capacity_flushes,level_flushes,tail_flushes,mean_batch,occupancy'
    do n=1,NPHASE
      write(unit,'(A)',advance='no') ','//trim(phaseNames(n))
    end do
    write(unit,*)
    smallest=0; average=0; occupancy=0
    if(profile%batches>0) then
      smallest=profile%minBatch
      average=real(profile%blocks,8)/profile%batches
      occupancy=average/capacity
    end if
    write(unit,'(I0,9(",",I0),2(",",ES26.17E3))',advance='no') step,rank,capacity,profile%batches, &
      profile%blocks,smallest,profile%maxBatch,profile%capacityFlush,profile%levelFlush,profile%tailFlush,average,occupancy
    do n=1,NPHASE
      write(unit,'(",",ES26.17E3)',advance='no') profile%seconds(n)
    end do
    write(unit,*)
    close(unit,iostat=ios)
    if(ios/=0) call Driver_abort('HY_GPU_AMR: timing CSV write failed')
  end subroutine
end module
