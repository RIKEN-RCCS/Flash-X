!!****if* source/Grid/GridMain/UG/gr_createDomain
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
!!  gr_createDomain
!!
!!
!! SYNOPSIS
!!
!!  gr_createDomain()
!!                  
!!
!!
!! DESCRIPTION
!!
!! Creates the Uniform Grid domain. The Uniform Grid expects to be
!! given the processor grid at runtime in the form of runtime parameters
!! iprocs,jprocs and  kprocs, where iprocs*jprocs*kprocs is the
!! number of processor in the run. If UG is running in fixed
!! blocksize mode, then the blocksizes NXB, NYB and NZB are specified
!! at compile time. One block is placed on each processor, so that
!! the global domain size is <NXB*iprocs, NYB*jprocs, NZB*kprocs>.
!! However, if UG is running in non-fixed blocksize mode, the
!! blocksize is determined at runtime. In this mode, UG expects to be
!! given the global domain size in form of runtime parameters
!! iGridSize,jGridSize and kGridSize, and the blocksize is <iGridSize
!!/iprocs,jGridSize/jprocs,kGridSize/kprocs>. As in fixedblocksize mode,
!! only one blocks is placed on each processor
!!  
!! This routine also creates directional communicators used in
!! guardcell exchanges, allocates space for storing physical
!! coordinates of grid points, and computes them.
!!
!!  
!! ARGUMENTS
!!
!!
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine gr_createDomain()

  use physicalData, ONLY : unk, facevarx,facevary,facevarz

  use Grid_data, ONLY : scratch,scratch_ctr,&
       &scratch_facevarx,scratch_facevary,scratch_facevarz
  
  use Grid_data, ONLY : gr_axisMe, gr_axisNumProcs,&
       gr_guard,gr_maxCells, &
       gr_gIndexSize, gr_blkCornerID, &
       gr_lIndexSize, gr_kmax, gr_kmin, &
       gr_jmax, gr_jmin, gr_imax, gr_imin, gr_delta, gr_blkBC, &
       gr_iCoords, gr_jCoords, gr_kCoords, gr_domainBC, gr_iLoGC, &
       gr_jLoGC, gr_kLoGC, gr_iHiGC, gr_jHiGC, gr_kHiGC,&
       gr_ilo,gr_ihi,gr_jlo,gr_jhi,gr_klo,gr_khi,&
       gr_iguard,gr_jguard,gr_kguard, gr_flxx, gr_flxy, gr_flxz

  
  implicit none

#include "constants.h"
#include "Simulation.h"
  include "Flashx_mpi.h"

  integer :: ierr
  integer :: color, key, range_b, range_e
  real :: halfDelta
  integer :: i,j




  !store local index size for each block
  gr_lIndexSize = gr_gIndexSize/gr_axisNumProcs
  gr_maxCells = max(gr_lIndexSize(IAXIS),gr_lIndexSize(JAXIS),gr_lIndexSize(KAXIS))

  !store lower left global index for each dim
  gr_blkCornerID = gr_axisMe*gr_lIndexSize+1  

  
  gr_iguard = gr_guard(IAXIS)
  gr_ilo = gr_blkCornerID(IAXIS)
  gr_iloGc = gr_ilo-gr_iguard

  gr_ihi = gr_ilo+gr_lIndexSize(IAXIS) -1
  gr_ihiGc = gr_ihi + gr_iguard


  if(NDIM>1)then

     gr_jguard = gr_guard(JAXIS)
     gr_jlo = gr_blkCornerID(JAXIS)
     gr_jloGc = gr_jlo-gr_jguard
     
     gr_jhi = gr_jlo+gr_lIndexSize(JAXIS) -1
     gr_jhiGc = gr_jhi + gr_jguard

  else
     gr_jloGc = 1
     gr_jlo = 1
     gr_jhi = 1
     gr_jhiGc = 1
     gr_jguard = 0
  end if


  if(NDIM>2)then
     gr_kguard = gr_guard(KAXIS)
     gr_klo = gr_blkCornerID(KAXIS)
     gr_kloGc = gr_klo-gr_kguard
     
     gr_khi = gr_klo+gr_lIndexSize(KAXIS) -1
     gr_khiGc = gr_khi + gr_kguard

  else
     gr_klo=1
     gr_kloGc=1
     gr_khi = 1
     gr_khiGc = 1
     gr_kguard = 0
  endif
  
  !! Now create the grid and coordinates etc
  allocate(gr_iCoords(3,gr_iloGc:gr_ihiGc,1))
  allocate(gr_jCoords(3,gr_jloGc:gr_jhiGc,1))
  allocate(gr_kCoords(3,gr_kloGc:gr_khiGc,1))

#ifndef FIXEDBLOCKSIZE
  
  allocate(unk(gr_iLoGc:gr_iHiGc, gr_jLoGc:gr_jHiGc,&
       gr_kLoGc:gr_kHiGc,UNK_VARS_BEGIN:UNK_VARS_END,1))
  
  allocate(gr_flxx(gr_ilo:gr_ihi+1, gr_jlo:gr_jhi, &
       gr_klo:gr_khi,NFLUXES))
  allocate(gr_flxy(gr_ilo:gr_ihi, gr_jlo:gr_jhi+K2D, &
       gr_klo:gr_khi,NFLUXES))
  allocate(gr_flxz(gr_ilo:gr_ihi, gr_jlo:gr_jhi, &
       gr_klo:gr_khi+K3D,NFLUXES))
  
#if(NFACE_VARS>0)
  
  allocate(facevarx(gr_iLoGc:gr_iHiGc+1, gr_jLoGc:gr_jHiGc,&
       gr_kLoGc:gr_kHiGc,NFACE_VARS,1))
  
  allocate(facevary(gr_iLoGc:gr_iHiGc, gr_jLoGc:gr_jHiGc+K2D,&
       gr_kLoGc:gr_kHiGc,NFACE_VARS,1))
  
  allocate(facevarz(gr_iLoGc:gr_iHiGc, gr_jLoGc:gr_jHiGc,&
       gr_kLoGc:gr_kHiGc+K3D,NFACE_VARS,1))
  
#else
  allocate(facevarx(1,1,1,1,1))
  allocate(facevarz(1,1,1,1,1))
  allocate(facevary(1,1,1,1,1))
#endif
#if NSCRATCH_GRID_VARS > 0
  allocate(scratch(gr_iLoGc:gr_iHiGc+1, gr_jLoGc:gr_jHiGc+1,&
       gr_kLoGc:gr_kHiGc+1,SCRATCH_GRID_VARS_BEGIN:SCRATCH_GRID_VARS_END,1))
#else
  allocate(scratch(1,1,1,1,1))
#endif

#if NSCRATCH_CENTER_VARS > 0
  allocate(scratch_ctr(gr_iLoGc:gr_iHiGc, gr_jLoGc:gr_jHiGc,&
       gr_kLoGc:gr_kHiGc,SCRATCH_CENTER_VARS_BEGIN:SCRATCH_CENTER_VARS_END,1))
#else
  allocate(scratch_ctr(1,1,1,1,1))
#endif

#if(NSCRATCH_FACEX_VARS>0)  
  allocate(scratch_facevarx(gr_iLoGc:gr_iHiGc+1, gr_jLoGc:gr_jHiGc,&
       gr_kLoGc:gr_kHiGc,SCRATCH_FACEX_VARS_BEGIN:SCRATCH_FACEX_VARS_END,1))
#else
  allocate(scratch_facevarx(1,1,1,1,1))
#endif

#if(NSCRATCH_FACEY_VARS>0)  
  allocate(scratch_facevary(gr_iLoGc:gr_iHiGc, gr_jLoGc:gr_jHiGc+K2D,&
       gr_kLoGc:gr_kHiGc,SCRATCH_FACEY_VARS_BEGIN:SCRATCH_FACEY_VARS_END,1))
#else
  allocate(scratch_facevary(1,1,1,1,1))
#endif  

#if(NSCRATCH_FACEZ_VARS>0)
  allocate(scratch_facevarz(gr_iLoGc:gr_iHiGc, gr_jLoGc:gr_jHiGc,&
       gr_kLoGc:gr_kHiGc+K3D,SCRATCH_FACEZ_VARS_BEGIN:SCRATCH_FACEZ_VARS_END,1))
#else
  allocate(scratch_facevarz(1,1,1,1,1))
#endif
  
#endif

#ifdef DEBUG_GRID
  write(6,*)'gr_gIndexSize', gr_gIndexSize
  write(6,*)'gr_lIndexSize',gr_lIndexSize
  write(6,*)'gr_blkCornerID',gr_blkCornerID
#endif
  gr_delta(IAXIS,1) = (gr_imax-gr_imin)/gr_gIndexSize(IAXIS)
  
  halfDelta = gr_delta(IAXIS,1)/2.0
  j = gr_blkCornerID(IAXIS)-gr_guard(IAXIS)-1
  
  do i = gr_iloGc,gr_ihiGc
     gr_iCoords(LEFT_EDGE,i,1) = gr_imin+j*gr_delta(IAXIS,1)
     gr_iCoords(CENTER,i,1) = gr_imin+j*gr_delta(IAXIS,1)+halfDelta
     j = j+1
     gr_iCoords(RIGHT_EDGE,i,1) = gr_imin+j*gr_delta(IAXIS,1)
  end do
  gr_blkBC = gr_domainBC
  if(gr_axisMe(IAXIS)/=0)gr_blkBC(LOW,IAXIS)=NOT_BOUNDARY
  if(gr_axisMe(IAXIS)/=(gr_axisNumProcs(IAXIS)-1))&
       gr_blkBC(HIGH,IAXIS)=NOT_BOUNDARY
  
  if(NDIM > 1) then
     gr_delta(JAXIS,1) = (gr_jmax-gr_jmin)/gr_gIndexSize(JAXIS)
     halfDelta = gr_delta(JAXIS,1)/2.0
     j = gr_blkCornerID(JAXIS)-gr_guard(JAXIS)-1
     do i = gr_jloGc,gr_jhiGc
        gr_jCoords(LEFT_EDGE,i,1) = gr_jmin+j*gr_delta(JAXIS,1)
        gr_jCoords(CENTER,i,1) = gr_jmin+j*gr_delta(JAXIS,1)+halfDelta
        j = j+1
        gr_jCoords(RIGHT_EDGE,i,1) = gr_jmin+j*gr_delta(JAXIS,1)
     end do
     if(gr_axisMe(JAXIS)/=0)gr_blkBC(LOW,JAXIS)=NOT_BOUNDARY
     if(gr_axisMe(JAXIS)/=(gr_axisNumProcs(JAXIS)-1))&
          gr_blkBC(HIGH,JAXIS)=NOT_BOUNDARY
  else
     gr_jCoords(LEFT_EDGE:RIGHT_EDGE,1,1)=gr_jmin
     gr_delta(JAXIS,1)=0.0
  end if
  
  if(NDIM > 2) then
     gr_delta(KAXIS,1) = (gr_kmax-gr_kmin)/gr_gIndexSize(KAXIS)
     halfDelta = gr_delta(KAXIS,1)/2.0
     j = gr_blkCornerID(KAXIS)-gr_guard(KAXIS)-1
     do i = gr_kloGc,gr_khiGc
        gr_kCoords(LEFT_EDGE,i,1) = gr_kmin+j*gr_delta(KAXIS,1)
        gr_kCoords(CENTER,i,1) = gr_kmin+j*gr_delta(KAXIS,1)+halfDelta
        j = j+1
        gr_kCoords(RIGHT_EDGE,i,1) = gr_kmin+j*gr_delta(KAXIS,1)
     end do
     if(gr_axisMe(KAXIS)/=0)gr_blkBC(LOW,KAXIS)=NOT_BOUNDARY
     if(gr_axisMe(KAXIS)/=(gr_axisNumProcs(KAXIS)-1))&
          gr_blkBC(HIGH,KAXIS)=NOT_BOUNDARY
  else
     gr_kCoords(LEFT_EDGE:RIGHT_EDGE,1,1) = gr_kmin
     gr_delta(KAXIS,1)=0.0
  end if
  
  
end subroutine gr_createDomain
