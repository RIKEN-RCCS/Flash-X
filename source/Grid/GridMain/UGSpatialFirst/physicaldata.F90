!!****if* source/Grid/GridMain/UG/physicaldata
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
!!  physicaldata
!!
!! SYNOPSIS
!!
!!  use physicaldata
!!
!! DESCRIPTION
!!  
!!  This is the data modules that defines data structures for the 
!!  physical variables on the discretized mesh for Uniform Grid. 
!!  The cell centered data structure is "unk", the three face centered
!!  ones are "facevarx", "facevary" and "facevarz". There is also
!!  a data strucure to provide global domain scope scratch space called
!!  "scratch"
!!
!!***

module physicalData
#include "Simulation.h"

#ifdef FIXEDBLOCKSIZE
  real, save, target :: &
       unk(GRID_ILO_GC:GRID_IHI_GC, GRID_JLO_GC:GRID_JHI_GC,&
       GRID_KLO_GC:GRID_KHI_GC,UNK_VARS_BEGIN:UNK_VARS_END,1)

#if(NFACE_VARS>0)

  real, save, target :: &
       facevarx(GRID_ILO_GC:GRID_IHI_GC+1, GRID_JLO_GC:GRID_JHI_GC,&
       GRID_KLO_GC:GRID_KHI_GC,NFACE_VARS,1)
  real, save, target :: &
       facevary(GRID_ILO_GC:GRID_IHI_GC, GRID_JLO_GC:GRID_JHI_GC+K2D,&
       GRID_KLO_GC:GRID_KHI_GC,NFACE_VARS,1)
  real, save, target :: &
       facevarz(GRID_ILO_GC:GRID_IHI_GC, GRID_JLO_GC:GRID_JHI_GC,&
       GRID_KLO_GC:GRID_KHI_GC+K3D,NFACE_VARS,1)

#else
  real, save,target,dimension(1,1,1,1,1)::facevarx,facevary,facevarz
#endif

#else
  real, save, target, allocatable :: unk(:,:,:,:,:)
  real, save, target, allocatable :: facevarx(:,:,:,:,:)
  real, save, target, allocatable :: facevary(:,:,:,:,:)
  real, save, target, allocatable :: facevarz(:,:,:,:,:)
#endif

end module physicalData
