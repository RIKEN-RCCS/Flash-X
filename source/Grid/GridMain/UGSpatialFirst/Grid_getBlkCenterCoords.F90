!!****if* source/Grid/GridMain/UG/Grid_getBlkCenterCoords
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
!!  Grid_getBlkCenterCoords
!!
!! SYNOPSIS
!!
!!  Grid_getBlkCenterCoords(type(block_metadta_t)(IN) :: block,
!!                         real(OUT)   :: blockCenter(MDIM))
!!                      
!!  
!!  
!! DESCRIPTION 
!!  Gets the coordinates of the center of the block identified by
!!  blockId.  Returns the coordinates in an array blockCenter  
!!
!!
!! ARGUMENTS
!!  blockId - local id number of the block. for UG always 1
!!  blockCenter - returned array of size MDIM holding the blockCenter coords  
!! 
!! Example
!!   In 2 dimensions, if physical coordinates are ...
!!    
!!     ________________(0.5 1.0)
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |_______________ |
!!  (-0.5, 0.0)
!!
!!  then the values returned in blockCenter are 
!!  blockCenter(IAXIS) = 0.0
!!  blockCenter(JAXIS) = 0.5
!!  blockCenter(KAXIS) = 0.0 since the dimension is not included  
!!
!!
!!***

#include "constants.h"
#include "Simulation.h"

subroutine Grid_getBlkCenterCoords(block,blockCenter)
  use Driver_interface, ONLY : Driver_abort
!  use Grid_interface, ONLY : Grid_getBlkBoundBox
  use Grid_tile, ONLY : Grid_tile_t

  implicit none

  type(Grid_tile_t), intent(in) :: block
  real,dimension(MDIM),intent(out) :: blockCenter

  blockCenter(:) = 0.0
  call Driver_abort("[Grid_getBlkCenterCoords] Implement for tiling")

!  real,dimension(2,MDIM) :: bndBox
!  call Grid_getBlkBoundBox(block, bndBox)
!
!  blockCenter(1) = bndBox(1,1) + (bndBox(2,1) - bndBox(1,1))/2
!  if(NDIM>1) blockCenter(2) = bndBox(1,2) + (bndBox(2,2) - bndBox(1,2))/2
!  if(NDIM>2) blockCenter(3) = bndBox(1,3) + (bndBox(2,3) - bndBox(1,3))/2
end subroutine Grid_getBlkCenterCoords

