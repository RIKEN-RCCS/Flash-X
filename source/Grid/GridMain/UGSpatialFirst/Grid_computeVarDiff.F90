!!****if* source/Grid/GridMain/UG/Grid_computeVarDiff
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
!!  Grid_computeVarDiff
!!
!! SYNOPSIS
!!  
!!  call Grid_computeVarDiff(integer(in) :: level,
!!                integer(in) :: gr_iRefSoln,
!!                integer(in) :: gr_iSoln,
!!                integer(in) :: ires)
!!
!! DESCRIPTION
!!  
!!  Compute the current difference of two variables.
!!  If level is -1, only the leaf block residuals are calculated;
!!  otherwise, residuals are computed for all blocks of the given level.
!!
!!
!! ARGUMENTS
!!
!!  level        - the level of blocks to take the residual for, or -1
!!                 for all LEAF blocks.
!!  gr_iRefSoln  - the "reference solution" variable
!!  gr_iSoln     - "the solution" variable
!!  ires         - the residual variable
!!
!! RESULT
!!
!!  The difference between the reference solution and solution in the blocks
!!  (either selected by level, or all LEAF blocks)
!!  is placed in the variable ires.
!!
!!***

!!REORDER(5): unk

subroutine Grid_computeVarDiff(level, gr_iRefSoln, gr_iSoln, ires)

!==============================================================================
#include "Simulation.h"
#include "constants.h"


  use physicaldata, ONLY : unk

  use Timers_interface, ONLY : Timers_start, Timers_stop
  
  
  implicit none

  integer, intent(in)          :: level, gr_iRefSoln, gr_iSoln, ires
  
  integer, dimension(2,MDIM)   :: blkLimitsGC, blkLimits
  integer                      :: i, j, k

  
  !=======================================================================

  call Timers_start("Grid_computeVarDiff")
  
  
  
  if (level==-1 .OR. level == 1) then
     unk(ires,:,:,:,1) = unk(gr_iSoln, :,:,:,1) - unk(gr_iRefSoln, :,:,:,1)
  endif


  call Timers_stop("Grid_computeVarDiff")

  !=======================================================================
  
  return
end subroutine Grid_computeVarDiff
