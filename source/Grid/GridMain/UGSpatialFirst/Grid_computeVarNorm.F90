!!****if* source/Grid/GridMain/UG/Grid_computeVarNorm
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
!!  Grid_computeVarNorm
!!
!! SYNOPSIS
!!  
!!  call Grid_computeVarNorm(integer(in)  :: level,
!!                           integer(in)  :: normType,
!!                           integer(in)  :: ivar,
!!                           real(out)    :: norm,
!!                           integer(in)  :: leafOnly)
!!
!! DESCRIPTION
!!
!!  Computes the L1 or L2 norm of the variable specified by ivar.  This
!!  can be done per-level, or on leaf or all nodes.  For multigrid, the
!!  L2 norm is used for convergence, but the L1 norm is incredibly useful
!!  for debugging purposes.
!!
!! ARGUMENTS
!!
!!  level     - If the norm is restricted to a given level; 0 is all
!!  normType - p in the Lp norm where choices of p are 1 or 2
!!  ivar      - the grid variable being normed; -1 for work
!!  norm      - the variable with which to return the norm
!!  leafOnly - if this isn't 0, compute the norm only on leaf nodes
!!
!! RESULT
!!
!!  The norm of ivar is in norm.
!!
!! NOTES
!!
!!  DEV: Currently only implemented for Paramesh4 and UG!
!!
!! EXAMPLE
!!  
!!  call gr_restrictTree()
!!  do i = 1, lrefine_max
!!    call Grid_computeVarNorm(i, 1, pdens, norm(i), 0)
!!  enddo
!!  do i = 1, lrefine_max
!!    if (norm(0) - norm(i) > 0.0000001) then
!!    call Driver_abort("restriction is highly nonconservatory!")
!!    endif
!!  enddo
!!
!!***

!!REORDER(5): unk

subroutine Grid_computeVarNorm (level, normType, ivar, norm, leafOnly)


#include "constants.h"

  use physicaldata, ONLY : unk
  use Grid_interface, ONLY : Grid_getCellVolumes, Grid_getTileIterator, Grid_releaseTileIterator
  use Driver_interface, ONLY : Driver_abort
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_data, ONLY : gr_meshComm
  use Grid_tile, ONLY : Grid_tile_t
  use Grid_iterator, ONLY : Grid_iterator_t

  implicit none

  include "Flashx_mpi.h"

  integer, intent(IN)  :: normType, level, ivar, leafOnly
  real, intent(OUT)    :: norm
  
  integer :: i, j, k, ierr
  real    :: lvol, lsum, bsum, sum
  real    :: cvol
  logical :: include_in_sum
  integer :: totalblockshere
  integer, dimension(LOW:HIGH,MDIM)   ::  blkLimits
  real, allocatable :: cellVolumes(:,:,:)
  type(Grid_tile_t) :: tileDesc
  type(Grid_iterator_t) :: itor
!===============================================================================

  call Timers_start("Grid_computeVarNorm")

  lvol = 0.
  lsum = 0.
  totalblockshere = 0

  if (normType /= 1 .and. normType /= 2) then
    call Driver_abort('only L1 and L2 norms supported in Grid_computeVarNorm!')
  endif

  include_in_sum = ((level == 1) .or. (level == 0))
  ! leafOnly must be ignored for UG
  if (include_in_sum) then
     call Grid_getTileIterator(itor, LEAF, level=1)
     do while (itor%isValid())
        call itor%currentTile(tileDesc)
        blkLimits(:,:) = tileDesc%limits
        allocate(cellVolumes(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))

        totalblockshere = totalblockshere + 1
        call Grid_getCellVolumes(tileDesc%level, &
                                 lbound(cellVolumes), ubound(cellVolumes), &
                                 cellVolumes)
        bsum = 0.
        if (ivar >= 0) then
           do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)  ! working on interior only
              do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
                 do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
                    cvol = cellVolumes(i,j,k)
                    bsum = bsum + abs(unk(ivar,i,j,k,1))**normType * cvol
                    lvol = lvol + cvol
                 enddo
              enddo
           enddo
!!        else
           ! DEV: Issue warning?
        endif
        lsum = lsum + bsum

        deallocate(cellVolumes)

        call itor%next()
     end do
     call Grid_releaseTileIterator(itor)
  endif

  call mpi_allreduce ( lsum, sum, 1, FLASH_REAL, & 
       MPI_SUM, gr_meshComm, ierr )
  !call mpi_allreduce ( lvol, vol, 1, FLASH_REAL,
  !                     MPI_SUM, FLASH_COMM, ierr )
  if (normType == 2) then
    norm = sqrt(sum)
  else if (normType == 1) then
    norm = sum
  endif

  call Timers_stop("Grid_computeVarNorm")

  !=================================================================
  
  return
end subroutine Grid_computeVarNorm
