!!****if* source/Grid/GridMain/UG/Grid_getSingleCellCoords
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
!!  Grid_getSingleCellCoords
!!
!! SYNOPSIS
!!
!!  Grid_getSingleCellCoords(integer(in)       ::  ind(MDIM), 
!!                   type(block_metadta_t)(IN) :: block,
!!                           integer(in)       ::  edge,
!!                           integer(in)       ::  beginCount,
!!                           real(OUT)         ::  coords(MDIM))
!!  
!! DESCRIPTION 
!!
!!  Returns the coordinates of a single cell of a given block.
!!
!!
!! ARGUMENTS
!!
!!  ind - array holding the indices of the cell whose coordinates to return
!! 
!!   block - block metadata
!!
!!  edge - indicates if user wants the left, center or right edge of the block
!!         options are LEFT_EDGE, RIGHT_EDGE, or CENTER.  These constants are
!!         defined in constants.h
!!
!!
!!  beginCount : tells the routine where to start index counting.  beginCount can
!!               be set to INTERIOR or EXTERIOR.  If INTERIOR is specified
!!               guardcell indices are not included and index 1 is the first interior cell. 
!!               If EXTERIOR is specified
!!               the first index, 1, is the left most guardcell.  See example
!!               below for more explanation.  (For most of the FLASH architecture code,
!!               we use EXTERIOR.  Some physics routines, however, find it helpful 
!!               only to work on the internal parts of the blocks (without
!!               guardcells) and wish to keep loop indicies  
!!               going from 1 to NXB without having to worry about finding 
!!               the correct offset for the number of guardcells.) 
!!               (INTERIOR and EXTERIOR are defined in constants.h)
!!
!!
!!    
!!  coords : returned coordinates of the specificed cell
!!
!!
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "constants.h"
#include "Simulation.h"

subroutine Grid_getSingleCellCoords(ind, level, edge, coords)
  use Driver_interface, ONLY : Driver_abort
  use Grid_data,        ONLY : gr_iCoords, gr_jCoords, gr_kCoords, &
                               gr_guard

  implicit none

  integer, intent(in)  :: ind(1:MDIM)
  integer, intent(in)  :: level
  integer, intent(in)  :: edge
  real,    intent(out) :: coords(1:MDIM)

  coords(:) = 0.0
  call Driver_abort("[Grid_getSingleCellCoords] DEPRECATED")

!#ifdef DEBUG_GRID
!  print*,' inside Grid_getSingleCellCoords', ind, level, edge, coords
!  if((level /= 1)) then
!     print*,"Grid_getSingleCellCoords: invalid level"
!     call Driver_abort("[Grid_getSingleCellCoords] invalid level")
!  end if
!  if((edge/=LEFT_EDGE).and.(edge/=CENTER).and.(edge/=RIGHT_EDGE))&
!       call Driver_abort('Grid_getSingleCellCoords : invalid edge')
!
!  print*, 'leaving the DEBUG_GRID statement'
!#endif
!
!  coords(IAXIS) = gr_iCoords(edge,ind(IAXIS)+gr_guard(IAXIS),1)
!  coords(JAXIS) = gr_jCoords(edge,ind(JAXIS)+gr_guard(JAXIS),1)
!  coords(KAXIS) = gr_kCoords(edge,ind(KAXIS)+gr_guard(KAXIS),1)
end subroutine Grid_getSingleCellCoords

