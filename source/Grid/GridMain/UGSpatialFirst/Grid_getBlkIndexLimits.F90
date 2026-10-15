!!****if* source/Grid/GridMain/UG/Grid_getBlkIndexLimits
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
!!  Grid_getBlkIndexLimits
!!
!! SYNOPSIS
!!
!!  call Grid_getBlkIndexLimits(integer(IN)  :: blockId,
!!                              integer(OUT) :: blkLimits(2,MDIM),
!!                              integer(OUT) :: blkLimitsGC(2,MDIM),
!!                     optional,integer(IN)  :: gridDataStruct)
!!  
!! DESCRIPTION 
!!  
!! Returns the indices of the upper and lower bounds of a block.
!! blkLimitsGC holds the entire dimension of the block (including guardcells)
!! blkLimits returns the indices of the interior of the block. 
!!
!! The first dimension holds the lower and upper bounds of the block.
!! The second dimension is of size MDIM to hold the upper and lower 
!! indices of each dimension.  
!!
!! This routine is also used to calculate the size of the internal 
!! index size of a block.  In FLASH2, this was simply NXB, NYB and NZB.
!! In FLASH3, since we are supporting blocksizes that are not fixed you 
!! should calculate on the size of the block on a per block basis.  See
!! example below.
!!  
!! Please note that the face centered variables have an extra data point 
!! along the corresponding dimension. For example variable centered on the
!! faces along IAXIS dimension have equivalent of NXB+1 data points in the
!! interior. If the optional argument gridDataStruct is specified, and it
!! is one of the face-centered data structures, the values returned in the
!! the two output arrays will reflect the additional data point. Since most
!! of FLASH applications use only cell centered variable, most users need not
!! concern themselves with this argument.
!!
!!
!! ARGUMENTS 
!!
!!  blockId - the local blockID
!!  blkLimits  -  array returned holding the upper and lower indices of
!!                the interior of a block (no guardcells!)
!!             
!!  blkLimitsGC  -  array returned holding the upper and lower indices
!!                  of an entire block, that is, including guardcells.
!!  gridDataStruct - Grid data structure being operated on, the valid
!!                   values are CENTER, CENTER_FACES, FACEX, FACEY, and FACEZ.
!!                   Specifying gridDataStructure= CENTER or CENTER_FACES is
!!                   equivalent to not including the optional argument in the
!!                   call.
!!             
!! EXAMPLE
!! 
!!  Take a 2 dimensional block of size NXB x NYB, with ghost cells along
!!  each dimension being GX and GY.  This block could be stored in an array of 
!!  size (NXB+2*GX,NYB+2*GY). For this array the returned values of blkLimits and blkLimitsGC
!!  are as follows.
!!
!!  blkLimitsGC(LOW, IAXIS) = 1 !lower bound index of the first cell in the entire block in xdir
!!  blkLimitsGC(LOW, JAXIS) = 1
!!  blkLimitsGC(LOW, KAXIS) = 1
!!
!!  blkLimitsGC(HIGH, IAXIS) = NXB+2*GX !upper bound index of the last cell in the entire block in xdir
!!  blkLimitsGC(HIGH, JAXIS) = NYB +2*GY
!!  blkLimitsGC(HIGH, KAXIS) = 1 !because there are only 2 dimensions
!!
!!
!!  blkLimits(LOW, IAXIS) = GX+1 !lower bound index of the first interior cell of a block in xdir
!!  blkLimits(LOW, JAXIS) = GY+1 
!!  blkLimits(LOW, KAXIS) = 1
!!
!!  blkLimits(HIGH, IAXIS) = NXB+GX !the upper bound index of the last interior cell of a block in xdir
!!  blkLimits(HIGH, JAXIS) = NYB+GY
!!  blkLimits(HIGH, KAXIS) = 1 !because there are only 2 dimensions
!!
!!  For the same block, if the Optional argument was present and was FACEX, then
!!  blkLimits(HIGH,IAXIS)=NXB+1, and blkLimitsGC(HIGH,IAXIS)=NXB+1+2*GX, and 
!!  all other values of blkLimits and blkLimitsGC are the same as those of the
!!  default case. Similarly if optional argument had a value FACEY, then
!!  blkLimits(HIGH,JAXIS)=NYB+1, and blkLimitsGC(HIGH,JAXIS)=NYB+1+2*GY, and
!!  all others are same as before. 
!!
!!
!!  To calculate NXB, NYB and NZB for a given block
!!
!!  NXB = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) +1
!!  NYB = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) +1
!!  NZB = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) +1
!!
!!
!! NOTES
!!  
!!   The #define constants IAXIS, JAXIS and KAXIS
!!   are defined in constants.h and are
!!   meant to ease the readability of the code.       
!!   instead of blkLimits(HIGH,3)  the code reads
!!   blkLimits(HIGH, KAXIS)
!!
!!
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

module Grid_getBlkIndexLimits_mod
contains
subroutine Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC, gridDataStruct)
  use Grid_data, ONLY : gr_ilo,gr_ihi,gr_jlo,gr_jhi,gr_klo,gr_khi,&
                        gr_iloGc,gr_ihiGc,gr_jloGc,gr_jhiGc,gr_kloGc,gr_khiGc,&
                        gr_meshMe
  use Driver_interface, ONLY : Driver_abort
  implicit none

#include "constants.h"
#include "Simulation.h"
  integer,intent(IN) :: blockID
  integer,dimension(LOW:HIGH,MDIM),intent(OUT) :: blkLimits,blkLimitsGC
  integer,optional, intent(IN) :: gridDataStruct
  integer,dimension(MDIM) :: faces

  blkLimits(LOW,IAXIS) = gr_ilo
  blkLimits(HIGH,IAXIS) = gr_ihi
  blkLimits(LOW,JAXIS) = gr_jlo
  blkLimits(HIGH,JAXIS) = gr_jhi
  blkLimits(LOW,KAXIS) = gr_klo
  blkLimits(HIGH,KAXIS) = gr_khi

  blkLimitsGC(LOW,IAXIS) = gr_iloGc
  blkLimitsGC(HIGH,IAXIS) = gr_ihiGc
  blkLimitsGC(LOW,JAXIS) = gr_jloGc
  blkLimitsGC(HIGH,JAXIS) = gr_jhiGc
  blkLimitsGC(LOW,KAXIS) = gr_kloGc
  blkLimitsGC(HIGH,KAXIS) = gr_khiGc
  if(present(gridDataStruct)) then
     faces = 0
     if((gridDataStruct==FACEX).or.(gridDataStruct==SCRATCH_FACEX)) then
        faces(IAXIS)=1
     elseif((gridDataStruct==FACEY).or.(gridDataStruct==SCRATCH_FACEY)) then
        faces(JAXIS)=1
     elseif((gridDataStruct==FACEZ).or.(gridDataStruct==SCRATCH_FACEZ)) then
        faces(KAXIS)=1
     elseif(gridDataStruct==SCRATCH) then
        !! This will have to be removed once we eliminate SCRATCH data struct
        faces = 1
     elseif((gridDataStruct/=CENTER).and.(gridDataStruct/=SCRATCH_CTR) .and. &
            (gridDataStruct/=CENTER_FACES)) then
        if(gr_meshMe == MASTER_PE)print*,'In BlkIndexLimits the provided data structure is',gridDataStruct
        call Driver_abort("called index limits with invalid gridDataStruct")
     end if
     blkLimitsGC(HIGH,:)=blkLimitsGC(HIGH,:)+faces
     blkLimits(HIGH,:)=blkLimits(HIGH,:)+faces
  end if
  
end subroutine Grid_getBlkIndexLimits

end module
