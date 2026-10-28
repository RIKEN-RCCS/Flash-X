!!****if* source/Grid/GridMain/UG/gr_createDataTypes
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
!!  gr_createDataTypes
!!
!! SYNOPSIS
!!
!!  call gr_createDataTypes()
!!  
!! DESCRIPTION 
!!  create the MPI derived Datatypes in the uniform grid for exchanging
!!  guardcells.  It is more efficient to make a datatype than to do multiple
!!  mpi send/recvs
!!  
!!  
!! ARGUMENTS 
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine gr_createDataTypes()

  use Grid_data, ONLY : gr_exch,&
       gr_numDataStruct,gr_gridDataStruct, gr_gridDataStructSize
  use Grid_interface, ONLY : Grid_getBlkIndexLimits

#include "Flashx_mpi_implicitNone.fh"
#include "constants.h"
#include "Simulation.h"


  integer :: ierr, i
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer, dimension(MDIM) :: guard, blkExtent

  do i=1,gr_numDataStruct
     call Grid_getBlkIndexLimits(1,blkLimits,blkLimitsGC,gr_gridDataStruct(i))
     guard(:) = blkLimits(LOW,:)-blkLimitsGC(LOW,:)
     blkExtent(:)=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1

     call MPI_TYPE_VECTOR(blkExtent(JAXIS)*blkExtent(KAXIS)*&
          gr_gridDataStructSize(i), guard(IAXIS), blkExtent(IAXIS), &
          FLASH_REAL, gr_exch(i,IAXIS), ierr)
     call MPI_TYPE_COMMIT(gr_exch(i,IAXIS),ierr)

     if(NDIM>1) then
        call MPI_TYPE_VECTOR(blkExtent(KAXIS)*gr_gridDataStructSize(i), &
             blkExtent(IAXIS)*guard(JAXIS), &
             blkExtent(IAXIS)*blkExtent(JAXIS), FLASH_REAL, &
             gr_exch(i,JAXIS), ierr)
        call MPI_TYPE_COMMIT(gr_exch(i,JAXIS),ierr)
     end if

     if(NDIM>2) then
        call MPI_TYPE_VECTOR(gr_gridDataStructSize(i), &
             blkExtent(IAXIS)*blkExtent(JAXIS)*guard(KAXIS), &
             blkExtent(IAXIS)*blkExtent(JAXIS)*blkExtent(KAXIS), &
             FLASH_REAL, gr_exch(i,KAXIS), ierr)
        call MPI_TYPE_COMMIT(gr_exch(i,KAXIS),ierr)
     end if
  end do
end subroutine gr_createDataTypes
   
   
   
