!!****if* source/Grid/GridMain/UG/gr_shiftData
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
!!  gr_shiftData
!!
!! SYNOPSIS
!!
!!  call gr_shiftData(integer (IN) :: comm,
!!                    integer (IN) :: dataType,
!!                    integer (IN) :: sendRight(4),
!!                    integer (IN) :: sendLeft(4),
!!                    integer (IN) :: recvRight(4),
!!                    integer (IN) :: recvLeft(4),
!!                    integer (IN) :: gridDataStruct)
!!  
!! DESCRIPTION 
!!  Given a communicator (comm) and derived mpi data type, "dataType", exchange 
!!  guardcell values with neighboring blocks by carrying out 
!!  out the shift data right and shift data left operations.
!!
!!  The arguments sendRight, sendLeft, recvRight, recvLeft hold the indices
!!  of the starting and receving points in the data structure specified by
!!  the gridDataStruct.  These values are closely tied to the dataType.
!!  
!!  
!! ARGUMENTS 
!!
!!  comm - the communicator. Here there will be one for each
!!          dimension of the process grid
!!  dataType - mpi datatype created in gr_createDataTypes
!!  sendRight - array holding the staring indices of the unk array for data sent to right neighbor
!!              (neighbor is up in 2nd dimension and back neighbor in 3rd dim)
!!  sendLeft - array holding the staring indices of the unk array for data sent to left neighbor
!!  recvRight - array holding the starting recv indices of the unk array for data sent to right neigh
!!  recvLeft - array holding the starting recv indices of the unk array for data sent to left neigh
!!  gridDataStruct - the datastructure in the grid, valid values are:
!!                   CENTER, FACEX, FACEY, FACEZ
!!
!!  EXAMPLE
!!   sendRight(1,1) in x direction holds the 1st unk dimension which is the variable ID
!!   sendRight(1,2) in x direction holds the 2nd unk dimension which is the i starting index
!!   sendRight(1,3) in x direction holds the 3rd unk dimension which is the j starting index
!!   sendRight(1,4) in x direction holds the 4th unk dimension which is the k starting index
!!
!!***
#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif



subroutine gr_shiftData(comm, dataType, sendRight, sendLeft, &
                        recvRight, recvLeft,gridDataStruct)

  use physicalData, ONLY : unk, facevarx, facevary, facevarz

#include "Flashx_mpi_implicitNone.fh"
#include "constants.h"
#include "Simulation.h"


  integer, intent(IN) :: comm
  integer, intent(IN) :: dataType  !The MPI_Datatype either a vector, 
                                   ! hvector or contig blk

  integer,dimension(MDIM+1),intent(IN) :: sendRight,sendLeft,recvRight,recvLeft
  integer, intent(IN) :: gridDataStruct
  
  integer :: localNumProcs, localMyPE, ierr
  integer :: dest, source
  integer,dimension(MPI_STATUS_SIZE) :: status


  !! For the current communicator, get the number of processors and
  !! the processor id of the PE.
  call MPI_Comm_Size(comm,localNumProcs,ierr)
  call MPI_Comm_rank(comm,localMyPE,ierr)
  
  !! First do the right shift, receive from the neighbor on the left and 
  !! send to the neighbor on the right

  
  if( localMyPE == 0) then
     source = localNumProcs-1
  else
     source = localMyPE-1
  end if
  
  if(localMyPE == localNumProcs-1) then
     dest = 0
  else
     dest = localMyPE+1
  end if
  
  
  !send guardcell data to neighbor on your right (up, or back)
  select case(gridDataStruct)
  case(CENTER)
     call MPI_SendRecv(unk(sendRight(1), sendRight(2), &
          sendRight(3), sendRight(4),1), &
          1, dataType, dest, 1,&
          unk(recvRight(1), recvRight(2), &
          recvRight(3), recvRight(4),1), &
          1, dataType,source, 1, comm, status,ierr)
  case(FACEX) 
     call MPI_SendRecv(facevarx(sendRight(1), sendRight(2), &
          sendRight(3), sendRight(4),1), &
          1, dataType, dest, 1,&
          facevarx(recvRight(1), recvRight(2), &
          recvRight(3), recvRight(4),1), &
          1, dataType,source, 1, comm, status,ierr)
     
  case(FACEY) 
     call MPI_SendRecv(facevary(sendRight(1), sendRight(2), &
          sendRight(3), sendRight(4),1), &
          1, dataType, dest, 1,&
          facevary(recvRight(1), recvRight(2), &
          recvRight(3), recvRight(4),1), &
          1, dataType,source, 1, comm, status,ierr)
     
  case(FACEZ) 
     call MPI_SendRecv(facevarz(sendRight(1), sendRight(2), &
          sendRight(3), sendRight(4),1), &
          1, dataType, dest, 1,&
          facevarz(recvRight(1), recvRight(2), &
          recvRight(3), recvRight(4),1), &
          1, dataType,source, 1, comm, status,ierr)
     
  end select
  !! Now do the left shift, receive from the neighbor on the right and 
  !! send to the neighbor on the left
  
  if( localMyPE == 0) then
     dest = localNumProcs-1
  else
     dest = localMyPE-1
  end if
  
  if(localMyPE == localNumProcs-1) then
     source = 0
  else
     source = localMyPE+1
  end if

  !  write(*,*)sendLeft
  !  write(*,*)recvLeft
  
  !send guardcell data to neighbor on your left (down or front)
  select case(gridDataStruct)
  case(CENTER)
     call MPI_SendRecv(unk(sendLeft(1), sendLeft(2), &
          sendLeft(3), sendLeft(4),1), &
          1, dataType, dest, 1,&
          unk(recvLeft(1), recvLeft(2), recvLeft(3), &
          recvLeft(4),1), &
          1, dataType, source, 1, comm, status,ierr)
  case(FACEX)
     call MPI_SendRecv(facevarx(sendLeft(1), sendLeft(2), &
          sendLeft(3), sendLeft(4),1), &
          1, dataType, dest, 1,&
          facevarx(recvLeft(1), recvLeft(2), recvLeft(3), &
          recvLeft(4),1), &
          1, dataType, source, 1, comm, status,ierr)
  case(FACEY)
     call MPI_SendRecv(facevary(sendLeft(1), sendLeft(2), &
          sendLeft(3), sendLeft(4),1), &
          1, dataType, dest, 1,&
          facevary(recvLeft(1), recvLeft(2), recvLeft(3), &
          recvLeft(4),1), &
          1, dataType, source, 1, comm, status,ierr)
  case(FACEZ)
     call MPI_SendRecv(facevarz(sendLeft(1), sendLeft(2), &
          sendLeft(3), sendLeft(4),1), &
          1, dataType, dest, 1,&
          facevarz(recvLeft(1), recvLeft(2), recvLeft(3), &
          recvLeft(4),1), &
          1, dataType, source, 1, comm, status,ierr)
  end select
  
end subroutine gr_shiftData


