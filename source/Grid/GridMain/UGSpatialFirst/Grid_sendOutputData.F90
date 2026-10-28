!!****if* source/Grid/GridMain/UG/Grid_sendOutputData
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
!!  Grid_sendOutputData
!!
!! SYNOPSIS
!!
!!  call Grid_sendOutputData()
!!  
!! DESCRIPTION 
!!   This routine prepares the Grid_IO data that is needed 
!!   in order to write the data to a checkpoint or plotfile.
!!   Data includes, ngid (nfaces + nchildren + 1 parent), globalOffset
!!   globalNumBlocks and other data that is specific to the grid type, 
!!   Paramesh, UG or other
!!
!!  ARGUMENTS  
!!
!!***

subroutine Grid_sendOutputData()

  use Grid_data, ONLY : gr_axisNumProcs, gr_str_geometry, &
       gr_domainBC,                                       &
       gr_globalNumBlocks, gr_gid,  gr_meshMe, gr_meshNumProcs
  use gr_specificData, ONLY : gr_globalOffset
  use IO_interface, ONLY : IO_setScalar

  implicit none
#include "constants.h"
#include "Simulation.h"


  integer :: leftNeigh, rightNeigh, blockID

  !Put NXB, NYB and NZB into a saved variable to prevent problems with
  !an xlf "feature."  These don't mean a whole lot in nofbs mode.
  integer,parameter :: local_nxb = NXB
  integer,parameter :: local_nyb = NYB
  integer,parameter :: local_nzb = NZB
  integer,parameter :: dimensionality = NDIM

  gr_globalNumBlocks = gr_meshNumProcs

  gr_globalOffset = mod(gr_meshMe, gr_globalNumBlocks)



  !set the scalars for the grid unit
  call IO_setScalar("nxb", local_nxb)
  call IO_setScalar("nyb", local_nyb)
  call IO_setScalar("nzb", local_nzb)
  call IO_setScalar("dimensionality", dimensionality)

  
  call IO_setScalar("globalNumBlocks", gr_globalNumBlocks)
  
  call IO_setScalar("geometry", gr_str_geometry)

  
  !!Write the global ID
  
  !-------------------------------------------------------------------------
  ! compute the global id -- this is a single array which stores the 
  ! neighbor block numbers, the parent, and the children of a given block
  !-------------------------------------------------------------------------
  

  !! get the neighbor blocks - the flash3 UG way is to just to number
  !! them in fortran ordering , gid(ngid, blockID)
  !! in UG there is always one block per proc so blockID = 1
  
  blockID = 1
  
  
  if (NDIM >= 1) then
     leftNeigh = gr_meshMe - 1 
     rightNeigh = gr_meshMe + 1
     !! gr_axisNumProcs is the number of procs in the x dim     
     if (mod(gr_meshMe, gr_axisNumProcs(1)) == 0) then
        leftNeigh = gr_domainBC(LOW,IAXIS)
     endif
     
     if (mod(gr_meshMe + 1, gr_axisNumProcs(1)) == 0) then
        rightNeigh = gr_domainBC(HIGH,IAXIS)
     endif


     gr_gid(1, blockID) = leftNeigh + max(0,sign(1,leftNeigh))
     gr_gid(2, blockID) = rightNeigh + max(0,sign(1,rightNeigh))
     !! non existent for parent and children
     gr_gid(3, blockID) = -1
     gr_gid(4, blockID) = -1
     gr_gid(5, blockID) = -1
     
  end if
  
  
#if NDIM >= 2

  leftNeigh = gr_meshMe - gr_axisNumProcs(1)
  rightNeigh = gr_meshMe + gr_axisNumProcs(1)

  if (mod(gr_meshMe/gr_axisNumProcs(1), gr_axisNumProcs(2)) == 0) then
     leftNeigh = gr_domainBC(LOW,JAXIS)
  end if

  if (mod(rightNeigh/gr_axisNumProcs(1), gr_axisNumProcs(2)) == 0) then
     rightNeigh = gr_domainBC(HIGH,JAXIS)
  end if

  gr_gid(3, blockID) = leftNeigh + max(0,sign(1,leftNeigh))
  gr_gid(4, blockID) = rightNeigh + max(0,sign(1,rightNeigh))

  !! set parents and children to non existent
  gr_gid(5, blockID) = -1
  gr_gid(6, blockID) = -1
  gr_gid(7, blockID) = -1
  gr_gid(8, blockID) = -1
  gr_gid(9, blockID) = -1

#endif

#if NDIM > 2

  leftNeigh = gr_meshMe - (gr_axisNumProcs(1) * gr_axisNumProcs(2))
  rightNeigh = gr_meshMe + (gr_axisNumProcs(1) * gr_axisNumProcs(2))

  if(leftNeigh < 0) then
     leftNeigh = gr_domainBC(LOW,KAXIS)
  end if
  
  if(rightNeigh >= gr_meshNumProcs) then
     rightNeigh = gr_domainBC(HIGH,KAXIS)
  end if

  gr_gid(5, blockID) = leftNeigh + max(0,sign(1,leftNeigh))
  gr_gid(6, blockID) = rightNeigh + max(0,sign(1,rightNeigh))

  !! set parents and children to non existent
  gr_gid(7, blockID) = -1
  gr_gid(8, blockID) = -1
  gr_gid(9, blockID) = -1
  gr_gid(10, blockID) = -1
  gr_gid(11, blockID) = -1
  gr_gid(12, blockID) = -1
  gr_gid(13, blockID) = -1
  gr_gid(14, blockID) = -1
  gr_gid(15, blockID) = -1

#endif

  





end subroutine Grid_sendOutputData
