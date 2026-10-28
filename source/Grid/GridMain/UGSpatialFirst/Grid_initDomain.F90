!!****if* source/Grid/GridMain/UG/Grid_initDomain
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
!!  Grid_initDomain
!!
!!
!! SYNOPSIS
!!
!!  call Grid_initDomain(logical(IN)  :: restart,
!!                       logical(INOUT) :: particlesInitialized)
!!
!!
!! DESCRIPTION
!! 
!!  Create the mesh, initialize all the mesh data structures
!!  and apply initial conditions
!!
!!
!! ARGUMENTS
!!
!!  restart : is true if the execution is starting from a checkpoint
!!            file, otherwise false.
!!  particlesInitialized : is true if particle positions were initialized before returning
!!                         from this routine
!!
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_initDomain( restart,particlesInitialized)

  use Grid_interface, ONLY : Grid_getTileIterator, Grid_releaseTileIterator, &
       Grid_fillGuardCells
  
  use Grid_data, ONLY : gr_gid, gr_eosModeInit, gr_eosMode, gr_eosModeNow, gr_globalMe,&
       gr_globalNumProcs
  use Eos_interface, ONLY : Eos_multiDim
  use Simulation_interface, ONLY : Simulation_initBlock, Simulation_initRestart
  use Grid_tile, ONLY : Grid_tile_t
  use Grid_iterator, ONLY : Grid_iterator_t
  implicit none

#include "Simulation.h"
#include "constants.h"
  type(Grid_tile_t) :: tileDesc
  type(Grid_iterator_t) :: itor
  
  logical, intent(in) :: restart
  logical, intent(inout) :: particlesInitialized

  integer :: blockID=1, ngid
  
  integer, dimension(LOW:HIGH, MDIM) :: blkLimitsGC
  real, pointer:: solnData(:,:,:,:)

  nullify(solnData)

  call gr_createDataTypes()

  if (NDIM == 1) then
        ngid = 5
     else if (NDIM == 2) then
        ngid = 9
     else if (NDIM == 3) then
        ngid = 15
  end if

  allocate(gr_gid(ngid,1)) !1 because in UG only 1 block per proc     

  if(.not.restart) then
     !  zero data in case we don't initialize all data in
     !  Simulation_initBlock... in particular the total vs. internal
     !  energies can cause problems in the eos call that follows
     call Grid_getTileIterator(itor, ALL_BLKS, tiling=.FALSE.)
     do while (itor%isValid())
        call itor%currentTile(tileDesc)
        call tileDesc%getDataPtr(solnData, CENTER)

        solnData(:,:,:,:) = 0.0
        call Simulation_initBlock(solnData, tileDesc)
        call Eos_multiDim(gr_eosModeInit, tileDesc%limits, tileDesc%blkLimitsGC(LOW,:), solnData)

        call tileDesc%releaseDataPtr(solnData, CENTER)
        call itor%next()
     end do
     call Grid_releaseTileIterator(itor)
  else
     ! Do no call the EOS here any more when restarting, since those calls
     ! may (depending on the Eos implementation, or whether Eos has
     ! been called explicitly before the checkpoint was written)
     ! introduce small data differences.
     ! As long as checkpoint data is written in a state where the
     ! solution data is thermodynamically consistent, there is no
     ! need for the call here. - KW

     ! Now give user code an opportunity to modify any values on restart.
     ! User implementation of Simulation_initRestart is responsible for
     ! calling Eos if that is needed to leave the solution data in a
     ! consistent state.

     call Simulation_initRestart()

  end if

  call Grid_fillGuardCells( CENTER_FACES, ALLDIR)

  gr_eosModeNow = gr_eosMode
  call gr_solversInit()


end subroutine Grid_initDomain
