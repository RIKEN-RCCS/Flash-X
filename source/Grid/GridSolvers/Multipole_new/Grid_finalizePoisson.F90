!!****if* source/Grid/GridSolvers/Multipole_new/Grid_solvePoisson
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
!!  Grid_finalizePoisson
!!
!! SYNOPSIS
!!
!!  Grid_finalizePoisson (integer(in)    :: iSoln,
!!                     integer(in)    :: iSrc, 
!!                     integer(6)(in) :: bcTypes,
!!                     real(2,6)(in)  :: bcValues,
!!                     real(inout)    :: poisFact)
!!
!! DESCRIPTION
!!
!!  Poisson solver routine.  This module implements the multipole
!!  summation method for isolated boundary problems.  Periodic
!!  problems are not supported; the multipole method should be
!!  used only for approximately spherical matter distributions.
!!
!! ARGUMENTS
!!
!!  iSoln          : index to variable containing potential
!!  iSrc           : index to variable containing density
!!  bcTypes        : boundary types along various faces,
!!                   only used in verifying that they are isolated
!!  bcValues       : the values to boundary conditions, currently not used
!!  poisFact       : factor to be used in calculation
!!
!!***

subroutine Grid_finalizePoisson (iSoln,                   &
                              iSrc,                    &
                              bcTypes,                 &
                              bcValues,                &
                              poisFact)

  use Grid_interface,    ONLY : GRID_PDE_BND_ISOLATED

  use Driver_interface,  ONLY : Driver_abort

  use Timers_interface,  ONLY : Timers_start,                  &
                                Timers_stop

  use gr_mpoleInterface, ONLY : gr_mpolePotentials,            &
                                gr_mpoleDumpMoments,           &
                                gr_mpoleDeallocateRadialArrays

  use gr_mpoleData,      ONLY : gr_mpoleMomentsDump,           &
                                gr_mpoleMultiThreading,        &
                                gr_mpoleRadialInfoPrint,       &
                                gr_mpoleMomentR,               &
                                gr_mpoleMomentI,               &
                                gr_mpoleScratch,               &
                                gr_mpoleMaxQ,                  &
                                gr_mpoleGravityConstant,       &
                                gr_mpoleFourPiInv,             &
                                gr_mpoleGeometry,              &
                                gr_mpoleRequest

  implicit none

#include "Simulation.h"
#include "constants.h"

  integer, intent(in)    :: iSoln, iSrc
  integer, intent(in)    :: bcTypes (6)
  real,    intent(in)    :: bcValues (2,6)
  real,    intent(inout) :: poisFact
!  
  call Timers_start ("Multipole Solver Finalize")
!  
  
  call Timers_start             ("gr_mpolePotentials")
  call gr_mpolePotentials       (iSoln, poisFact,               &
                                   gr_mpoleMomentR,              &
                                   gr_mpoleMomentI,              &
                                   gr_mpoleScratch,              &
                                   gr_mpoleMaxQ,                 &
                                   gr_mpoleGravityConstant,      &
                                   gr_mpoleFourPiInv,            &
                                   gr_mpoleGeometry,             &
                                   gr_mpoleRequest)
  call Timers_stop              ("gr_mpolePotentials")
  
  !
  !
  !    ...Dump the moments if requested by the user.
  !
  !
  if (gr_mpoleMomentsDump) then
     call gr_mpoleDumpMoments ()
  end if
  !
  !
  !    ...Final chores.
  !
  !
  call gr_mpoleDeallocateRadialArrays ()
  !  
  !
  !    ...End timer.
  !
  !
  call Timers_stop ("Multipole Solver Finalize")
!
!
!    ...Check the potentials obtained. This routine is only internal
!       and should only be activated by the author of the code for
!       debugging and accuracy check purposes.
!
!
!  call gr_mpolePotential_exact   (iSrc,iSoln, poisFact)
!  
!
!    ...Ready!
!
!
  return
end subroutine Grid_finalizePoisson
