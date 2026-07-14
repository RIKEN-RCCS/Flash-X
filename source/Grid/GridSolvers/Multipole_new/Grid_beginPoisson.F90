!!****if* source/Grid/GridSolvers/Multipole_new/Grid_beginPoisson
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
!!  Grid_beginPoisson
!!
!! SYNOPSIS
!!
!!  Grid_beginPoisson (integer(in)    :: iSoln,
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

subroutine Grid_beginPoisson (iSoln,                   &
                              iSrc,                    &
                              bcTypes,                 &
                              bcValues,                &
                              poisFact)

  use Grid_interface,    ONLY : GRID_PDE_BND_ISOLATED

  use Driver_interface,  ONLY : Driver_abort

  use Timers_interface,  ONLY : Timers_start,                  &
                                Timers_stop

  use gr_mpoleInterface, ONLY : gr_mpoleCenterOfExpansion,     &
                                gr_mpoleRadialSampling,        &
                                gr_mpoleSetRadialBinData,      &
                                gr_mpolePrintRadialInfo,       &
                                gr_mpoleMoments,               &
                                gr_mpoleCollectMoments

  use gr_mpoleData,      ONLY : gr_mpoleMomentsDump,           &
                                gr_mpoleMultiThreading,        &
                                gr_mpoleRadialInfoPrint,       &
                                gr_mpoleMomRI,                 &
                                gr_mpoleMomentR,               &
                                gr_mpoleMomentI,               &
                                gr_mpoleScratch,               &
                                gr_mpoleQRadii,                &
                                gr_mpoleQDampingR,             &
                                gr_mpoleQDampingI,             &
                                gr_mpoleQused,                 &
                                gr_mpoleMaxLM,                 &
                                gr_mpoleMaxQ,                  &
                                gr_mpoleAllocateRadialArrays

  implicit none

#include "Simulation.h"
#include "constants.h"

  integer, intent(in)    :: iSoln, iSrc
  integer, intent(in)    :: bcTypes (6)
  real,    intent(in)    :: bcValues (2,6)
  real,    intent(inout) :: poisFact
!  
!
!    ...Check the grid boundaries. Abort, if not isolated
!
!
  if (any (bcTypes (1:2*NDIM) /= GRID_PDE_BND_ISOLATED) ) then
      call Driver_abort ("FATAL: Multipole Poisson solver requires isolated boundaries")
  end if
!  
!
!    ...Start timer.
!
!
  call Timers_start ("Multipole Solver Moments")
!  
!
!    ...Preliminary chores for setting up the run.
!
  !
   call gr_mpoleCenterOfExpansion    (iSrc)
   call gr_mpoleRadialSampling       ()
   call gr_mpoleAllocateRadialArrays (gr_mpoleMomRI,        &
                                       gr_mpoleMomentR,      &
                                       gr_mpoleMomentI,      &
                                       gr_mpoleScratch,      &
                                       gr_mpoleQRadii,       &
                                       gr_mpoleQDampingR,    &
                                       gr_mpoleQDampingI,    &
                                       gr_mpoleQused,        &
                                       gr_mpoleMaxLM,        &
                                       gr_mpoleMaxQ)
   call gr_mpoleSetRadialBinData     ()
!
!
!     ...Print radial info if requested by user.
!
  !
  
  if (gr_mpoleRadialInfoPrint) then
      call gr_mpolePrintRadialInfo ()
  end if
  call Timers_start             ("gr_mpoleMoments")
  call gr_mpoleMoments          (iSrc)
  call Timers_stop              ("gr_mpoleMoments")
  
  call Timers_start             ("gr_mpoleCollectMoments")   ! see "Note" above
  call gr_mpoleCollectMoments ()
  call Timers_stop              ("gr_mpoleCollectMoments")   ! see 'Note' above
  
  call Timers_stop ("Multipole Solver Moments")

  return
end subroutine Grid_beginPoisson
