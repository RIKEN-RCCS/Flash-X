!!****if* source/Grid/localAPI/gr_mpolePotentials
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
!!  gr_mpolePotentials
!!
!! SYNOPSIS
!!
!!  gr_mpolePotentials  (integer, intent(in) :: ipotvar,
!!                       real,    intent(in) :: Poisson_factor,
!!                       real,           intent(inout) :: gr_mpoleMomentR (:,:),
!!                       real,           intent(inout) :: gr_mpoleMomentI (:,:),
!!                       real,           intent(inout) :: gr_mpoleScratch (:,:,:),
!!                       integer,           intent(in)    :: gr_mpoleMaxQ,
!!                       real,              intent(out)   :: gr_mpoleGravityConstant,
!!                       real,              intent(in)    :: gr_mpoleFourPiInv,
!!                       integer,           intent(in)    :: gr_mpoleGeometry,
!!                       integer,           intent(inout) :: gr_mpoleRequest )
!!
!! DESCRIPTION
!!
!!  Computes the potential field using the mass moments already
!!  calculated. On output tha variable indexed by ipotvar contains
!!  the potential. The calculations are entirely local to each
!!  processor, since each processor has a local copy of the moments.
!!
!!  This routine calls the appropriate subroutines according to
!!  the geometry specified.
!!
!! ARGUMENTS
!!
!!  ipotvar        : index to variable containing the potential
!!  Poisson_factor : the name says it all 
!!  gr_mpoleMomentR  : regular moment array (allocated externally)
!!  gr_mpoleMomentI  : irregular moment array (allocated externally)
!!  gr_mpoleScratch  : scratch array (allocated externally)
!!  gr_mpoleMaxQ     : maximum number of radial bins
!!  gr_mpoleGravityConstant : the gravitational constant (output)
!!  gr_mpoleFourPiInv      : inverse of 4*pi (input)
!!  gr_mpoleGeometry       : geometry handle (input)
!!  gr_mpoleRequest        : MPI request handle (input/output)
!!
!!***

subroutine gr_mpolePotentials (ipotvar, Poisson_factor,  &
                               gr_mpoleMomentR,          &
                               gr_mpoleMomentI,          &
                               gr_mpoleScratch,          &
                               gr_mpoleMaxQ,             &
                               gr_mpoleGravityConstant,  &
                               gr_mpoleFourPiInv,        &
                               gr_mpoleGeometry,         &
                               gr_mpoleRequest)

  implicit none
    
  integer, intent (in) :: ipotvar
  real,    intent (in) :: Poisson_factor
  real,           intent (inout) :: gr_mpoleMomentR (1:,0:)
  real,           intent (inout) :: gr_mpoleMomentI (1:,1:)
  real,           intent (inout) :: gr_mpoleScratch (1:,1:,1:)
  integer, intent (in) :: gr_mpoleMaxQ
  real,    intent (out) :: gr_mpoleGravityConstant
  real,    intent (in)  :: gr_mpoleFourPiInv
  integer, intent (in)  :: gr_mpoleGeometry
  integer, intent (inout) :: gr_mpoleRequest

  return
end subroutine gr_mpolePotentials
