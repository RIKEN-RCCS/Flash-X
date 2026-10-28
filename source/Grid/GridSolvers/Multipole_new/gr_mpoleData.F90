!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleData
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
!!  gr_mpoleData
!!
!! SYNOPSIS
!!
!!  use gr_mpoleData
!!
!! DESCRIPTION
!!
!!  Data module for Multipole Solver
!!  --------------------------------
!!   
!!   Legend: (P) means data that is set as (P)arameters
!!           (G) means data that is (G)et from other units (driver, physical constants, etc ...)
!!           (R) means data that is supplied as input through (R)untime parameters
!!           (C) means data that is (C)alculated internally by the multipole solver code
!!
!!   gr_mpoleDomainX/Y/Z/R/Phi/ThetaMax (R,C) : Largest domain boundaries of the problem
!!   gr_mpoleDomainX/Y/Z/R/Phi/ThetaMin (R,C) : Lowest domain boundaries of the problem
!!   gr_mpoleDr                           (C) : Smallest (atomic) radial size for the outer zone
!!   gr_mpoleDrInnerZone                  (C) : Smallest (atomic) radial size for the inner zone
!!   gr_mpoleDrInnerZoneInv               (C) : Inverse of the inner zone smallest radial size
!!   gr_mpoleDrInv                        (C) : Inverse of the outer zone smallest radial size
!!   gr_mpoleEbase                        (C) : Value of 'e'
!!   gr_mpoleEbaseInv                     (C) : Value of '1/e'
!!   gr_mpoleGeometry                     (C) : Internal handle representing the geometry of the problem
!!   gr_mpoleGravityConstant              (C) : Internal copy of the gravitational constant
!!   gr_mpoleIgnoreInnerZone              (R) : If true, only the outer (statistical) zone(s) are considered
!!   gr_mpoleInnerZoneDrRadii             (C) : Maximum radial value (in atomic radial size) of each inner zone bin
!!   gr_mpoleInnerZoneExists              (C) : If true, the inner zone will be calculated
!!   gr_mpoleInnerZoneMaxR                (C) : Maximum radial value for the inner zone
!!   gr_mpoleInnerZoneQmax                (C) : Maximum number of radial bins for the inner zone
!!   gr_mpoleInnerZoneQlower/upper        (C) : Radial bin boundaries for efficient inner zone bin search
!!   gr_mpoleInnerZoneResolution          (R) : Width of radial bins for the inner zone
!!   gr_mpoleInnerZoneResolutionInv       (C) : Inverse of the width of radial bins for the inner zone
!!   gr_mpoleInnerZoneSize                (R) : Size of the inner zone (in units of 'atomic' inner zone length)
!!   gr_mpoleMax2L                        (C) : Twice the largest angular momentum L
!!   gr_mpoleMaxL                         (R) : The largest angular momentum L specified for the run
!!   gr_mpoleMaxLM                        (C) : The overall dimension of the moments
!!   gr_mpoleMaxM                         (C) : The largest magnetic number
!!   gr_mpoleMaxQ                         (C) : Maximum number of radial bins for the entire problem
!!   gr_mpoleMaxR                         (C) : Maximum radial value for the entire problem
!!   gr_mpoleMaxRadialZones               (R) : Maximum number of outer (statistical) radial zones
!!   gr_mpoleMinRadialZone                (C) : Minimum outer (statistical) radial zone index
!!   gr_mpoleMomentI                      (C) : Irregular moment array (2-dimensional)
!!   gr_mpoleMomentR                      (C) : Regular moment array (2-dimensional)
!!   gr_mpoleMomentsDump                  (R) : If true, the moments will be dumped at each time step
!!   gr_mpoleMultiThreading               (R) : If true, the code will run in multithreaded mode
!!   gr_mpoleNumberInv                    (C) : Array containing inverses of whole numbers
!!   gr_mpoleOuterZoneExists              (C) : If true, the outer zone will be calculated
!!   gr_mpoleOuterZoneQshift              (C) : Shift value connecting smoothly the inner/outer zones radial bins
!!   gr_mpole/Two/Half/...Pi              (C) : Several fixed constants related to Pi
!!   gr_mpoleQ                            (C) : Contians a list of radial bin indices
!!   gr_mpoleQdataCells1D/2D/3D           (C) : Data (1D/2D/3D geometries) for each cell for all radial bin indices
!!   gr_mpoleQDampingR/I                  (C) : Radial damping factor parts associated with radial bin value
!!   gr_mpoleQRadii                       (C) : Upper radial limit value for each radial bin
!!   gr_mpoleQused                        (C) : Monitoring array for the use of each radial bin
!!   gr_mpoleQnumberOfCells               (C) : Associates the number of cells with each radial bin
!!   gr_mpoleRad2Deg                      (C) : Angular conversion factor radians -> degrees
!!   gr_mpoleRadialInfoPrint              (R) : If true, detailed radial info will be printed at each time step
!!   gr_mpoleScratch                      (C) : Scratch array (2-dimensional, same size as moment arrays)
!!   gr_mpoleSymmetryAxis3D               (R) : If true, the code assumes a symmetry axis along the z-axis
!!   gr_mpoleSymmetryPlane2D              (R) : If true, the code assumes a symmetry axis along the radial axis
!!   gr_mpoleTotalMass                    (C) : The total mass of the problem
!!   gr_mpoleTotalNrCosineMoments         (C) : The total number of moments associated with a cosine component
!!   gr_mpoleX/Y/Z/R/Phi/ThetaCenter      (C) : Location coordinates of the multipole center of expansion
!!   gr_mpoleZoneExponent                 (R) : Outer zones exponential values
!!   gr_mpoleZoneExponentInv              (C) : Inverses of outer zones exponential values (for speedup)
!!   gr_mpoleZoneLogNorm                  (C) : Internal factors (speedup) for logarithmic outer zones
!!   gr_mpoleZoneLogNormInv               (C) : Inverses of internal factors (for speedup)
!!   gr_mpoleZoneMaxRadiusFraction        (R) : Fraction (maximum) of entire radial domain defining each outer zone
!!   gr_mpoleZoneQmax                     (C) : Maximum radial bin index (global) for each outer radial zone
!!   gr_mpoleZoneRmax                     (C) : Maximum radial value (global) for each outer radial zone
!!   gr_mpoleZoneScalar                   (R) : Outer zones scalar factors (for adjusting radial bin sizes) 
!!   gr_mpoleZoneScalarInv                (R) : Inverses of outer zones scalar factors (for speedup)
!!   gr_mpoleZoneType                     (C) : Type of outer radial zone (exponential/logarithmic)
!!
!!
!!  Explanation of the array types used to accumulate data for evaluation of moments
!!  --------------------------------------------------------------------------------
!!
!!   type cellData1D:
!!
!!     mass           (C) : The mass associated with the cell
!!     radius         (C) : The radial distance of the cell from the center of expansion
!!
!!   type cellData2D:
!!
!!     coord1         (C) : The 1st coordinate of the cell location
!!     mass           (C) : The mass associated with the cell
!!     radius         (C) : The radial distance of the cell from the center of expansion
!!
!!   type cellData3D:
!!
!!     coord1         (C) : The 1st coordinate of the cell location
!!     coord2         (C) : The 2nd coordinate of the cell location
!!     coord3         (C) : The 3rd coordinate of the cell location
!!     mass           (C) : The mass associated with the cell
!!     radius         (C) : The radial distance of the cell from the center of expansion
!!
!!
!!***

#include "constants.h"
#include "Simulation.h"

module gr_mpoleData

  implicit none
  
  logical, save :: gr_mpoleIgnoreInnerZone
  logical, save :: gr_mpoleInnerZoneExists
  logical, save :: gr_mpoleMomentsDump
  logical, save :: gr_mpoleMultiThreading
  logical, save :: gr_mpoleThreadTileList
  logical, save :: gr_mpoleThreadWithinTile
  logical, save :: gr_mpoleOuterZoneExists
  logical, save :: gr_mpoleRadialInfoPrint
  logical, save :: gr_mpoleSymmetryAxis3D
  logical, save :: gr_mpoleSymmetryPlane2D

  integer, save :: gr_mpoleInnerZoneQmax
  integer, save :: gr_mpoleInnerZoneSize
  integer, save :: gr_mpoleGeometry
  integer, save :: gr_mpoleMax2L
  integer, save :: gr_mpoleMaxL
  integer, save :: gr_mpoleMaxLM
  integer, save :: gr_mpoleMaxM
  integer, save :: gr_mpoleMaxQ
  integer, save :: gr_mpoleMaxRadialZones
  integer, save :: gr_mpoleMinRadialZone
  integer, save :: gr_mpoleOuterZoneQshift
  integer, save :: gr_mpoleTotalNrCosineMoments
  integer, save :: gr_mpoleRequest  

  real,    save :: gr_mpoleDr
  real,    save :: gr_mpoleDrInv
  real,    save :: gr_mpoleDrInnerZone
  real,    save :: gr_mpoleDrInnerZoneInv
  real,    save :: gr_mpoleGravityConstant
  real,    save :: gr_mpoleInnerZoneResolution
  real,    save :: gr_mpoleInnerZoneResolutionInv
  real,    save :: gr_mpoleInnerZoneMaxR
  real,    save :: gr_mpoleMaxR
  real,    save :: gr_mpoleTotalMass
  real,    save :: gr_mpoleXcenter
  real,    save :: gr_mpoleYcenter
  real,    save :: gr_mpoleZcenter
  real,    save :: gr_mpoleRcenter
  real,    save :: gr_mpolePhiCenter
  real,    save :: gr_mpoleThetaCenter
  real,    save :: gr_mpoleDomainXmax
  real,    save :: gr_mpoleDomainYmax
  real,    save :: gr_mpoleDomainZmax
  real,    save :: gr_mpoleDomainXmin
  real,    save :: gr_mpoleDomainYmin
  real,    save :: gr_mpoleDomainZmin
  real,    save :: gr_mpoleDomainRmax
  real,    save :: gr_mpoleDomainRmin
  real,    save :: gr_mpoleDomainPhiMax
  real,    save :: gr_mpoleDomainPhiMin
  real,    save :: gr_mpoleDomainThetaMax
  real,    save :: gr_mpoleDomainThetaMin
  real,    save :: gr_mpolePi
  real,    save :: gr_mpoleTwoPi
  real,    save :: gr_mpoleHalfPi
  real,    save :: gr_mpoleThirdPi
  real,    save :: gr_mpoleSixthPi
  real,    save :: gr_mpoleFourPi
  real,    save :: gr_mpoleFourPiInv
  real,    save :: gr_mpoleEbase
  real,    save :: gr_mpoleEbaseInv
  real,    save :: gr_mpoleRad2Deg

  real,    save :: gr_mpoleXdens2CoM
  real,    save :: gr_mpoleYdens2CoM
  real,    save :: gr_mpoleZdens2CoM
  real,    save :: gr_mpoleXcenterOfMass
  real,    save :: gr_mpoleYcenterOfMass
  real,    save :: gr_mpoleZcenterOfMass

  type cellData1D
    real    :: cellMass
    real    :: radius
  end type cellData1D

  type cellData2D
    real    :: coord1
    real    :: cellMass
    real    :: radius
  end type cellData2D

  type cellData3D
    real    :: coord1
    real    :: coord2
    real    :: coord3
    real    :: cellMass
    real    :: radius
  end type cellData3D

  integer,           allocatable, save :: gr_mpoleInnerZoneQlower          (:)
  integer,           allocatable, save :: gr_mpoleInnerZoneQupper          (:)
  integer,           allocatable, save :: gr_mpoleQ                        (:)
  integer,           allocatable, save :: gr_mpoleQnumberOfCells           (:)
  integer,           allocatable, save :: gr_mpoleQused                    (:)
  integer,           allocatable, save :: gr_mpoleZoneQmax                 (:)
  integer,           allocatable, save :: gr_mpoleZoneType                 (:)

  real,              allocatable, save :: gr_mpoleNumberInv                (:)
  real,              allocatable, save :: gr_mpoleQDampingR                (:)
  real,              allocatable, save :: gr_mpoleQDampingI                (:)
  real,              allocatable, save :: gr_mpoleQRadii                   (:)
  real,              allocatable, save :: gr_mpoleZoneScalar               (:)
  real,              allocatable, save :: gr_mpoleZoneExponent             (:)
  real,              allocatable, save :: gr_mpoleZoneLogNorm              (:)
  real,              allocatable, save :: gr_mpoleZoneScalarInv            (:)
  real,              allocatable, save :: gr_mpoleZoneExponentInv          (:)
  real,              allocatable, save :: gr_mpoleZoneLogNormInv           (:)
  real,              allocatable, save :: gr_mpoleZoneMaxRadiusFraction    (:)
  real,              allocatable, save :: gr_mpoleZoneRmax                 (:)
  real,              allocatable, save :: gr_mpoleInnerZoneDrRadii         (:)
  real,              allocatable, save :: gr_mpoleMomRI                    (:,:,:)
  real,              allocatable, save :: gr_mpoleMomentR                  (:,:)
  real,              allocatable, save :: gr_mpoleMomentI                  (:,:)
  real,              allocatable, save :: gr_mpoleScratch                  (:,:,:)

  type (cellData1D), allocatable, save :: gr_mpoleQdataCells1D             (:,:)
  type (cellData2D), allocatable, save :: gr_mpoleQdataCells2D             (:,:)
  type (cellData3D), allocatable, save :: gr_mpoleQdataCells3D             (:,:)

  public :: gr_mpoleAllocateRadialArrays

contains

!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleAllocateRadialArrays
! NAME
!
!  gr_mpoleAllocateRadialArrays
!
! SYNOPSIS
!
!  gr_mpoleAllocateRadialArrays (gr_mpoleMomRI,        &
!                                 gr_mpoleMomentR,      &
!                                 gr_mpoleMomentI,      &
!                                 gr_mpoleScratch,      &
!                                 gr_mpoleQRadii,       &
!                                 gr_mpoleQDampingR,    &
!                                 gr_mpoleQDampingI,    &
!                                 gr_mpoleQused,        &
!                                 gr_mpoleMaxLM,        &
!                                 gr_mpoleMaxQ)
!
! DESCRIPTION
!
!  Allocates the radial multi-bin arrays. Before this routine can be called,
!  the maximum number of radial bins for the current iteration must have
!  been established.
!
!***

subroutine gr_mpoleAllocateRadialArrays (gr_mpoleMomRI,        &
                                         gr_mpoleMomentR,      &
                                         gr_mpoleMomentI,      &
                                         gr_mpoleScratch,      &
                                         gr_mpoleQRadii,       &
                                         gr_mpoleQDampingR,    &
                                         gr_mpoleQDampingI,    &
                                         gr_mpoleQused,        &
                                         gr_mpoleMaxLM,        &
                                         gr_mpoleMaxQ)

  use Driver_interface,  ONLY : Driver_abort

  implicit none

  real,    allocatable, intent (inout) :: gr_mpoleMomRI     (:,:,:)
  real,    allocatable, intent (inout) :: gr_mpoleMomentR   (:,:)
  real,    allocatable, intent (inout) :: gr_mpoleMomentI   (:,:)
  real,    allocatable, intent (inout) :: gr_mpoleScratch   (:,:,:)
  real,    allocatable, intent (inout) :: gr_mpoleQRadii    (:)
  real,    allocatable, intent (inout) :: gr_mpoleQDampingR (:)
  real,    allocatable, intent (inout) :: gr_mpoleQDampingI (:)
  integer, allocatable, intent (inout) :: gr_mpoleQused     (:)

  integer,              intent (in)    :: gr_mpoleMaxLM
  integer,              intent (in)    :: gr_mpoleMaxQ

  integer :: status
!
!
!       ...Allocate the radial multi-bin moment arrays and the corresponding
!          scratch array. Note the extra outer radial bin for the irregular
!          moment array and the extra inner radial bin for the regular moment
!          array!
!
!
  allocate (gr_mpoleMomRI(1:gr_mpoleMaxLM,1:gr_mpoleMaxQ,2))
  allocate (gr_mpoleMomentR (1:gr_mpoleMaxLM,0:gr_mpoleMaxQ  ), stat = status)

  if (status > 0) then
      call Driver_abort ('[gr_mpoleAllocateRadialArrays] ERROR: gr_mpoleMomentR allocate failed')
  end if

  allocate (gr_mpoleMomentI (1:gr_mpoleMaxLM,1:gr_mpoleMaxQ+1), stat = status)

  if (status > 0) then
      call Driver_abort ('[gr_mpoleAllocateRadialArrays] ERROR: gr_mpoleMomentI allocate failed')
  end if
!!$  nullify(gr_mpoleMomentR)
!!$  nullify(gr_mpoleMomentI)
!!$  gr_mpoleMomentR => gr_mpoleMomRI(:,:,1)
!!$  gr_mpoleMomentI => gr_mpoleMomRI(:,:,2)
  allocate (gr_mpoleScratch (1:gr_mpoleMaxLM,1:gr_mpoleMaxQ,2  ), stat = status)

  if (status > 0) then
      call Driver_abort ('[gr_mpoleAllocateRadialArrays] ERROR: gr_mpoleScratch allocate failed')
  end if
!
!
!       ...Allocate the radial multi-bin radii array.
!
!
  allocate (gr_mpoleQRadii (0:gr_mpoleMaxQ), stat = status)

  if (status > 0) then
      call Driver_abort ('[gr_mpoleAllocateRadialArrays] ERROR: gr_mpoleQRadii allocate failed')
  end if
!
!
!       ...Allocate the radial multi-bin damping arrays.
!
!
  allocate (gr_mpoleQDampingR (0:gr_mpoleMaxQ  ), stat = status)

  if (status > 0) then
      call Driver_abort ('[gr_mpoleAllocateRadialArrays] ERROR: gr_mpoleQDampingR allocate failed')
  end if

  allocate (gr_mpoleQDampingI (1:gr_mpoleMaxQ+1), stat = status)

  if (status > 0) then
      call Driver_abort ('[gr_mpoleAllocateRadialArrays] ERROR: gr_mpoleQDampingI allocate failed')
  end if
!
!
!       ...Allocate the radial multi-bin monitoring array.
!
!
  allocate (gr_mpoleQused (1:gr_mpoleMaxQ), stat = status)

  if (status > 0) then
      call Driver_abort ('[gr_mpoleAllocateRadialArrays] ERROR: gr_mpoleQused allocate failed')
  end if
!
!
!       Done.
!
!
  return
end subroutine gr_mpoleAllocateRadialArrays

end module gr_mpoleData
