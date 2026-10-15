!!****if* source/Grid/GridMain/UG/Grid_data
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
!!  Grid_data
!!
!! SYNOPSIS
!!
!!  use Grid_data
!!
!! DESCRIPTION 
!!  
!!  This includes the global integer identifier for a block, the grid geometry information
!!  
!!  Defining data structures for Uniform Grid
!!
!!***


Module Grid_data
  implicit none

#include "Simulation.h"
#include "constants.h"


  real, save, target, allocatable :: scratch(:,:,:,:,:)
  real, save, target, allocatable :: scratch_ctr(:,:,:,:,:)
  real, save, target, allocatable :: scratch_facevarx(:,:,:,:,:)
  real, save, target, allocatable :: scratch_facevary(:,:,:,:,:)
  real, save, target, allocatable :: scratch_facevarz(:,:,:,:,:)


  integer, parameter :: lnblocks = 1
  integer, parameter :: SAVED_VARS=1
  integer, save :: gr_oneRefLev = 1
  !stores the global number of cells in simulation (w/out guardcells) in each dimension
  integer, save, dimension(MDIM) :: gr_gIndexSize
 
  integer, save, dimension(MDIM) :: gr_guard
  integer, save, dimension(MDIM) :: gr_lIndexSize !stores number of cells in local block (no gcells)

  !stores lower left hand global cell index for each local block in each dim
  integer, save, dimension(MDIM) :: gr_blkCornerID

  integer, save :: gr_meshComm, gr_meshMe, gr_meshNumProcs
  integer, save ::  gr_meshAcrossComm, gr_meshAcrossMe, gr_meshAcrossNumProcs
  integer, save :: gr_globalComm, gr_globalMe, gr_globalNumProcs
  integer, save, dimension(MDIM) :: gr_axisComm, gr_axisMe, gr_axisNumProcs
  integer, save, dimension(NDATATYPES,MDIM) :: gr_offset
  integer, save, dimension(NDATATYPES,MDIM) :: gr_exch
  integer, save :: gr_numDataStruct
  integer, save, dimension(NDATATYPES) :: gr_gridDataStruct,gr_gridDataStructSize

  integer, save, dimension(MDIM) :: gr_bndOrder
  logical, save :: gr_allPeriodic
  integer, save, dimension(1)::gr_blkList
  integer, save, dimension(LOW:HIGH,MDIM) :: gr_domainBC, gr_blkBC
  real ,save :: gr_imin,gr_imax,gr_jmin,gr_jmax,gr_kmin,gr_kmax

  integer ,save :: gr_iGridSize,gr_jGridSize,gr_kGridSize
  character(len=MAX_STRING_LENGTH), save :: gr_str_geometry
  integer, save :: gr_geometry
  integer, save, dimension(MDIM) :: gr_dirGeom
  logical, save, dimension(MDIM) :: gr_dirIsAngular
  logical, save :: gr_geometryOverride
  integer, save :: gr_eosMode
  integer, save :: gr_eosModeInit, gr_eosModeNow
  logical, save :: gr_justExchangedGC, gr_useParticles, gr_isolatedBoundaries

  logical, save :: gr_useEnergyDeposition

  real, save, dimension(MDIM,1) :: gr_delta
  real, save, dimension(LOW:HIGH,MDIM) :: gr_globalDomain

  integer, save :: gr_iguard,gr_jguard,gr_kguard
  integer, save :: gr_iloGc,gr_ihiGc,gr_jloGc,gr_jhiGc,gr_kloGc,gr_khiGc
  integer, save :: gr_ilo,gr_ihi,gr_jlo,gr_jhi,gr_klo,gr_khi
  real, save,target, allocatable, dimension(:,:,:) ::  gr_iCoords,gr_jCoords,gr_kCoords
  integer, save:: gr_lrefineMax=1, gr_maxRefine=1

  logical, save :: gr_compute_grid_size !used in reordered UG
  real, save :: gr_minCellSize
  real, save, dimension(MDIM) :: gr_minCellSizes

  !below values needed to make data structures for IO output
!!  integer,save :: gr_globalOffset !stores beginning blk offset for a proc
  integer,save,target,allocatable, dimension(:,:) :: gr_gid  !holds neigh, child, parent info for checkpoint files
  integer, save :: gr_globalNumBlocks !
  integer, save, allocatable :: gr_nToLeft(:) !array holding local blocks on each proc, used to calc blk offset for IO
  real, save :: gr_smallx
  real, save :: gr_smalle, gr_smallrho
  logical, save :: gr_refineOnParticleCount=.false.

  real, save, dimension(LOW:HIGH,MDIM) :: gr_region
  integer, save :: surr_blks(3,3,1+2*K2D,1+2*K3D,1)

  logical, save :: gr_bcEnableApplyMixedGds

  real,save,target,allocatable,dimension(:,:,:,:) :: gr_flxx, gr_flxy, gr_flxz
  
  logical, save :: gr_enableTiling
  logical, save :: gr_useTiling
  integer, save :: gr_tileSize(1:MDIM)
  integer, save :: gr_maxCells
end Module Grid_data
