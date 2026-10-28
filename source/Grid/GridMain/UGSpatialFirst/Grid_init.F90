!!****if* source/Grid/GridMain/UG/Grid_init
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
!!  Grid_init
!!
!!
!! SYNOPSIS
!!
!!  call Grid_init()
!!
!! DESCRIPTION
!!
!!  Initialize the runtime parameters of the Grid unit.
!!
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!  iguard [INTEGER] 
!!   number of guard cells along IAXIS
!!  jguard [INTEGER] 
!!   number of guard cells along JAXIS
!!  kguard [INTEGER] 
!!   number of guard cells along KAXIS
!!  iGridSize [INTEGER]
!!   the global number of interior grid cells along IAXIS
!!   only needed when operating in NONFIXEDBLOSIZE mode
!!  jGridSize [INTEGER]
!!   the global number of interior grid cells along JAXIS
!!   only needed when operating in NONFIXEDBLOSIZE mode
!!  kGridSize [INTEGER]
!!   the global number of interior grid cells along KAXIS
!!   only needed when operating in NONFIXEDBLOSIZE mode
!!   
!! flux_correct [BOOLEAN] 
!!    turns flux correction on or off 
!!    should always be false in UG since all blocks are on same level
!! compute_grid_size [BOOLEAN]
!!    compute grid size in the case of non-fixed-block size,
!!    non fixed block size mode means block dims are not specified at compile time
!! smalle [REAL]  
!!   Cutoff value for energy (used in some boundary condition handling)
!! smallx [REAL]
!!   Cutoff value for abundances
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_init()

  use physicalData
  use Driver_interface, ONLY : Driver_abort, Driver_getMype, &
    Driver_getNumProcs, Driver_getComm
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt

  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs, Driver_getComm
  use Grid_data
  use gr_specificData, ONLY : gr_globalOffset
  use gr_tilePolicyData,ONLY : gr_initTilePolicy
  use gr_sbInterface, ONLY : gr_sbInit

  implicit none

#include "Simulation.h"
#include "constants.h"


  character(len=MAX_STRING_LENGTH) :: xl_bcString,xr_bcString
  character(len=MAX_STRING_LENGTH) :: yl_bcString,yr_bcString
  character(len=MAX_STRING_LENGTH) :: zl_bcString,zr_bcString
  character(len=MAX_STRING_LENGTH) :: eosModeString, grav_boundary_type


  integer :: i, localnxb, localnyb, localnzb

  !! Get all of the parallel environment information from the driver

  call Driver_getMype(GLOBAL_COMM, gr_globalMe)
  call Driver_getNumProcs(GLOBAL_COMM, gr_globalNumProcs)
  call Driver_getComm(GLOBAL_COMM, gr_globalComm)

  call Driver_getMype(MESH_COMM, gr_meshMe)
  call Driver_getNumProcs(MESH_COMM, gr_meshNumProcs)
  call Driver_getComm(MESH_COMM, gr_meshComm)

  call Driver_getMype(MESH_ACROSS_COMM, gr_meshAcrossMe)
  call Driver_getNumProcs(MESH_ACROSS_COMM, gr_meshAcrossNumProcs)
  call Driver_getComm(MESH_ACROSS_COMM, gr_meshAcrossComm)

  do i = 1, MDIM
     call Driver_getMype(AXIS_COMM, gr_axisMe(i),i)
     call Driver_getNumProcs(AXIS_COMM, gr_axisNumProcs(i),i)
     call Driver_getComm(AXIS_COMM, gr_axisComm(i),i)
  end do

  ! Initialization of gr_geometry etc is done in gr_initGeometry, called below.

  call RuntimeParameters_get("bndPriorityOne",gr_bndOrder(1))
  call RuntimeParameters_get("bndPriorityTwo",gr_bndOrder(2))
  call RuntimeParameters_get("bndPriorityThree",gr_bndOrder(3))

  call RuntimeParameters_get("iGridSize", gr_gIndexSize(IAXIS))
  call RuntimeParameters_get("jGridSize", gr_gIndexSize(JAXIS))
  call RuntimeParameters_get("kGridSize", gr_gIndexSize(KAXIS))   

  call RuntimeParameters_get("compute_grid_size",gr_compute_grid_size)


  !get the boundary conditions stored as strings in the flash.par file
  call RuntimeParameters_get("xl_boundary_type", xl_bcString)
  call RuntimeParameters_get("xr_boundary_type", xr_bcString)
  call RuntimeParameters_get("yl_boundary_type", yl_bcString)
  call RuntimeParameters_get("yr_boundary_type", yr_bcString)
  call RuntimeParameters_get("zl_boundary_type", zl_bcString)
  call RuntimeParameters_get("zr_boundary_type", zr_bcString)

  !map the string boundary conditions to integer constants defined in constants.h
  call RuntimeParameters_mapStrToInt(xl_bcString,gr_domainBC(LOW,IAXIS))
  call RuntimeParameters_mapStrToInt(xr_bcString,gr_domainBC(HIGH,IAXIS))
  call RuntimeParameters_mapStrToInt(yl_bcString,gr_domainBC(LOW,JAXIS))
  call RuntimeParameters_mapStrToInt(yr_bcString,gr_domainBC(HIGH,JAXIS))
  call RuntimeParameters_mapStrToInt(zl_bcString,gr_domainBC(LOW,KAXIS))
  call RuntimeParameters_mapStrToInt(zr_bcString,gr_domainBC(HIGH,KAXIS))

  ! Initialization of gr_globalDomain and of gr_imin,...,gr_kmax (which contain
  ! the same information) is done in gr_initGeometry, called below.

  call RuntimeParameters_get('smalle', gr_smalle)
  call RuntimeParameters_get('smallx', gr_smallx)
  call RuntimeParameters_get('smlrho', gr_smallrho)

! Initialize gr_geometry from the "geometry" runtime parameter,
! and perform some other initializations related to the geometry.
! Determine the geometries of the individual dimensions, and scale
! angle value parameters that are expressed in degrees to radians.
! This call must be made after boundary conditions in gr_domainBC
! have been set based on the corresponding runtime parameters.
! After this call, gr_geometry and gr_{i,j,k}{min,max} will have
! been set correctly.
  call gr_initGeometry()


#ifdef FLASH_EOS
  call RuntimeParameters_get("eosMode", eosModeString)
  call RuntimeParameters_mapStrToInt(eosModeString,gr_eosMode)
  call RuntimeParameters_get("eosModeInit", eosModeString)
  call RuntimeParameters_mapStrToInt(eosModeString,gr_eosModeInit)
#else 
  gr_eosMode = 1
  gr_eosModeInit = 1
#endif
  gr_eosModeNow = gr_eosModeInit ! may change after initialization is done

#ifdef FLASH_PARTICLES
  call RuntimeParameters_get('useParticles',gr_useParticles)
#else
  gr_useParticles=.false.
#endif



  gr_globalNumBlocks = gr_meshNumProcs
  gr_globalOffset = gr_meshMe


  gr_justExchangedGC = .false.

#ifdef FIXEDBLOCKSIZE
  call Driver_abort("Uniform Grid is not supported in fixedblocksize mode")
#else

  if (any((gr_gIndexSize / gr_axisNumProcs) <=  0)) then
     if(gr_meshMe == MASTER_PE) print*,'the local block size is:', gr_gIndexSize / gr_axisNumProcs
     call Driver_abort("[Grid_init] Must set runtime parameters iGridSize, jGridSize, kGridSize so that" //  &
          " a block contains at least one cell.")
  end if

#endif

  call gr_create_surr_blks()


  gr_guard = NGUARD
  gr_guard(JAXIS) = gr_guard(JAXIS)*K2D
  gr_guard(KAXIS) = gr_guard(KAXIS)*K3D

  gr_allPeriodic = .true.
  do i = 1,NDIM
     if(gr_domainBC(LOW,i)/=PERIODIC)gr_allPeriodic=.false.
     if(gr_domainBC(HIGH,i)/=PERIODIC)gr_allPeriodic=.false.
  end do

  !Check if there are gravitational isolated boundary conditions
  !in order to determine which solvers to intialize.
  call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)
  gr_isolatedBoundaries = (grav_boundary_type=="isolated")


  !! Now find the minimum cell size (gr_minCellSize)
  !! - only considering non-angle coordinates

  gr_minCellSizes(IAXIS) = (gr_imax - gr_imin) / gr_gIndexSize(IAXIS)
  gr_minCellSize = gr_minCellSizes(IAXIS)

  if (NDIM >= 2) then
     gr_minCellSizes(JAXIS) = (gr_jmax - gr_jmin) / gr_gIndexSize(JAXIS)
     if (.not.gr_dirIsAngular(JAXIS)) then
        gr_minCellSize = min(gr_minCellSize,gr_minCellSizes(JAXIS))
     end if
  end if

  if (NDIM == 3) then
     gr_minCellSizes(KAXIS) = (gr_kmax - gr_kmin) / gr_gIndexSize(KAXIS)
     if (.not. gr_dirIsAngular(KAXIS)) then
        gr_minCellSize = min(gr_minCellSize,gr_minCellSizes(KAXIS))
     end if
  end if



  call gr_setDataStructInfo()
  gr_offset(:,:)=0
  gr_offset(FACEX_DATATYPE,IAXIS)=1
  gr_offset(FACEY_DATATYPE,JAXIS)=1
  gr_offset(FACEZ_DATATYPE,KAXIS)=1
  call gr_createDomain()
  call gr_ptInit()
  call gr_bcInit()
  call gr_ptMapInit()

  call gr_sbInit()

  gr_region=0.0

  call RuntimeParameters_get("gr_enableTiling", gr_enableTiling)
  call RuntimeParameters_get("gr_useTiling", gr_useTiling)
  if (gr_useTiling .AND. .NOT. gr_enableTiling) then
     if (gr_meshMe == MASTER_PE) then
        write(*,*) "[Grid_init] WARNING: Tiling is disabled by gr_enableTiling=.FALSE."
        write(*,*) "                     Therefore gr_useTiling will not take effect."
     end if
     ! DEV: TODO Should we print this to the log as well?
  end if
  call RuntimeParameters_get("gr_tileSizeX", gr_tileSize(IAXIS))
  call RuntimeParameters_get("gr_tileSizeY", gr_tileSize(JAXIS))
  call RuntimeParameters_get("gr_tileSizeZ", gr_tileSize(KAXIS))

  call gr_initTilePolicy

end subroutine Grid_init

subroutine gr_create_surr_blks
  use Grid_data, ONLY : surr_blks, gr_axisNumProcs, gr_meshMe, gr_domainBC, &
       gr_meshNumProcs, gr_axisMe
  implicit none
  integer :: surr_dist(-1:1,-K2D:K2D,-K3D:K3D)
  integer :: iMe, jMe, kMe, iNeigh, jNeigh, kNeigh
  integer :: i, j, k, n , lx, ly, lz
  logical :: outsideI, outsideJ, outsideK

  lx = gr_axisNumProcs(1)
  if (NDIM < 2) then
     ly = 1
  else
     ly = gr_axisNumProcs(2)
  end if
  if (NDIM < 3) then
     lz = 1
  else
     lz = gr_axisNumProcs(3)
  end if

  do k = -K3D, K3D
     kMe = gr_axisMe(3)
     kNeigh = kMe + k
     outsideK = (kNeigh < 0 .OR. kNeigh .GE. lz)
     if (outsideK .AND. k==(-1)) then
        if (gr_domainBC(LOW,KAXIS) == PERIODIC) then
           kNeigh = kNeigh + lz
        else
           kNeigh = gr_domainBC(LOW,KAXIS)
        end if
     else if (outsideK .AND. k==1) then
        if (gr_domainBC(LOW,KAXIS) == PERIODIC) then
           kNeigh = kNeigh - lz
        else
           kNeigh = gr_domainBC(HIGH,KAXIS)
        end if
     end if
     do j = -K2D, K2D
        jMe = gr_axisMe(2)
        jNeigh = jMe + j
        outsideJ = (jNeigh < 0 .OR. jNeigh .GE. ly)
        if (outsideJ .AND. j==(-1)) then
           if (gr_domainBC(LOW,JAXIS) == PERIODIC) then
              jNeigh = jNeigh + ly
           else
              jNeigh = gr_domainBC(LOW,JAXIS)
           end if
        else if (outsideJ .AND. j==1) then
           if (gr_domainBC(LOW,JAXIS) == PERIODIC) then
              jNeigh = jNeigh - ly
           else
              jNeigh = gr_domainBC(HIGH,JAXIS)
           end if
        end if
        do i = -1, 1
           iMe = gr_axisMe(1)
           iNeigh = iMe + i
           outsideI = (iNeigh < 0 .OR. iNeigh .GE. ly)
           if (outsideI .AND. i==(-1)) then
              if (gr_domainBC(LOW,IAXIS) == PERIODIC) then
                 iNeigh = iNeigh + lx
              else
                 iNeigh = gr_domainBC(LOW,IAXIS)
              end if
           else if (outsideI .AND. i==1) then
              if (gr_domainBC(LOW,IAXIS) == PERIODIC) then
                 iNeigh = iNeigh - lx
              else
                 iNeigh = gr_domainBC(HIGH,IAXIS)
              end if
           end if
           if (iNeigh < 0) then
              n = iNeigh
           else if (jNeigh < 0) then
              n = jNeigh
           else if (kNeigh < 0) then
              n = kNeigh
           else
              n = iNeigh + (jNeigh + kNeigh * ly) * lx
           end if
           surr_dist(i,j,k) = n
        end do
     end do
  end do

  surr_blks(1,:,:,:,1) = 1
  surr_blks(2,:,:,:,1) = surr_dist(:,:,:)
  surr_blks(3,:,:,:,1) = 1

#if NDIM == 1
#define THREE_POW 3
#elif NDIM == 2
#define THREE_POW 9
#else
#define THREE_POW 27
#endif
999     format("Neighbors on PE",I5,", lb=",I2,",surr_blks=", &
             THREE_POW("(",I4,",",I4,",",I2,")"))

  !print 999,gr_meshMe,1,surr_blks(1:3,:,1:1+2*K2D,1:1+2*K3D,1)
end subroutine gr_create_surr_blks
