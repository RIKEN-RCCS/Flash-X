!!****ih* source/Grid/GridMain/UG/Grid_tile
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
!!
!!****

!!REORDER(5): unk, facevar[xyz]
!!REORDER(4): dataPtr

#include "constants.h"
#include "Simulation.h"

module Grid_tile
    implicit none

    private

    type, public :: Grid_tile_t
        integer, public :: id
        integer, public :: cid(MDIM)
        integer, public :: stride(MDIM)
        integer, public :: level
        integer, public :: limits(LOW:HIGH, MDIM)
        integer, public :: grownLimits(LOW:HIGH, MDIM)
        integer, public :: blkLimitsGC(LOW:HIGH, MDIM)
    contains
        procedure, public :: deltas
        procedure, public :: boundBox
        procedure, public :: physicalSize
        procedure, public :: faceBCs
        procedure, public :: getDataPtr
        procedure, public :: releaseDataPtr
        procedure, public :: enclosingBlock
    end type Grid_tile_t

contains

    subroutine deltas(this, dx)
        use Grid_data, ONLY : gr_delta

        class(Grid_tile_t), intent(IN)  :: this
        real,               intent(OUT) :: dx(1:MDIM)

        dx = gr_delta(:, 1)
    end subroutine deltas

    subroutine boundBox(this, box)
        use Grid_data, ONLY : gr_delta
        use Grid_data, ONLY : gr_iCoords, gr_jCoords, gr_kCoords
        use Grid_data, ONLY : gr_ilo, gr_ihi, &
                              gr_jlo, gr_jhi, &
                              gr_klo, gr_khi

        class(Grid_tile_t), intent(IN)  :: this
        real,               intent(OUT) :: box(LOW:HIGH, 1:MDIM)
        integer         :: i

        box(:, :) = 0.0
   
        box(LOW,  IAXIS) = gr_iCoords(LEFT_EDGE,  gr_ilo, 1)
        box(HIGH, IAXIS) = gr_iCoords(RIGHT_EDGE, gr_ihi, 1)
#if NDIM > 1
        box(LOW,  JAXIS) = gr_jCoords(LEFT_EDGE,  gr_jlo, 1)
        box(HIGH, JAXIS) = gr_jCoords(RIGHT_EDGE, gr_jhi, 1)
#endif
#if NDIM > 2
        box(LOW,  KAXIS) = gr_kCoords(LEFT_EDGE,  gr_klo, 1)
        box(HIGH, KAXIS) = gr_kCoords(RIGHT_EDGE, gr_khi, 1)
#endif
        do i = 1,NDIM
           box(LOW,i)  = box(LOW,i)  + (this%limits(LOW,i)  - this%blkLimitsGC(LOW,i)  - NGUARD) * gr_delta(i, 1)
           box(HIGH,i) = box(HIGH,i) + (this%limits(HIGH,i) - this%blkLimitsGC(HIGH,i) + NGUARD) * gr_delta(i, 1)
        end do
    end subroutine boundBox

    subroutine physicalSize(this, tileSize) 
        use Grid_data, ONLY:  gr_delta
        use Grid_data, ONLY : gr_iCoords, gr_jCoords, gr_kCoords
        use Grid_data, ONLY : gr_ilo, gr_ihi, &
                              gr_jlo, gr_jhi, &
                              gr_klo, gr_khi

        class(Grid_tile_t), intent(IN)  :: this
        real,               intent(OUT) :: tileSize(1:MDIM)
        integer         :: i

        tileSize(:) = 0.0
        tileSize(IAXIS) =   gr_iCoords(RIGHT_EDGE, gr_ihi, 1) &
                          - gr_iCoords(LEFT_EDGE,  gr_ilo, 1)
#if NDIM > 1
        tileSize(JAXIS) =   gr_jCoords(RIGHT_EDGE, gr_jhi, 1) &
                          - gr_jCoords(LEFT_EDGE,  gr_jlo, 1)
#endif
#if NDIM > 2
        tileSize(KAXIS) =   gr_kCoords(RIGHT_EDGE, gr_khi, 1) &
                          - gr_kCoords(LEFT_EDGE,  gr_klo, 1)
#endif
        do i = 1,NDIM
           tileSize(i) = tileSize(i) + gr_delta(i, 1) * &
                ( this%limits(HIGH,i) - this%limits(LOW,i) + this%blkLimitsGC(LOW,i) - this%blkLimitsGC(HIGH,i) + 2*NGUARD )
        end do
    end subroutine physicalSize
    
    subroutine faceBCs(this, faces, onBoundary)
        use Grid_data, ONLY : gr_blkBC

        class(Grid_tile_t), intent(IN)            :: this
        integer,            intent(OUT)           :: faces(LOW:HIGH, 1:MDIM)
        integer,            intent(OUT), optional :: onBoundary(LOW:HIGH, 1:MDIM)

        faces = gr_blkBC
        where (faces == PERIODIC)
           faces = NOT_BOUNDARY
        end where 

        if (present(onBoundary)) then     
           onBoundary = gr_blkBC
        end if
    end subroutine faceBCs

    pure function enclosingBlock(this)

        class(Grid_tile_t), intent(IN)  :: this
        type(Grid_tile_t)               :: enclosingBlock

        integer         :: i
        integer,parameter,dimension(MDIM) :: guardv = (/(NGUARD,i=1,NDIM),(0,i=NDIM+1,MDIM)/)

        select type(this)
        type is (Grid_tile_t)
           enclosingBlock = this
        class default
           enclosingBlock%id          = this%id
           enclosingBlock%cid         = this%cid
           enclosingBlock%stride      = this%stride
           enclosingBlock%level       = this%level
           enclosingBlock%limits      = this%limits
           enclosingBlock%blkLimitsGC = this%blkLimitsGC
        end select

        enclosingBlock%limits(LOW, 1:NDIM) = &
             this%blkLimitsGC(LOW,  1:NDIM) + guardv(1:NDIM)
        enclosingBlock%limits(HIGH,1:NDIM) = &
             this%blkLimitsGC(HIGH, 1:NDIM) - guardv(1:NDIM)

        enclosingBlock%grownLimits(:,:) = this%blkLimitsGC(:,:)

    end function enclosingBlock

    ! DEV: If the client code requests a pointer to data that is not 
    ! included in the problem, this routine will return a null pointer
    ! without indicating an error.
    !
    ! This gives the client code the possibility to use either preprocessor
    ! checks to avoid calling this routine needlessly or to do runtime checks
    ! of pointers.
    !
    ! DEV: For now, the localFlag parameter might be useful as we pull in more
    !      units from FLASH4.4.  Try to get rid of it along the way.
    subroutine getDataPtr(this, dataPtr, gridDataStruct, localFlag)
        use physicaldata, ONLY : unk, &
                                 facevarx, facevary, facevarz
        use Grid_data,    ONLY : gr_flxx, gr_flxy, gr_flxz, &
                                 gr_ilo,   gr_jlo,   gr_klo, &
                                 gr_iloGc, gr_jloGc, gr_kloGc, &
                                 scratch_ctr

        class(Grid_tile_t), intent(IN),  target   :: this
        real,                            pointer  :: dataPtr(:, :, :, :)
        integer,            intent(IN)            :: gridDataStruct
        logical,            intent(IN),  optional :: localFlag

        ! Avoid possible memory leaks
        if (associated(dataPtr)) then
            call Driver_abort("[getDataPtr] Given data pointer must be NULL")
        end if

        if (present(localFlag)) then
            call Driver_abort("[getDataPtr] localFlag not implemented")
        end if

        if (this%level /= 1) then
            call Driver_abort("[getDataPtr] Level must be one")
        end if

        select case (gridDataStruct)
        case(CENTER)
           dataPtr(1:, gr_iloGc:, gr_jloGc:, gr_kloGc:) => unk(:,:,:,:,1)
        case(FACEX)
           dataPtr => facevarx(1:,gr_iloGc:,gr_jloGc:,gr_kloGc:,1)
        case(FACEY)
           dataPtr => facevary(1:,gr_iloGc:,gr_jloGc:,gr_kloGc:,1)
        case(FACEZ)
           dataPtr => facevarz(1:,gr_iloGc:,gr_jloGc:,gr_kloGc:,1)
        case(FLUXX)
           dataPtr(1:, gr_ilo:, gr_jlo:, gr_klo:) => gr_flxx(:,:,:,:)
        case(FLUXY)
           dataPtr(1:, gr_ilo:, gr_jlo:, gr_klo:) => gr_flxy(:,:,:,:)
        case(FLUXZ)
           dataPtr(1:, gr_ilo:, gr_jlo:, gr_klo:) => gr_flxz(:,:,:,:) 
       case(SCRATCH_CTR)
          dataPtr(1:, gr_iloGc:, gr_jloGc:, gr_kloGc:) => scratch_ctr(:,:,:,:,this%id)
        case DEFAULT
            call Driver_abort("[getDataPtr] Unknown grid data structure")
        end select
    end subroutine getDataPtr

    subroutine releaseDataPtr(this, dataPtr, gridDataStruct)
        class(Grid_tile_t), intent(IN)            :: this
        real,               intent(OUT), pointer  :: dataPtr(:, :, :, :)
        integer,            intent(IN)            :: gridDataStruct

        nullify(dataPtr)
    end subroutine releaseDataPtr

end module Grid_tile

